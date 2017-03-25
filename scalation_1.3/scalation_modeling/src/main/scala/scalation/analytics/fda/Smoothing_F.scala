
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell, John Miller, Dong-Yu Yu
 *  @version 1.3
 *  @date    Fri Mar 24 20:24:40 2017
 *  @see     LICENSE (MIT style license file).
 *
 *  @see open.uct.ac.za/bitstream/item/16664/thesis_sci_2015_essomba_rene_franck.pdf?sequence=1
 *  @see www.jstatsoft.org/article/view/v051i04/v51i04.pdf
 *  @see Functional Data Analysis, Second Edition, Chapter 4
 *  @see http://link.springer.com/book/10.1007%2Fb98888
 */

package scalation.analytics.fda

import scalation.analytics.RidgeRegression
import scalation.linalgebra.{MatrixD, VectoD, VectorD}
import scalation.plot.{FPlot, Plot}
import scalation.util.Error

object SmoothingMethod extends Enumeration
{
    type SmoothingMethod = Value
    val ROUGHNESS, RIDGE, WLS, OLS = Value
} // SmoothingMethod

import SmoothingMethod._

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_F` class fits a time-dependent data vector 'y' to B-Splines.
 *  <p>
 *      y(t(i)) = x(t(i)) + ε(t(i))
 *      x(t) = cΦ(t)
 *  <p>
 *  where 'x' is the signal, 'ε' is the noise, 'c' is a coefficient vector and
 *  'Φ(t)' is a vector of basis functions. 
 *-----------------------------------------------------------------------------
 *  @param y       the (raw) data points/vector
 *  @param t       the data time points/vector
 *  @param τ       the time points/vector for the knots
 *  @param ord     the order (degree+1) of B-Splines (2, 3, 4, 5 or 6)
 *  @param method  the smoothing method
 *  @param lambda  the regularization parameter (>= 0 or -1 to use GCV) 
 */
class Smoothing_F (y: VectorD, t: VectorD, private var τ: VectorD = null, ord: Int = 4, method: SmoothingMethod = RIDGE, lambda: Double = 1E-5)
      extends Error
{
    private val DEBUG = true                   // debug flag
    private val GAP   = 5                      // gap between time points and knots
    private val m     = t.dim                  // number of data time points
    if (τ == null) τ  = makeKnots
    private val n     = τ.dim                  // number of time points for the knots
    private val bs    = new B_Spline (τ, ord)  // use B-Spline basis functions
    private val ns    = bs.size ()
    private val Φ     = bs.phi ()(t)           // matrix: jth spline at time ti
    private val Σ     = bs.penalty ()(t)       // penalty matrix
    private val I     = MatrixD.eye(m)         // identity matrix
    private val W     = MatrixD.eye(m)         // weight matrix

    private var λopt  = lambda                 // regularization parameter    
    private var c: VectoD = null               // coefficient vector
    private var e: VectoD = null               // residual/error vector    
    private var sse   = 0.0                    // sum of squared error

    if (y.dim != m) flaw ("constructor", "require # data points == # data time points")
    if (n > m)      flaw ("constructor", "require # knot points <= # data time points")

    private def makeKnots: VectorD = VectorD.range (0, m/GAP) / ((m-1)/GAP) * t(t.dim-1)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the estimated coefficients for the model using one of several
     *  training methods.
     *  @param λ  the regularization parameter
     */
    private def findc (λ: Double): MatrixD = method match {
        case ROUGHNESS => ((Φ.t * W * Φ) + (Σ * λ)).inverse * Φ.t * W
        case RIDGE     => ((Φ.t * W * Φ) + (I * λ)).inverse * Φ.t * W
        case WLS       => (Φ.t * W * Φ).inverse * Φ.t * W
        case OLS       => (Φ.t * Φ).inverse * Φ.t
    } // findc

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Computes the "hat" matrix.
     *  @param λ  the regularization parameter
     */
    private def H (λ: Double) : MatrixD = Φ * findc (λ)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the model, i.e., determine the optimal coeifficient 'c' for the
     *  basis functions by finding optimal Lamdba to minimize gcv.
     */
    def train (): VectoD =
    {
        if (λopt == -1) useGCV ()
        c   = findc (λopt) * y      //
        e   = y - Φ * c
        sse = e dot e
        c
    } // train

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the degrees of freedom based on the trace of hat matrix.
     *  @param λ  the regularization parameter
     */
    private def df (λ: Double): Double = H (λ).trace

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the Generalized Cross Validation (GCV) score.
     *  @param λ  the regularization parameter
     */
    private def gcv (λ: Double): Double =
    {
        (ns / (ns - df (λ))) * (sse / (ns - df (λ)))
    } // gcv

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Use the Generalized Cross Validation (GCV) method to find the optimal
     *  λ value. It finds the λ that minimizes the GCV score. 
     */
    private def useGCV ()
    {
        import scalation.minima.GoldenSectionLS
        def f (l: Double): Double =
        {
            c   = findc (l) * y
            e   = y - Φ * c
            sse = e dot e
            gcv (l)
        } // f
        val gs   = new GoldenSectionLS (f)
        val step = 100
        λopt     = gs.search (step)
        if (DEBUG) println (s"λopt = $λopt")
    } // useGCV

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Predict the y-value at time point 'tt'.
     *  @param tt  the given time point
     */
    def predict (tt: Double): Double =
    {
        var sum = 0.0
        for (j <- bs.range ()) sum += c(j) * bs (ord) (j, tt)
        sum
    } // predict

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Predict the y-values at all time points in vector 'tv'.
     *  @param tv  the given vector of time points
     */
    def predict (tv: VectorD): VectorD = tv.map (predict (_))

} // Smoothing_F class

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest
 */
object Smoothing_FTest extends App
{
    import scalation.random.Normal

    val normal = Normal ()                                           // normal random variate generator
    val t = VectorD.range (0, 100) / 100.0                           // time points
    val y = t.map ((x: Double) => 3.0 + 2.0 * x * x + normal.gen)    // raw data points

    for (ord <- 2 to 10) {
        val τ   = null                                               // let `Smoothing_F` nake the knots
        val moo = new Smoothing_F (y, t, τ, ord)                    // smoother
        val c   = moo.train ()                                       // train -> set coefficients
        println (s"y = $y \nt = $t \nc = $c")
        val x = moo.predict (t)                                      // predict for all time points
        new Plot (t, y, x, s"B-Spline Fit: ord = $ord")
    } // for

} // Smoothing_FTest object


