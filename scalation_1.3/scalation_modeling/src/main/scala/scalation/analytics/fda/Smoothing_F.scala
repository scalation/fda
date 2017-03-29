
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
    private val DEBUG = false                   // debug flag
    private val GAP   = 2                      // gap between time points and knots
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

    def summary
    {
        println (s"size of λΣ  = ${λopt * Σ.norm}")
        println (s"size of Φ'Φ = ${(Φ.t * Φ).norm}")
        println (s"ratio       = ${(λopt * Σ.norm) / (Φ.t * Φ).norm}")
    } // summary

    def getLambda = λopt

    if (y.dim != m) flaw ("constructor", "require # data points == # data time points")
    if (n > m)      flaw ("constructor", "require # knot points <= # data time points")

    private def makeKnots: VectorD = VectorD.range (0, m/GAP) / ((m-1)/GAP) * t(t.dim-1)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the partial fit for the model using one of several training 
     *  methods. If you multitply the training data by this `MatrixD`, then you
     *  get the symmetric "hat" matrix. If you multitply this `MatrixD` by the
     *  response vector, then you will get the vector estimated model
     *  coefficients.
     *  @param λ  the regularization parameter
     */
    private def pfit (λ: Double): MatrixD = method match {
        case ROUGHNESS => ((Φ.t * W * Φ) + (Σ * λ)).inverse * Φ.t * W
        case RIDGE     => ((Φ.t * W * Φ) + (I * λ)).inverse * Φ.t * W
        case WLS       => (Φ.t * W * Φ).inverse * Φ.t * W
        case OLS       => (Φ.t * Φ).inverse * Φ.t
    } // pfit

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Computes the "hat" matrix.
     *  @param λ  the regularization parameter
     */
    private def H (λ: Double) : MatrixD = Φ * pfit (λ)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the model, i.e., determine the optimal coeifficient 'c' for the
     *  basis functions by finding optimal Lamdba to minimize gcv.
     */
    def train (): VectoD =
    {
        if (λopt == -1) useGCV ()
        c   = pfit (λopt) * y      //
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
        val obj = (m / (m - df (λ))) * (sse / (m - df (λ)))
        //println (s"λ = $λ; ns = $ns; m = $m; df = ${df(λ)}; sse = $sse; obj = $obj")
        // (ns / (ns - df (λ))) * (sse / (ns - df (λ)))
        // (m / (m - df (λ))) * (sse / (m - df (λ)))
        obj 
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
            c   = pfit (l) * y
            e   = y - Φ * c
            sse = e dot e
            gcv (l)
        } // f
        val gs   = new GoldenSectionLS (f)
        val step = 10
        λopt     = gs.search (step)
        //if (DEBUG) println (s"λopt = $λopt")
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

    for (ord <- 4 to 4) {
        val τ   = null                                               // let `Smoothing_F` nake the knots
        val moo = new Smoothing_F (y, t, τ, ord, lambda = -1)                    // smoother
        val c   = moo.train ()                                       // train -> set coefficients
        println (s"y = $y \nt = $t \nc = $c")
        val x = moo.predict (t)                                      // predict for all time points
        new Plot (t, y, x, s"B-Spline Fit: ord = $ord; λ = ${moo.getLambda}")
        moo.summary
    } // for

} // Smoothing_FTest object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest2` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest2
 */
object Smoothing_FTest2 extends App
{
    import scalation.analytics.clusterer.{GapStatistic, KMeansPPClusterer}

    val ord  = 4                            // b-spline order
    val kMax = 30                           // maximum number of clusters

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* LOAD DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    println ("Loading observed data...")
    val data = MatrixD ("../data/gene_expression.csv")  // observed data

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PREPARE DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    println ("Preparing data...")
    val x = data //.slice (0, 1000)
    val z = new MatrixD (x.dim1, x.dim2)    // matrix for smoothed sample
    val t = VectorD.range (0, x.dim2)       // vector for time points
    val c = new MatrixD (x.dim1, x.dim2)    // matrix for smoothing spline coefficients

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM SMOOTHING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    println ("Smoothing...")
    for (i <- x.range1) {
        val y   = x(i)
        val moo = new Smoothing_F (y, t, null, ord, lambda = -1) // smoother
        c(i)    = moo.train ()                                   // train -> set coefficients
        z(i)    = moo.predict (t)                                // predict for all time points
        //new Plot (t, y, z, s"Smoothing_F ord = $ord; λ = ${moo.getLambda}", lines = true)
    } // for

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM KMEANS++ CLUSTERING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    val sseObs = new VectorD (kMax)
    val sseSmo = new VectorD (kMax)
    val kVals  = VectorD.range (1, kMax + 1)

    println ("Clustering... ")
    for (k <- 1 to kMax) {
        print (s"k = $k; ")
        //val (cl1, cls1) = KMeansPPClusterer (x, k)
        val cl1  = new KMeansPPClusterer (x, k)
        val cls1 = cl1.cluster ()
        val sse1 = cl1.sse ()
        sseObs (k-1) = sse1
        print (s"observed SSE = $sse1; ")
        //val (cl2, cls2) = KMeansPPClusterer (z, k)
        val cl2  = new KMeansPPClusterer (z, k)
        val cls2 = cl2.cluster ()
        val sse2 = cl2.sse ()
        sseSmo (k-1) = sse2        
        print (s"smoothed SSE = $sse2")
        println ()
    } // for

    new Plot (kVals, sseObs, sseSmo, "k-means++ SSEs: Observed vs. Smoothed", true)
    new Plot (kVals, sseObs.map(math.log _), sseSmo.map(math.log _), "k-means++ log SSEs: Observed vs. Smoothed", true)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM TIGHT CLUSTERING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    // TODO implement

} // Smoothing_FTest2
