
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.2
 *  @date    Tue Oct 11 16:12:54 EDT 2016
 *  @see     LICENSE (MIT style license file).
 *
 *  @see orfe.princeton.edu/~jqfan/papers/07/WuFanMueller1.pdf
 *  @see www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&cad=rja&uact=8&ved=0ahUKEwizk5f2q9PPAhWCSCYKHVF2Be8QFggxMAM&url=http%3A%2F%2Fanson.ucdavis.edu%2F~mueller%2Fhandbook-final1.pdf&usg=AFQjCNHS96onDE2qFFynU1L1xAx27wh0lA&sig2=PLLilZsXqGaI-GV8g5njQA
 *  @see Functional Data Analysis, Second Edition, Chapter 12
 *  @see http://link.springer.com/book/10.1007%2Fb98888
 */

package scalation.analytics.fda

import scalation.analytics.RidgeRegression
import scalation.linalgebra.{MatrixD, VectoD, VectorD}
import scalation.plot.Plot

import scalation.calculus.Integral.∫
import scalation.calculus.functionS2S2Hilbert
//import scalation.calculus.Hilbert.*

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Regression_F2S` class performs functional linear regression with
 *  scaler response and functional covariates.
 *  <p>
 *      y = a + <b(t),x(t)> + ε
 *  <p>
 *  where `<b, x>` denotes the inner product of `b` and `x`. 
 *  @param y    the response vector
 *  @param x    the covariate matrix - treated as functional
 *  @param t    the time vector
 *  @param τ    the knot vector
 *  @param ord  the order (degree+1) of the B-Splines (2 to 6)
 */
class Regression_F2S (y: VectorD, x: MatrixD, t: VectorD, τ: VectorD, ord: Int = 4)
{
    private val DEBUG = true                             // debug flag

    private var b: VectoD = null                         // regression coefficients
    private var c: VectoD = null                         // smoothing coefficients

    private val bs  = new B_Spline (τ)                   // use B-Spline basis functions
    bs.plot()(t)


    private val xmoo = for (i <- y.range) yield new Smoothing_F (x(i), t, τ, ord)
    private val z    = new MatrixD (y.dim, x.dim2 + 1)    // 2-column data matrix [1, xs]

    z.setCol (0, VectorD.one(y.dim))                              // column of all ones
/**
    for (1 until z.dim2) z.setCol
    xx.setCol (1, sxmooth (x))                            // column of sxmoothed x, i.e., xs
  */

    if (DEBUG) println ("data matrix z = " + z)

    def Φ (k: Int) (tt: Double): Double =
    {
        bs.bs(ord) (k, tt)
    } // Φ

/*
    def xi (i: Int, tt: Double): Double =
    {
        xmoo(i)
    } // xi
 */

    def x_i (i: Int): VectorD =
    {
        def xi (tt: Double) = xmoo(i).predict (tt)
        val (a, b) = (t(0), t(t.dim-1))
        VectorD(for (k <- bs.range(ord)) yield ∫ ((a, b), Φ(k) _ * xi _))
    } // x_i

/*    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the model using the sxmoothed data to find the regression coefficients 'b'.
     */
    def train (): VectoD =
    {
/*
        b = (xx.t * xx).inverse * xx.t * y               // direct solution, may produce NaN
*/
        val lambda = 0.01
        val rrg = new RidgeRegression (xx, y, lambda)
        rrg.train ()
        b = rrg.coefficient
        b
    } // train

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Predict the y-value at time point 'tt'.
     *  @param tt  the given time point
     */
    def predict (tt: Double): Double =
    {
        var sum = 0.0
p        println ("b = " + b)
        println ("c = " + c)
        val xt = VectorD (1.0, xmoo.predict (tt))
        for (j <- xt.indices) sum += b(j) * xt(j)        // c(j) * bs.b3 (j, tt)
        sum
    } // predict

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Smooth the data vector 'x' using B-Spline expansion
     *  @param x  the data vector to sxmooth
     */
    private def sxmooth (x: VectorD): VectorD =
    {
        val xs = new VectorD (x.dim)                     // sxmoothed version of x
        c      = xmoo.train ()                            // sxmoothing coefficients
        if (DEBUG) println ("c = " + c)
        for (j <- t.range) xs(j) = xmoo.predict (t(j))
        xs
    } // smooth
 */

} // Regression_F2S class

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Regression_F2STest` object is used to test the `Regression_F2S` class.
 *  > run-main scalation.analytics.fda.Regression_F2STest
 */
object Regression_F2STest extends App
{
    import scalation.random.Normal

    val normal = Normal (0.0, 0.2)
    val n  = 100 
    val t  = VectorD.range (0, n)                            // time-points for plotting
    val y = t.map ((x: Double) => 3.0 + 5.0 * x * x + normal.gen)
    val x = new MatrixD (y.dim, t.dim)

    println (s"y = $y \nx = $x \n t = $t")

    val ord = 4
    val τ  = VectorD (0.0, 20.0, 40.0, 60.0, 80.0, 100.0)    // knot time-points    
    val rgf = new Regression_F2S (y, x, t, τ, ord)
    rgf.x_i(1)

} // Regression_F2STest object

/*

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Regression_F2STest2` object is used to test the `Regression_F2S` class.
 *  > run-main scalation.analytics.fda.Regression_F2STest2
 */
object Regression_F2STest2 extends App
{
    val ord = 4
    val y   = VectorD ( 2.1,  4.3,  5.9,  7.7, 10.3, 11.8, 14.1, 15.9, 18.1, 20.0, 
                       22.1, 24.3, 25.9, 27.7, 30.3, 31.8, 34.1, 35.9, 38.1, 40.0) 
    val x   = VectorD.range (1, 21)
    val t   = VectorD.range (0, 20) / 20.0
    val τ   = VectorD.range (0, 10 + ord) / 10.0

    println (s"y = $y \nx = $x \n t = $t")

    val rgf = new Regression_F2S (y, x, t, τ)
    println ("b = " + rgf.train ())
    val yp = t.map (rgf.predict (_))
    new Plot (t, y, yp, "Regression - ord = default")

} // Regression_F2STest2 object

*/
