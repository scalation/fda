
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller, Michael Cotterell, Hao Peng, Dong-Yu Yu
 *  @version 1.2
 *  @date    Thu Sep 22 21:45:58 EDT 2016
 *  @see     LICENSE (MIT style license file).
 *
 *  @see en.wikipedia.org/wiki/B-spline
 *  @see cran.r-project.org/web/packages/crs/vignettes/spline_primer.pdf
 *  @see http://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html
 */

package scalation.analytics.fda

import scalation.linalgebra.VectorD
import scalation.math.double_exp
import scalation.util.Error
import scalation.math.ExtremeD.TOL

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `B_Spline` class provides B-Spline basis functions for various orders 'm',
 *  where the order is one more than the degree.  A spline function is a piecewise
 *  polynomial function where the pieces are stitched together at knots with the
 *  goal of maintaining continuity and differentability.  B-Spline basis functions
 *  form a popular form of basis functions used in Functional Data Analysis.
 *  @see open.uct.ac.za/bitstream/item/16664/thesis_sci_2015_essomba_rene_franck.pdf?sequence=1
 *-----------------------------------------------------------------------------
 *  @param ττ    the time-points of the original knots in the time dimension
 *  @param mMax  the maximum order, allowing splines orders from 1 to mMax
 */
class B_Spline (ττ: VectorD, mMax: Int = 4)
      extends Error
{
    private val DEBUG = true                                 // debug flag
    private val l     = ττ.dim - 1                           // the number of intervals
//    ττ(0)        -= TOL
//    ττ(ττ.dim-1) += TOL
    private val ante  = VectorD.fill (mMax-1)(ττ(0))         // "before" knots
    private val post  = VectorD.fill (mMax-1)(ττ(l))         // "after" knots
    private val τ     = ante ++ ττ ++ post                   // augment to accomodate knots

    if (mMax < 1 || mMax > 10) flaw ("constructor", "B_Spline order restricted to 1 thru 10")

    if (DEBUG) {
        println (s"B_Spline (ττ = $ττ, mMax = $mMax)")
        println ("l = " + l)
        println ("τ = " + τ)
    } // if

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Range of "usable" splines when using the `bs` function.
     *  @param m  the order of the spline
     */
    def range (m: Int) = 0 to l+m-2

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Order one 'm = 1' B-Spline basis functions (flat functions).
     *  @param j  indicates which spline function 0 to l-1
     *  @param t  the time parameter
     */
    def b1 (j: Int, t: Double): Double =
    {
        val k = j + mMax - 1                                    // index shift
        if (τ(k) <= t && t < τ(k+1)) 1.0 else 0.0
    } // b1

    private def bb1 (k: Int, t: Double): Double =
    {
        if (τ(k) <= t && t < τ(k+1)) 1.0 else 0.0
    } // bb_1

/*

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Order two 'm = 2' B-Spline basis functions (linear functions).
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def b2 (j: Int, t: Double): Double =
    {
        if (mMax < 2) flaw ("b2", s"mMax = $mMax can't be less than m = 2")
        val k  = j + mMax - 2                                   // index shift
        val d1 = τ(k+1) - τ(k)
        val d2 = τ(k+2) - τ(k+1)
        val a  = if (d1 =~ 0.0) 0.0 else ((t  -  τ(k)) / d1) * bb1 (k, t)
        val b  = if (d2 =~ 0.0) 0.0 else ((τ(k+2) - t) / d2) * bb1 (k+1, t)
        a + b
    } // b2

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Order three 'm = 3' B-Spline basis functions (quadratic functions).
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def b3 (j: Int, t: Double): Double =
    {
        if (mMax < 3) flaw ("b3", s"mMax = $mMax can't be less than m = 3")
        val k  = j + mMax - 3                                   // index shift
        val d1 = τ(k+2) - τ(k)
        val d2 = τ(k+3) - τ(k+1)
        val a  = if (d1 =~ 0.0) 0.0 else ((t  -  τ(k)) / d1) * bb (2)(k, t)
        val b  = if (d2 =~ 0.0) 0.0 else ((τ(k+3) - t) / d2) * bb (2)(k+1, t)
        a + b
    } // b3

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Order four 'm = 4' B-Spline basis functions (cubic functions).
     *  FIX: missing first interesting spline, see B_SplineTest
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def b4 (j: Int, t: Double): Double =
    {
        if (mMax < 4) flaw ("b4", s"mMax = $mMax can't be less than m = 4")
        val k  = j + mMax - 4                                   // index shift
        val d1 = τ(k+3) - τ(k)
        val d2 = τ(k+4) - τ(k+1)
        val a  = if (d1 =~ 0.0) 0.0 else ((t  -  τ(k)) / d1) * bb (3)(k, t)
        val b  = if (d1 =~ 0.0) 0.0 else ((τ(k+4) - t) / d2) * bb (3)(k+1, t)
        a + b
    } // b4

*/

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Order 'm' B-Spline basis functions (general recurrence).
     *  @param m  the order of the spline function (degree = order - 1)
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def bb (m: Int)(j: Int, t: Double): Double =
    {
        if (mMax < m) flaw ("bb", s"mMax = $mMax can't be less than m = $m")
        val EPS = TOL // 0.001
        val tt = t - EPS
        val k = j
        if (m == 1) return bb1 (k, tt)
        val km = k + m
        val n1 = tt  - τ(k)
        val n2 = τ(km) - tt
        val d1 = τ(km-1) - τ(k) + EPS
        val d2 = τ(km) - τ(k+1) + EPS
        //val a  = if (d1 =~ 0.0) 0.0 else ((tt  - τ(k)) / d1) * bb (m-1)(k, t)
        //val b  = if (d2 =~ 0.0) 0.0 else ((τ(km) - tt) / d2) * bb (m-1)(k+1, t)
        val a  = (n1 / d1) * bb (m-1)(k, t)
        val b  = (n2 / d2) * bb (m-1)(k+1, t)
        a + b
    } // bb

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Adjusted order 'm' B-Spline basis functions (general recurrence). These
     *  are adjusted so that the first "usable" spline is at `j = 0`. The valid
     *  range of usable splines is defined in `range`. 
     *  @param m  the order of the spline function (degree = order - 1)
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def bs (m: Int) (j: Int, t: Double): Double = bb (m)(j+mMax-m, t)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** First derivatives of order 'm' B-Spline basis functions (general 
     *  recurrence).
     *  @param m  the order of the spline function (degree = order - 1)
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def d1bb (m:Int) (j: Int, t: Double): Double =
    {
        if (mMax < m) flaw ("d1bb", s"mMax = $mMax can't be less than m = $m")
        val EPS = 0.001
        val tt  = t - EPS
        val k   = j
        if (m == 1) return 0.0
        val km = k + m
        val n1 = tt  - τ(k)
        val n2 = τ(km) - tt
        val d1 = τ(km-1) - τ(k) + EPS
        val d2 = τ(km) - τ(k+1) + EPS
        val a  = ( bb (m-1)(k, t)   + n1 * d1bb (m-1)(k, t)) / d1
        val b  = (-bb (m-1)(k+1, t) + n2 * d1bb (m-1)(k+1, t)) / d2
        a + b
    } // d1bb

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Adjusted derivatives of order 'm' B-Spline basis functions (general 
     *  recurrence). These are adjusted so that the first "usable" spline is at
     *  `j = 0`. The valid range of usable splines is defined in `range`.
     *  @param m  the order of the spline function (degree = order - 1)
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def d1bs (m: Int) (j: Int, t: Double): Double = d1bb (m)(j+mMax-m, t)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Second derivatives of order 'm' B-Spline basis functions (general 
     *  recurrence).
     *  @param m  the order of the spline function (degree = order - 1)
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def d2bb (m:Int) (j: Int, t: Double): Double ={
        if (mMax < m) flaw ("d2bb", s"mMax = $mMax can't be less than m = $m")
        val EPS = 0.001
        val tt = t - EPS
        val k = j
        if (m == 2) return 0.0
        val km = k + m
        val n1 = tt  - τ(k)
        val n2 = τ(km) - tt
        val d1 = τ(km-1) - τ(k) + EPS
        val d2 = τ(km) - τ(k+1) + EPS
        val a  = ( 2.0 * d1bb (m-1)(k, t)   + n1 * d2bb (m-1)(k, t)) / d1
        val b  = (-2.0 * d1bb (m-1)(k+1, t) + n2 * d2bb (m-1)(k+1, t)) / d2
        a + b
    } // d2bb

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    /** Adjusted seconf derivatives of order 'm' B-Spline basis functions 
     *  (general recurrence). These are adjusted so that the first "usable" 
     *  spline is at `j = 0`. The valid range of usable splines is defined in 
     *  `range`.
     *  @param m  the order of the spline function (degree = order - 1)
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def d2bs (m: Int) (j: Int, t: Double): Double = d2bb (m)(j+mMax-m, t)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Adjusted order 'm' B-Spline basis functions (general recurrence). These
     *  are adjusted so that the first "usable" spline is at `j = 0`. The valid
     *  range of usable splines is defined in `range`. 
     *  @param m  the order of the spline function (degree = order - 1)
     *  @param j  indicates which spline function
     *  @param t  the time parameter
     */
    def apply (m: Int) (j: Int, t: Double): Double = bs (m)(j, t)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Plot the B-spline basis functions and as well as their frist and second
     *  derivatives. 
     *  @param m  the order of the spline function (degree = order - 1)
     *  @param t  the time parameter
     */
    def plot (m: Int = mMax) (t: VectorD)
    {
        import scalation.linalgebra.MatrixD
        import scalation.plot.PlotM
        val n   = t.dim
        val d0y = new MatrixD (τ.dim + mMax, n)                    // matrix to hold initial B-Splines
        val d1y = new MatrixD (τ.dim + mMax, n)
        val d2y = new MatrixD (τ.dim + mMax, n)
        for (i <- 0 until n; j <- range(m)) {
            if (m > 0) d0y(j, i) =   bs (m)(j, t(i))
            if (m > 1) d1y(j, i) = d1bs (m)(j, t(i))
            if (m > 2) d2y(j, i) = d2bs (m)(j, t(i))
        } // for
        if (m > 0) new PlotM (t, d0y, null, "B-Spline order " + m)
        if (m > 1) new PlotM (t, d1y, null, "First Derivative of B-Spline order " + m)
        if (m > 2) new PlotM (t, d2y, null, "Second Derivative of B-Spline order " + m)
    } // plot

} // B_Spline class

/*

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `B_SplineTest` object is used to test the `B_Spline` class.
 *  It tests the B-Spline functions for specific orders.
 *  > run-main scalation.analytics.fda.B_SplineTest
 */
object B_SplineTest extends App
{
    import scalation.plot.Plot

    val mM = 4                                               // maximum order to test
    val τ  = VectorD (0.0, 20.0, 40.0, 60.0, 80.0, 100)      // knot time-points
    val bs = new B_Spline (τ, mM)                            // B-Spline generator
    val n  = 100                                             // number of time-points for plotting
    val t  = VectorD.range (0, n)                            // time-points for plotting

    //-------------------------------------------------------------------------
    // order 1 B-Splines (flat functions, degree 0)

    val y1 = new VectorD (n)
    val y2 = new VectorD (n)
    val k = 1

    for (i <- 0 until n) {                                   // order 1 B-Splines
        y1(i) = bs.b1 (k,   t(i))                                 // first "interesting" B-Spline
        y2(i) = bs.b1 (k+1, t(i))                               // next B-Spline
    } // for
    new Plot (t, y1, y2, "B-Spline order " + 1)

    //-------------------------------------------------------------------------
    // order 2 B-Splines (linear functions, degree 1)

    val y3 = new VectorD (n)
    val y4 = new VectorD (n)

    for (i <- 0 until n) {                                   // order 2 B-Splines
        y3(i) = bs.b2 (k, i)                                 // first "interesting" B-Spline
        y4(i) = bs.b2 (k+1, i)                               // next B-Spline
    } // for
    new Plot (t, y3, y4, "B-Spline order " + 2)

    //-------------------------------------------------------------------------
    // order 3 B-Splines (quadratic functions, degree 2)

    val y5 = new VectorD (n)
    val y6 = new VectorD (n)

    for (i <- 0 until n) {                                   // order 3 B-Splines
        y5(i) = bs.b3 (k, i)                                 // first "interesting" B-Spline
        y6(i) = bs.b3 (k+1, i)                               // next B-Spline
    } // for
    new Plot (t, y5, y6, "B-Spline order " + 3)

    //-------------------------------------------------------------------------
    // order 4 B-Splines (cubic functions, degree 3)

    val y7 = new VectorD (n)
    val y8 = new VectorD (n)

    for (i <- 0 until n) {                                   // order 4 B-Splines
        y7(i) = bs.b4 (k, i)                                 // first "interesting" B-Spline
        y8(i) = bs.b4 (k+1, i)                               // next B-Spline
    } // for
    new Plot (t, y7, y8, "B-Spline order " + 4)

    println("tes"+τ(1))

} // B_SplineTest object

*/

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `B_SplineTest2` object is used to test the `B_Spline` class.
 *  It tests the B-Spline functions using the general recurrence.
 *  > run-main scalation.analytics.fda.B_SplineTest2
 */
object B_SplineTest2 extends App
{
    import scalation.plot.Plot

    val mM = 4                                               // maximum order to test
    val τ  = VectorD (0.0, 20.0, 40.0, 60.0, 80.0, 100.0)    // knot time-points
    val bs = new B_Spline (τ, mM)                            // B-Spline generator
    val n  = 100                                             // number of time-points for plotting
    val t  = VectorD.range (0, n)                            // time-points for plotting

    for (m <- 1 to mM) {

        //---------------------------------------------------------------------
        // order m B-Splines (polynomial functions)

        val y1 = new VectorD (n)
        val y2 = new VectorD (n)

        for (i <- 0 until n) {
            y1(i) = bs.bb (m)(mM-m+1, t(i))                     // first "interesting" B-Spline
            y2(i) = bs.bb (m)(mM-m+2, t(i))                     // next B-Spline
        } // for
        new Plot (t, y1, y2, "B-Spline order " + m)

   } // for

} // B_SplineTest2 object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `B_SplineTest3` object is used to test the `B_Spline` class.
 *  It tests the B-Spline functions using the general recurrence and plots
 *  several basis functions using `PlotM`.
 *  > run-main scalation.analytics.fda.B_SplineTest3
 */
object B_SplineTest3 extends App
{
    import scalation.linalgebra.MatrixD
    import scalation.plot.PlotM

    val mM = 4                                               // maximum order to test
    val τ  = VectorD (0.0, 25.0, 50.0, 75.0, 100.0)          // knot time-points
    val bs = new B_Spline (τ, mM)                            // B-Spline generator
    val n  = 100                                             // number of time-points for plotting
    val t  = VectorD.range (0, n)                            // time-points for plotting
    val k  = 0                                               // index for initial B-Spline

    for (m <- 1 to mM) {
        //---------------------------------------------------------------------
        // order m B-Splines (polynomial functions)
        val y = new MatrixD (τ.dim + mM, n)                    // matrix to hold initial B-Splines
        for (i <- 0 until n; j <- 0 to bs.range(m).last) {
            y(j, i) = bs (m)(j, t(i))
        } // for
        new PlotM (t, y, null, "B-Spline order " + m)
   } // for

} // B_SplineTest3 object

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `B_SplineTest4` object is used to test the `B_Spline` class.
  *  It tests the B-Spline functions using the general recurrence.
  *  > run-main scalation.analytics.fda.B_SplineTest4
  */
object B_SplineTest4 extends App
{
    import scalation.linalgebra.MatrixD
    import scalation.plot.PlotM

    val mM = 4                                               // maximum order to test
    val n  = 100
    val τ  = VectorD (0.0, 0.5 * n , n) / (0.5 * n)          // knot time-points
    val bs = new B_Spline (τ, mM)                            // B-Spline generator
    val t  = VectorD.range (0, n) / (0.5 * n)                // time-points for plotting
    val k  = 0                                               // index for initial B-Spline

    for (m <- 2 to mM) {
        //---------------------------------------------------------------------
        // order m B-Splines (polynomial functions)
        val d0y = new MatrixD (τ.dim + mM, n)                    // matrix to hold initial B-Splines
        val d1y = new MatrixD (τ.dim + mM, n)
        val d2y = new MatrixD (τ.dim + mM, n)
        for (i <- 0 until n; j <- 0 to bs.range(m).last) {
            d0y(j, i) = bs (m)(j, t(i))
            d1y(j, i) = bs.d1bs (m)(j, t(i))
            d2y(j, i) = bs.d2bs (m)(j, t(i))
        } // for
        new PlotM (t, d0y, null, "B-Spline order " + m)
        new PlotM (t, d1y, null, "First Derivative of B-Spline order " + m)
        new PlotM (t, d2y, null, "Second Derivative of B-Spline order " + m)
    } // for

} // B_SplineTest4 object

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `B_SplineTest5` object is used to test the `B_Spline` class.
  *  It tests the B-Spline functions using the general recurrence.
  *  > run-main scalation.analytics.fda.B_SplineTest5
  */
object B_SplineTest5 extends App
{
    import scalation.linalgebra.MatrixD
    import scalation.plot.PlotM

    val mM = 4                                               // maximum order to test
    val n  = 100
    val τ  = VectorD (0.0, 0.5 * n , n) / (0.5 * n)          // knot time-points
    val bs = new B_Spline (τ, mM)                            // B-Spline generator
    val t  = VectorD.range (0, n) / (0.5 * n)                // time-points for plotting
    for (m <- 2 to mM) bs.plot(m)(t)

} // B_SplineTest5 object

