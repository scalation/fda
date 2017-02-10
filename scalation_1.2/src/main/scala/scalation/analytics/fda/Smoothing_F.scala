
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.2
 *  @date    Thu Sep 22 21:45:58 EDT 2016
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
import scala.util
import java.io.BufferedReader
import java.io.FileReader
import java.io.IOException

import scalation.random.Uniform
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_F` class fits a time-dependent data vector 'y' to B-Splines.
 *  <p>
 *      y(t(i)) = x(t(i)) + ε(t(i))
 *      x(t) = cΦ(t)
 *  <p>
 *  where 'x' is the signal, 'ε' is the noise, 'c' is a coefficient vector and
 *  'Φ(t)' is a vector of basis functions. 
 *-----------------------------------------------------------------------------
 *  @param y    the (raw) data points/vector
 *  @param t    the data time points/vector
 *  @param τ the time points/vector for the knots
 *  @param ord  the order (degree+1) of B-Splines (2, 3, 4, 5 or 6)
 */
class Smoothing_F (y: VectorD, t: VectorD, private var τ: VectorD = null, ord: Int = 4)
      extends Error
{
    private val DEBUG = true                   // debug flag
    private val GAP   = 5                      // gap between time points and knots
    private val m     = t.dim                  // number of data time points
    if (τ == null) τ  = makeKnots
    println("makeKnots="+τ)
    private val n     = τ.dim                  // number of time points for the knots

    if (y.dim != m) flaw ("constructor", "require # data points == # data time points")
    if (n > m)      flaw ("constructor", "require # knot points <= # data time points")


    private var c: VectoD = null                // coefficient vector
    private val bs = new B_Spline (τ, ord)     // use B-Spline basis functions
    private val phi = new MatrixD (m, bs.range(ord).length)       // form a Matrix to record the value of jth spline at time ti

    if (DEBUG) println (s"m = $m, n = $n")

   // def makeKnots: VectorD = VectorD.range (0, (m/GAP + ord)) / (m/GAP.toDouble)
   def makeKnots: VectorD =
   {
       val knots = (VectorD.range (0, 20) / (20))*t(t.dim-1)
       knots(knots.dim-1) = t(t.dim-1)
       knots
   }
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the model, i.e., determine the optimal coeifficient 'c' for the
     *  basis functions by finding optimal Lamdba to minimize gcv.
     */
    def train (): VectoD =
    {
        form_phi ()
        //new Plot(t,phi.col(bs.range(ord).last-1))
        import scalation.minima.GoldenSectionLS                         //Using GoldenSectionSearch to find minimal GCV in a range of Lambda

        /** f is the function gcv = f (Lambda) to put as a parameter of GoldenSectionLS
          */
        def f (l: Double): Double =
        {
            val rrg = new RidgeRegression (phi, y, l)
            rrg.train()
            gcv(l,rrg)
        }
        val gs = new GoldenSectionLS(f)
        val step = 1                                                     // step is the search space (between 0 to step) for f
        val optimalLambda = gs.search(step)                              // find the optimal Lambda to form the minimal gcv
        val minRrg = new RidgeRegression (phi, y, optimalLambda)
        minRrg.train()
// un-note to print out the value of mingcv and optimal Lambda
//        val minGCV = gcv(optimalLambda, minRrg)
//        println (s"mingcv($optimalLambda)="+minGCV)
        c = minRrg.coefficient
        c
    } // train
    def printτ()={

        println(τ.toString)


    }
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count degree of freedom of ridge regression by the trace of hat matrix H(Lambda)
      * where H(Lambda)= X*inverse(transpose(X)*X+Lambda*I)*transpose(X), where X=the matrix phi
      */
    def df (lam: Double):Double =
    {
        import scalation.linalgebra.MatrixD.eye
        (phi*(phi.t*phi+eye(phi.dim2)*lam).inverse*phi.t).trace
    } // df

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /**Count gcv. GCV is a model for reducing under smoothing the data by CV method.
      The formula is (n/(n-df(Lambda)))*(SSE/(n-df(Lambda)))
      */
    def gcv(lam: Double, rrg: RidgeRegression): Double =
    {
        val n = phi.dim1
        (n/(n-df(lam))) *(rrg.sse/(n-df(lam)))
    } // gcv

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Predict the y-value at time point 'tt'.
     *  @param tt  the given time point
     */
    def predict (tt: Double): Double =
    {
        var sum = 0.0
        for (j <- bs.range(ord)) {
            sum += c(j) * bs.bs (ord) (j, tt)
          //  println(s"c($j) = ${c(j)}; bs.bs ($ord) ($j, $tt) = ${bs.bs (ord) (j, tt)}")
        }
        //for (j <- 0 until n) sum += c(j) * bs.b3  (j, tt)
        sum
    } // predict

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Predict the y-values at all time points in vector 'tv'.
     *  @param tv  the given vector of time points
     */
    def predict (tv: VectorD): VectorD = tv.map (predict (_))

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Form the phi matrix by evaluating the basis functions at the time points.
     */
    private def form_phi ()
    {   println("t.indices="+t.indices)
        for (i <- t.indices; j <- bs.range(ord)) phi (i, j) = bs.bs(ord) (j, t(i))
        //for (i <- t.indices; j <- 0 until n) phi (i, j) = bs.bb (ord) (j, t(i))
    } // form_phi

    def print_phi(): Unit =
    {
        println(phi.toString)
    }

} // Smoothing_F class

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest
 */
object Smoothing_FTest extends App
{
    import scalation.random.Normal

    val normal = Normal ()                                           // normal random variate generator
    val t = VectorD.range (0, 101) / 100.0                           // time points
    val y = t.map ((x: Double) => 3.0 + 2.0 * x * x + normal.gen)    // raw data points

    for (ord <- 2 to 10) {
//      val τ   = VectorD.range (0, 20 + ord) / 20.0                 // time points for knots
        val τ   = null                                               // let `Smoothing_F` nake the knots
        val moo = new Smoothing_F (y, t, τ, ord)                    // smoother
        val c   = moo.train ()                                       // train -> set coefficients

        println (s"y = $y \nt = $t \nc = $c")

        val x = moo.predict (t)// predict for all time points

        println("last ="+moo.predict(t(t.dim-1)))
        new Plot (t, y, x, s"B-Spline Fit: ord = $ord")
    } // for

} // Smoothing_FTest object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest2` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest2
 */
object Smoothing_FTest2 extends App
{
    import scalation.random.Normal

    val normal = Normal ()
    val t  = VectorD.range (0, 100) / 100.0
    val t2 = VectorD.range (0, 1000) / 1000.0
    val y  = t.map ((x: Double) => 3.0 + 2.0 * x * x + normal.gen)

    for (ord <- 2 to 6) {
        val τ   = VectorD.range (0, 20 + ord) / 20.0
        val moo = new Smoothing_F (y, t, τ, ord)
        val c   = moo.train ()

        println (s"y = $y \nt = $t \nc = $c")
       // moo.print_phi()
        new FPlot (t, y, t2, moo.predict, s"B-Spline Fit: ord = $ord")
    } // for

} // Smoothing_FTest2 object

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest3` is used to test the `Smoothing_F` class.
  *  > run-main scalation.analytics.fda.Smoothing_FTest3
  */
object Smoothing_FTest3 extends App
{
    import scalation.random.Normal

    val normal = Normal ()
    val t  = VectorD.range (0, 100) / 100.0
    val t2 = VectorD.range (0, 1000) / 1000.0
    val y  = t.map ((x: Double) => 3.0 + 2.0 * x * x + normal.gen)

    for (ord <- 3 to 3) {
        val τ   = VectorD.range (0, 20 + ord) / 20.0
        val moo = new Smoothing_F (y, t, τ, ord)
      //
        print(moo.train())

    } // for

} // Smoothing_FTest3 object
object Smoothing_FTest5 extends App {



    import scala.io.Source
    import Array._
    import scala.util
    var csvFile = "C:\\Users\\Owner\\Desktop\\Life with Divine\\scale.csv"
    var line = ""
    var cvsSplitBy = ","
    var ord =3
    val τ = null
    val rs = new MatrixD(9122,ord)
    val rs1 = new MatrixD(9122,4)

    var row =0
    try
    {

        for (line <- Source.fromFile(csvFile).getLines() ) {

            // use comma as separator
            val value: Array[String] = line.split(cvsSplitBy)
            val t = VectorD.range(0, 4) / 4.0
            val t2 = VectorD.range(0, 1000) / 1000.0
            val y = VectorD.apply(value(0).toDouble, value(1).toDouble, value(2).toDouble, value(3).toDouble)

            rs1.set(row, y)   //rs1 is the raw data

            val moo = new Smoothing_F(y, t, τ, ord)

            val c = moo.train()

            //var i =0
            println(s"y = $y \nt = $t \nc = $c")
            //    for (i <-0 to c.dim-1)  println("test"+c(i))
            println("testc")
          //  moo.print_phi()
                println("test")
              //  moo.printτ()

            rs.set(row, c)   //rs is the coefiiect matrix
            row += 1

            // new FPlot (t, y, t2, moo.predict, s"B-Spline Fit: ord = $ord")


        }

    } catch {
        case ex: IOException => println("Couldn't find that file.")
    }
    import scalation.analytics.KMeansClustering

    val cl = new KMeansClustering (rs.sliceCol(2,3), 10, 0)
    val c2= new MatrixD(2,4)




    for (i<-0 to 2) {
        c2(0, i) = rs.sliceCol(i, i + 1).min()
        c2(1, i) = rs.sliceCol(i, i + 1).max()
    }




    val c3 = new MatrixD(9122,ord)  //generate uniforom population for counting num S
    import scalation.random.Variate
    var uni0 = new Uniform(c2(0,0),c2(1,0))
    var uni1 = new Uniform(c2(0,1),c2(1,1))
    var uni2 = new Uniform(c2(0,2),c2(1,2))
    var uni3 = new Uniform(c2(0,3),c2(1,3))

    for (i<-0 to 9121)
    c3(i,0)= uni0.gen
    for (i<-0 to 9121)
        c3(i,1)= uni1.gen
    for (i<-0 to 9121)
        c3(i,2)= uni2.gen
   // for (i<-0 to 9121)
      ///  c3(i,3)= uni3.gen

   // println(c3.sliceCol(0,1).toString)
    println(rs.toString)    //print the coefiicence
   // println(rs1.toString)   //print the raw records
  // println(c3.toString)    //print the uniform
  //  println ("--- final cluster = " + cl.cluster ().deep + "\n")
/*
    val ca = new KMeansClustering (rs, 10, 0).w()
    val cb = new KMeansClustering (c3, 10, 0).w() //for uniform

    val gap = cb-ca



    print("gap"+gap)

    */



}

