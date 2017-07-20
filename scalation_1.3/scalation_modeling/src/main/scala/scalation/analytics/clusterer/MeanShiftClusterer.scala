
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell, John Miller
 *  @version 1.3
 *  @date    Tue Mar  7 22:10:21 2017
 *  @see     LICENSE (MIT style license file).
 */

package scalation.analytics.clusterer

import scala.collection.mutable.{Map, Set}
import scala.util.control.Breaks.{breakable, break}

import scalation.linalgebra.{MatrixD, VectorD, VectorI}
import scalation.random.{Discrete, Randi, Uniform, PermutedVecI, RandomVecI}
import scalation.util.{banner, Error}

object Kernel extends Enumeration
{
    type Kernel = Value
    val GAUSSIAN, FLAT = Value
} // Kernel

import Kernel._

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `MeanShiftClusterer` class cluster several vectors/points using 
 *  the mean shift clustering technique.  
 *-----------------------------------------------------------------------------
 *  @see https://en.wikipedia.org/wiki/Mean_shift
 *  @see https://doi.org/10.1109/34.400568
 *-----------------------------------------------------------------------------
 *  @param x     the vectors/points to be clustered stored as rows of a matrix
 *  @param h     the bandwidth parameter 
 *  @param s     the random number stream (to vary the clusters made)
 */
class MeanShiftClusterer (x: MatrixD, h: Double = 0.5, kernel: Kernel = GAUSSIAN, s: Int = 0)
    extends Clusterer with Error
{
    protected val DEBUG    = true                                // debug flag
    protected val MAX_ITER = 200                                 // the maximum number of iterations
    protected val EPS      = 1E-12                               // epsilon threshold
    protected val clustr   = Map [VectorD, Set [Int]] ()         // clusters

    val ker: (VectorD => Double) = kernel match {
        case GAUSSIAN => { v => math.exp (- v.normSq / (2 * h * h)) }
        //case GAUSSIAN => { v => math.exp (- v.normSq / (hxc)) }
        case FLAT     => { v => if (v.normSq <= h) 1.0 else 0.0 }
    } // ker

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the weighted mean of the density.
     *  @param v the point
     */
    def m (v: VectorD): VectorD =
    {
        var sum1 = new VectorD (x.dim2)
        var sum2 = 0.0
        for (i <- x.range1 if ker (x(i) - v) != 0) sum1 += x(i) * ker (x(i) - v)
        for (i <- x.range1 if ker (x(i) - v) != 0) sum2 += ker (x(i) - v)
        sum1 / sum2
    } // m

    def centroids(): scalation.linalgebra.MatrixD = ???
    def csize(): scalation.linalgebra.VectorI = ???

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a set of points/vectors, put them in clusters, returning the cluster
     *  assignment vector.
     */
    def cluster (): Array [Int] =
    {
        if (DEBUG) println (x)
        breakable { for (l <- 1 to MAX_ITER) {
            var done = true
            if (DEBUG) println (s"######## iter l = $l")
            for (i <- x.range1) {
                val v = x(i)
                x(i) = m(v)
                if ((x(i) - v).normSq > EPS) done = false
            } // for
            if (DEBUG) println (x)

            if (done) break
/*            if (done) {
                for (i <- x.range1; j <- x.range1) {
                    if ((x(i) - x(j)).normSq <= EPS) {
                        val v = x(i)
                        clustr.getOrElseUpdate (v, Set.empty [Int])
                        clustr(v) += j
                    } // if
                } // for
                break
            } // if
 */
        }} // for

        if (DEBUG) {
            println ("######## result")
            println (clustr)
        } // if
        null
    } // cluster

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a new point/vector y, determine which cluster it belongs to.
     *  @param y  the vector to classify
     */
    def classify (y: VectorD): Int = 0

} // MeanShiftClusterer

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `MeanShiftClustererTest` object is used to test the 
 *  `MeanShiftClusterer` class.
 *  > run-main scalation.analytics.clusterer.MeanShiftClustererTest
 */
object MeanShiftClustererTest extends App
{
    val x  = new MatrixD ((6, 2), 1.0, 2.0,
                                  2.0, 1.0,
                                  5.0, 4.0,
                                  4.0, 5.0,
                                  9.0, 8.0,
                                  8.0, 9.0)

    //val mean = x.mean
    //for (i <- x.range1) x(i) = x(i) - mean

    val msc = new MeanShiftClusterer (x, kernel = GAUSSIAN)
    // val msc = new MeanShiftClusterer (x, h = 2.0, kernel = FLAT)
    msc.cluster ()

} // MeanShiftClustererTest


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `MeanShiftClustererTest2` object is used to test the 
 *  `MeansShiftCluterer` class.
 *  > run-main scalation.analytics.clusterer.MeanShiftClustererTest2
 */
object MeanShiftClustererTest2 extends App
{
    import scalation.random.{Normal, Bernoulli}
    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 1.0)
    val dist2 = Normal (8.0, 1.0)
    val v     = new MatrixD (10, 2)
    for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen)

    val msc = new MeanShiftClusterer (v, h = 1.0, kernel = FLAT)
    msc.cluster ()

} // MeanShiftClustererTest2


