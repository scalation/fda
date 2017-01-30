
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell
 *  @version 1.2
 *  @date    Mon Jan 30 15:00:32 EDT 2017
 *  @see     LICENSE (MIT style license file).
 */

package scalation.analytics

import scala.collection.mutable.Set
import scala.util.control.Breaks.{breakable, break}

import scalation.linalgebra.{MatrixD, VectorD}
import scalation.random.{Randi, Uniform}
import scalation.util.{banner, Error}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `TightClustering` class uses k-means to cluster points into tight 
 *  clusters.
 *  @see http://doi.wiley.com/10.1111/j.0006-341X.2005.031032.x
 *  @param x        the vectors/points to be clustered stored as rows of a matrix
 *  @param k        the number of clusters to make
 *  @param s        the random number stream (to vary the clusters made)
 *  @param primary  true indicates use the primary technique for initiating the clustering
 */
class TightClustering (x: MatrixD, k: Int)
      extends Error
{
    if (k >= x.dim1) flaw ("constructor", "k must be less than the number of vectors")

    private val DEBUG    = true                        // debug flag
    private val MAX_ITER = 200                         // the maximum number of iterations
    private val α        = 0.4                         // constant close to zero
    private val B        = 5                           // the number of random k-means samples
    private val cent     = new MatrixD (k, x.dim2)     // the k centroids of tight clusters
    private val comember = new MatrixD (x.dim1, k)     // comembership matrix for k tight clusters

    private def meanComember (): MatrixD =
    {
        val sco = new MatrixD (x.dim1, k)              // sum/total comembership matrix
        for (b <- 0 until B) {
            val kmeans = new KMeansClustering (x, k, b)
            val clustr = kmeans.cluster ()             // randomly cluster
            for (i <- 0 until clustr.length) {         // update sums
                sco(i, clustr(i)) += 1.0
            } // for
        } // for
        sco / B                                        // return average comembership matrix
    } // meanComember

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Iteratively compute tight clusters.
     */
    def cluster (): MatrixD =
    {
        val mco = meanComember()
        println ("average comembership matrix = ")
        println (mco)

        for (i <- mco.range1) {
            for (j <- mco.range2) {
                if (mco(i, j) >= 1.0 - α) println (s"mco($i, $j) >= 1.0 - α")
            } // for
        } // for

        comember                                       // return the cluster assignment vector
    } // cluster

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a new point/vector 'y', determine which cluster it belongs to,
     *  i.e., the cluster whose centroid it is closest to.
     *  @param y  the vector to classify
     */
    def classify (y: VectorD): Int =
    {
        /*
        var dist = distance (y, cent(0))               // calc distance to centroid 0
        var clus = 0                                   // assign y to cluster 0
        for (c <- 1 until k) {
            val newDist = distance (y, cent(c))        // calc distance to centroid c
            if (newDist < dist) {                      // is it closer than old distance
                dist = newDist                         // make it the new distance
                clus = c                               // assign y to cluster c
            } // if
        } // for
        clus                                           // return cluster y belongs to
         */
        0
    } // classify

} // TightClustering class


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `TightClusteringTest` object is used to test the `TightClustering` class.
 *  > run-main scalation.analytics.TightClusteringTest
 */
object TightClusteringTest extends App
{
    val v = new MatrixD ((6, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 5.0, 4.0,
                                 4.0, 5.0,
                                 9.0, 8.0,
                                 8.0, 9.0)
    println ("v = " + v)
    println ("----------------------------------------------------")

    banner ("TightClustering for stream")
    val tcl = new TightClustering (v, 3)
    println ("--- final cluster = " + tcl.cluster () + "\n")

} // TightClusteringTest object
