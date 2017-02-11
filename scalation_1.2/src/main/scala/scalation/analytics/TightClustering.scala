
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell
 *  @version 1.2
 *  @date    Mon Jan 30 15:00:32 EDT 2017
 *  @see     LICENSE (MIT style license file).
 */

package scalation.analytics

import scala.collection.mutable.Set
import scala.util.control.Breaks.{breakable, break}

import scalation.linalgebra.{MatrixD, SparseMatrixD, VectorD}
import scalation.random.{Randi, Uniform}
import scalation.util.{banner, Error}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `TightClustering` class uses k-means to cluster points into tight 
 *  and stable clusters without forcing all points into clusters via the
 *  Tseng-Wong algorithm.
 *  @see http://doi.wiley.com/10.1111/j.0006-341X.2005.031032.x
 *  @param x  the vectors/points to be clustered stored as rows of a matrix
 *  @param k  the number of clusters to make
 *  @param α  the tightness of clusters; a constant close to zero (default = 0.0)
 *  @param β  the stableness of clusters; a constant close to one (default = 0.7)
 *  @param B  the number of random k-means samples (default = 10)
 */
class TightClustering (x: MatrixD, k: Int, α: Double = 0.0, β: Double = 0.7, B: Int = 10)
      extends Error
{
    if (k >= x.dim1) flaw ("constructor", "k must be less than the number of vectors")
    if (α < 0)       flaw ("constructor", "α must be greater than or equal to 0")
    if (β > 1)       flaw ("constructor", "β must be less than or equal to 1")

    private val DEBUG    = true                         // debug flag
    private val MAX_ITER = 200                          // the maximum number of iterations
    private val cent     = new MatrixD (k, x.dim2)      // the k centroids of tight clusters

    if (DEBUG) {
        banner (s"TightClustering(x, k=$k, α=$α, β=$β, B=$B)")
        println (s"x.dim1 = ${x.dim1}")
        println (s"x.dim2 = ${x.dim2}")
    } // if

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the average or mean comembership matrix of `B`-many
     *  independent random subsample clusterings via the K-Means clustering
     *  algorithm.
     */
    private def meanComember (): SparseMatrixD =
    {
        if (DEBUG) banner ("TightClustering: meanComember()")
        val co = new SparseMatrixD (x.dim1, x.dim1)     // comembership matrix
        for (b <- 0 until B) {
            val kmeans = new KMeansClustering (x, k, b) // k-means clusterer
            val clustr = kmeans.cluster ()              // randomly cluster
            for (i <- co.range1; j <- co.range2) {      // update comembership matrix
                if (clustr(i) == clustr(j)) co(i, j) += 1.0 
            } // for
        } // for
        co /= B                                         // average comembership matrix
        co
    } // meanComember

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Similarity measure of two sets of points. `sim(u, v)=1` if and only if 
     *  sets `u` and `v` are identical.
     *  @param u  first set of points
     *  @param v  second set of points
     */    
    private def sim (u: Set[VectorD], v: Set[VectorD]): Double = (u & v).size / (u | v).size

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Search for a set of points `V = {v1, ... , vm} ⊂ {1, ... ,n}` such that
     *  `dbar(i, j) >= 1 − α`, for all `i`, `j` where `α` is a constant close
     *  to 0. Order sets with this property by size to obtain candidates of 
     *  tight clusters.
     */
    private def tightCandidates (dbar: SparseMatrixD): Set[Set[VectorD]] = 
    {
        val sets = Set.empty[Set[VectorD]]
        println (dbar)
        for (i <- dbar.range1) {
            val set = Set(x(i))
            for (j <- dbar.range2) if (dbar(i, j) >= 1.0 - α) set += x(j)
            sets += set
        } // for
        if (DEBUG) for (set <- sets) println (s"tight cluster candidate = $set")
        sets
    } // tightCandidates

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Iteratively compute tight clusters.
     */
    def cluster (): MatrixD =
    {
        val dbar = meanComember ()
        val v    = tightCandidates (dbar)
        null
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

    val tcl = new TightClustering (v, 3, 0.1)
    println ("--- final cluster = " + tcl.cluster () + "\n")

} // TightClusteringTest object

