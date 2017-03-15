
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell, John Miller
 *  @version 1.3
 *  @date    Wed Mar  8 11:19:09 2017
 *  @see     LICENSE (MIT style license file).
 */

package scalation.analytics.clusterer

import scala.util.control.Breaks.{breakable, break}
import math.{min, max}

import scalation.linalgebra.{MatrixD, VectorD, VectorI}
import scalation.random.{Discrete, Randi, Uniform, RandomVecD, RandomVecI}
import scalation.util.{banner, Error}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `AffinityPropagationClusterer` class cluster several vectors/points 
 *  using the Adaptive Affinity Propagation clustering technique.  
 *-----------------------------------------------------------------------------
 *  @see https://en.wikipedia.org/wiki/Affinity_propagation
 *  @see https://arxiv.org/abs/0805.1096
 *-----------------------------------------------------------------------------
 *  @param x        the vectors/points to be clustered stored as rows of a matrix
 *  @param k        the number of clusters to make
 *  @param λ        the damping factor (default = 0.5)
 *  @param s        the random number stream (to vary the clusters made)
 */
class AffinityPropagationClusterer (x: MatrixD, k: Int = -1, λ: Double = 0.5, s: Int = 0)
    extends Clusterer with Error
{
    protected val DEBUG    = true                                // debug flag
    protected val n        = x.dim1
    protected val r        = new MatrixD (n, n)                  // "responsibility" matrix
    protected val a        = new MatrixD (n, n)                  // "availability" matrix

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the centroids. Should only be called after `cluster ()`. 
     */
    def centroids (): MatrixD = null

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the sizes of the centroids. Should only be called after 
     *  `cluster ()`. 
     */
    def csize (): VectorI = null

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute a distance metric (distance squared) between vectors/points 'u' and 'v'.
     *  @param u  the first vector/point
     *  @param v  the second vector/point
     */
    override def distance (u: VectorD, v: VectorD): Double =
    {
        (u - v).normSq       // squared Euclidean norm used for efficiency, may use other norms
    } // distance

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the similarity between any two points using the negative of the
     *  distance metric (usually the squared Euclidean distance).
     *  @param i  the index of the first vector/point
     *  @param j  the index of the second vector/point
     */
    private def sim (i: Int, j: Int): Double =
    {
        -distance (x(i), x(j))
    } // sim

    private def updateR ()
    {
        for (i <- 0 until n) {
            var max1e = Double.NegativeInfinity
            var max2e = Double.NegativeInfinity
            var jmax  = 0
            for (j <- 0 until n) {
                for (jp <- 0 until n) {
                    val e = a(i, jp) + sim (i, jp)
                    if (e > max1e)      { max2e = max1e; max1e = e; jmax = j }
                    else if (e > max2e) { max2e = e }
                } // for
            } // for
            for (j <- 0 until n) {
                val maxe = if (j == jmax) max2e else max1e
                r(i, j) = (λ * r(i, j)) + ((1.0 - λ) * (sim(i, j) - maxe))
            } // for
        } // for
        if (DEBUG) println (s"updated r = $r")
    } // updateR

    private def updateA ()
    {
        for (i <- 0 until n) {
            val rp  = new VectorD (n)
            var sum = 0.0
            for (j <- 0 until n) {
                rp(j) = if (r(j, i) < 0 && j != i) 0 else r(j, i)
                sum += rp(j)
            } // for
            for (j <- 0 until n) {
                var newa = sum - rp(j)
                if (newa > 0 & j != i) newa = 0
                a(j, i) = (λ * a(j, i)) + ((1.0 - λ) * newa)
            } // for
        } // for
        if (DEBUG) println (s"updated a = $a")
    } // updateA

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a set of points/vectors, put them in clusters, returning the cluster
     *  assignment vector.  A basic goal is to minimize the sum of the distances
     *  between points within each cluster.
     */
    def cluster (): Array [Int] =
    {
        updateR ()
        updateA ()
        null
    } // cluster

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a new point/vector y, determine which cluster it belongs to.
     *  @param y  the vector to classify
     */
    def classify (y: VectorD): Int = 0
    
} // AffinityPropagationClusterer

object AffinityPropagationClustererTest extends App
{
    val v = new MatrixD ((6, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 5.0, 4.0,
                                 4.0, 5.0,
                                 9.0, 8.0,
                                 8.0, 9.0)

    val k = 3

    println ("v = " + v)
    println ("k = " + k)
    println ("----------------------------------------------------")

    val cl = new AffinityPropagationClusterer (v, k)
    cl.cluster ()

} // AffinityPropagationClustererTest

