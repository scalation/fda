
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell, John Miller, Hao Peng
 *  @version 1.3
 *  @date    Thu Mar  9 15:08:30 2017
 *  @see     LICENSE (MIT style license file).
  */

package scalation.analytics.clusterer

import math.log

import scalation.linalgebra.{MatrixD, SVD, VectorD}
import scalation.random.RandomVecD

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `GapStatistic` object is used to help determine the optimal number
 *  of clusters for a clusterer by comparing results to a reference 
 *  distribution.
 *-----------------------------------------------------------------------------
 *  @see https://web.stanford.edu/~hastie/Papers/gap.pdf 
 */
object GapStatistic
{

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute a reference distribution based on a set of points.
     *  @param x        the vectors/points to be clustered stored as rows of a matrix
     *  @param useSVD   use SVD to account for the shape of the points (default = true)
     *  @param s        the random number stream (to vary the clusters made)
     */
    def reference (x: MatrixD, useSVD: Boolean = true, stream: Int = 0): MatrixD =
    {
        var ref = new MatrixD (x.dim1, x.dim2)
        if (useSVD) {
            val mean  = x.mean
            val xzero = x - mean
            val svd   = new SVD (xzero)
            val (u, s, vt) = svd.factor ()
            val xp    = xzero * vt.t
            val zp    = new MatrixD (x.dim1, x.dim2)
            for (i <- zp.range2) {
                val ci = xp.col(i)
                zp.setCol(i, RandomVecD (zp.dim1, ci.max, ci.min, stream = (stream + i) % 1000).gen)
            } // for
            ref = (zp * vt) + mean
        } else {
            for (i <- ref.range2) {
                val ci = x.col(i)
                ref.setCol(i, RandomVecD (ref.dim1, ci.max, ci.min, stream = (stream + i) % 1000).gen)
            } // for
        } // if
        ref 
    } // reference


    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute a sum of pairwise distances between points in each cluster (in
     *  one direction). 
     *  @param x        the vectors/points to be clustered stored as rows of a matrix
     *  @param cl       the `Clusterer` use to compute the distance metric 
     *  @param clustr   the cluster assignments
     *  @param k        the number of clusters
     */
    def cumDistance (x: MatrixD, cl: Clusterer, clustr: Array[Int], k: Int): VectorD =
    {
        val sums   = new VectorD (k)
        for (i <- 0 until x.dim1-1; j <- i+1 until x.dim1 if clustr(i) == clustr(j)) {
            sums(clustr(j)) += cl.distance (x(i), x(j))
        } // for
        sums
    } // cumDistance

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the within sum of squares (in one direction).
     *  @param x        the vectors/points to be clustered stored as rows of a matrix
     *  @param cl       the `Clusterer` use to compute the distance metric 
     *  @param clustr   the cluster assignments
     *  @param k        the number of clusters
     */
    def withinSSE (x: MatrixD, cl: Clusterer, clustr: Array[Int], k: Int): Double =
    {
        println (cl.csize())
        (cumDistance (x, cl, clustr, k) / cl.csize().toDouble).sum
    } // withinSSE

} // GapStatistic

object GapStatisticTest extends App
{

    val v  = new MatrixD ((6, 2), 1.0, 2.0,
                                  2.0, 1.0,
                                  5.0, 4.0,
                                  4.0, 5.0,
                                  9.0, 8.0,
                                  8.0, 9.0)

    v += 100

/*    
    import scalation.linalgebra.VectorD
    import scalation.random.{Normal, Bernoulli}
    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 0.1)
    val dist2 = Normal (8.0, 0.1)
    val v     = new MatrixD (50, 2)
    for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen)
 */

    val ref1 = GapStatistic.reference (v, useSVD = false)
    val ref2 = GapStatistic.reference (v, useSVD = true)    

    println (s"   v = $v")
    println (s"ref1 = $ref1")
    println (s"ref2 = $ref2")    

    import scalation.plot.Plot

    new Plot (v.col(0), v.col(1), _title = "Original")
    new Plot (ref1.col(0), ref1.col(1), _title = "Reference (no SVD)")
    new Plot (ref2.col(0), ref2.col(1), _title = "Reference (with SVD)")        

} // GapStatisticTest


object GapStatisticTest2 extends App
{
    import Algorithm._

/*    
    val v  = new MatrixD ((6, 2), 1.0, 2.0,
                                  2.0, 1.0,
                                  5.0, 4.0,
                                  4.0, 5.0,
                                  9.0, 8.0,
                                  8.0, 9.0)

 */

    import scalation.linalgebra.VectorD
    import scalation.random.{Normal, Bernoulli}
    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 0.1)
    val dist2 = Normal (8.0, 0.1)
    val v     = new MatrixD (50, 2)
    for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen)
    

    val kMax  = 5

    val awk = new VectorD (kMax)
    val rwk = new VectorD (kMax)
    val gap = new VectorD (kMax)
    val kv  = VectorD.range(1, kMax+1)

    for (k <- 0 until kMax) {
        val ref        = GapStatistic.reference (v, useSVD = false)
        val (acl, acls)  = KMeansPlusPlusClusterer.run (v, k+1, HARTIGAN)
        val (rcl, rcls)  = KMeansPlusPlusClusterer.run (ref, k+1, HARTIGAN)        
        awk(k) = log(GapStatistic.withinSSE (v, acl, acls, k+1))
        rwk(k) = log(GapStatistic.withinSSE (ref, rcl, rcls, k+1))
        gap(k) = rwk(k) - awk(k)
        if (k != 0) {
            if (gap(k-1) >= gap(k) - gap(k-1)*0.1) println (s"optimal k = ${k}")
        } // if
    } // for

    println (s"log_awk = $awk")
    println (s"log_rwk = $rwk")
/*
    val log_awk = awk.map(log(_))
    val log_rwk = rwk.map(log(_))
    val gap     = log_rwk - log_awk
 */
    println (s"gap = $gap")    

    import scalation.plot.Plot
    new Plot (kv, awk, rwk)

} // GapStatisticTest2


