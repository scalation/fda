
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell, John Miller, Hao Peng
 *  @version 1.3
 *  @date    Tue Mar  7 22:10:21 2017
 *  @see     LICENSE (MIT style license file).
 *
 *  @see webee.technion.ac.il/people/koby/publications_all/pdfs//conf_ijcai_SlonimAC13.pdf
 *  @see cseweb.ucsd.edu/~avattani/papers/hartigan.pdf
 */

package scalation.analytics.clusterer

import scala.collection.mutable.Set
import scala.util.control.Breaks.{breakable, break}

import scalation.linalgebra.{MatrixD, VectorD, VectorI}
import scalation.random.{Discrete, Randi, Uniform, PermutedVecI, RandomVecI}
import scalation.util.{banner, Error}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Algorithm` object specifies which algorithm to use.
 */
object Algorithm extends Enumeration
{
    type Algorithm = Value
    val HARTIGAN, LLOYD = Value
} // Algorithm

import Algorithm._

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansPlusPlusClusterer` class cluster several vectors/points using 
 *  the k-means++ clustering technique.  
 *-----------------------------------------------------------------------------
 *  @see http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf 
 *-----------------------------------------------------------------------------
 *  @param x        the vectors/points to be clustered stored as rows of a matrix
 *  @param k        the number of clusters to make
 *  @param s        the random number stream (to vary the clusters made)
 */
class KMeansPlusPlusClusterer (x: MatrixD, k: Int, algo: Algorithm = LLOYD, s: Int = 0)
    extends Clusterer with Error
{
    protected val DEBUG    = false                               // debug flag
    protected val MAX_ITER = 200                                 // the maximum number of iterations
    protected var _k       = 0                                   // last centroid chosen 
    protected val cent     = new MatrixD (k, x.dim2)             // the k centroids of clusters
    protected val pdf      = new VectorD (x.dim1)                // pdf for choosing centroids
    protected val sizes    = new VectorI (k)                     // the cluster sizes
    protected val clustr   = Array.ofDim [Int] (x.dim1)          // assignment of vectors to clusters
//    protected val dist     = new VectorD (x.dim1)                // distance to closest centroid
    protected val raniv    = PermutedVecI (VectorI.range (0, x.dim1), s)
    protected val live     = Set((0 until k).toSeq:_*)           // live set of clusters

//    dist.set (Double.PositiveInfinity)

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the sizes of the centroids. Should only be called after 
     *  `cluster ()`. 
     */
    def csize (): VectorI = sizes

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the centroids. Should only be called after `cluster ()`. 
     */
    def centroids (): MatrixD = cent    

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the shortest distance from a data point 'u' to the closest 
     *  centroid we have already chosen.
     *  @param u  the vector/point
     */
    private def shortestDistance (u: VectorD): Double =
    {
        var min = Double.PositiveInfinity
        for (i <- 0 to _k) {
            val d = distance (u, cent(i))
            if (d < min) min = d
        } // for
        min
    } // shortestDistance

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the sum of the shortest distances from each data point to 
     *  its closest centroid.
     */
    private def pseudoSSE (): Double =
    {
        var sum = 0.0
        for (i <- x.range1) sum += shortestDistance (x(i))
        sum
    } // sse

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Update the probability distribution for choosing the next centroid.
     */
    private def updatePDF ()
    {
        val psse = pseudoSSE ()
        for (i <- x.range1) pdf(i) = shortestDistance (x(i)) / psse
        if (DEBUG) println (s"updated pdf = $pdf")
    } // updatePDF

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Initialize the centroids according to the k-means++ technique.
     */
    private def initCentroids ()
    {
        val ran1 = new Randi (0, x.dim1-1, s)                    // uniform integer distribution 
        cent(0)  = x(ran1.igen)                                  // pick first centroid uniformly at random
        for (i <- 1 until k) {                                   // pick remaining centroids
            updatePDF ()                                         // update probability distribution
            val ran2 = Discrete (pdf, stream = (s+i) % 1000)
            cent(i)  = x(ran2.igen)                              // pick next centroid according to pdf
            _k       = i                                         // update last centroid chosen
        } // for
        if (DEBUG) println (s"initial cent = $cent")
    } // initCentroids

    /*
    private def reassign (first: Boolean = false): Boolean =
    {
        println (s"reassign(first = $first")
        var done = true                                          // done indicates no changes
        for (i <- x.range1) {
            val v = x(i)                                         // let v be the ith vector
            val curDist = if (first) Double.PositiveInfinity
                          else       distance (v, cent(clustr(i))) // current distance from centroid
            println (s"curDist = $curDist; dist($i) = ${dist(i)}")
            for (c <- 0 until k) {
                val newDist = distance (v, cent(c))              // calc distance to centroid c

//                if (newDist < dist(i)) {                         // is it closer than old distance
//                    dist(i)           = newDist                  // make it the new distance
                if (newDist < curDist) {
                    sizes(clustr(i)) -= 1                        // decrement size of previous cluster
                    clustr(i)         = c                        // reassign vector x(i) to cluster c
                    sizes(c)         += 1                        // increment size of cluster c
                    done              = false                    // changed clusters => not done
                } // if
            } // for
        } // for
        //for (c <- 0 until k) sizes(c) = 0
        //for (i <- x.range1) sizes(clustr(i)) += 1                // TODO fix
        done                                                     // return whether there were no changes
    } // reassign
     */

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reassign each vector/point to the cluster with the closest centroid.
     *  Indicate done, if no points changed clusters (for stopping rule).
      */
    private def reassign (first: Boolean = false): Boolean =
    {
        var done = true
        for (i <- x.range1) {
            if (first) {
                val cc     = closest (x(i))
                clustr(i)  = cc
                sizes(cc) += 1                
                done       = false
            } else {
                val c0 = clustr(i)
                if (sizes(c0) > 1) {
                    val cc  = closest (x(i))
                    if (cc != c0) {
                        clustr(i)  = cc
                        sizes(c0) -= 1
                        sizes(cc) += 1
                        done       = false
                    } // if
                } // if
            } // if
        } // for
        for (c <- 0 until k if sizes(c) == 0) {
            assert (false, s"reassign (first = $first): empty cluster c = $c; sizes = $sizes")
        } // if
        done
    } // reassign

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Calculate the centroids based on current assignment of points to clusters.
     */
    private def calcCentroids ()
    {
        val cx = new MatrixD (k, x.dim2)                         // to hold sum of vectors for each cluster
        val cs = new VectorD (k)                                 // to hold number of vectors in each cluster
        for (i <- x.range1) {
            val ci  = clustr(i)                                  // x(i) currently assigned to cluster ci
            cx(ci)  = cx(ci) + x(i)                          // add the next vector in cluster
//            cx(ci) += x(i)                          // add the next vector in cluster            
            cs(ci) += 1.0                                    // add 1 to number in cluster
        } // for
        for (c <- 0 until k) cent(c) = cx(c) / cs(c)             // divide to get averages/means
    } // calcCentroids
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Cluster the points using Llooyd's algorithm.
     */
    private def clusterLloyd (): Array [Int] =
    {
        breakable { for (l <- 1 to MAX_ITER) {
            if (reassign (l==1)) break                           // reassign points to clusters (no change => break)
            calcCentroids ()                                     // re-calculate the centroids
            if (DEBUG) {
                println ("(" + l + ") clustr = " + clustr.deep)
                println ("(" + l + ") cent   = " + cent)
            } // if
        }} // for
        if (DEBUG) println (s"clustr = ${clustr.deep}")
        clustr
    } // clusterLloyd

    /** Cluster the points using a simplified version of the Hartigan-Wong 
     *  algorithm.
     *  @see http://www.tqmp.org/RegularArticles/vol09-1/p015/p015.pdf
     */
    def clusterHartigan (): Array [Int] =
    {
        reassign (true)
        //for (c <- 0 until k if sizes(c) == 0) assert (false, s"1!!!empty cluster c = $c; sizes = $sizes")
        calcCentroids ()
        breakable { for (l <- 1 to MAX_ITER) {
            if (reassign2 ()) break                              // reassign points to clusters (no change => break)
            if (DEBUG) {
                println ("(" + l + ") clustr = " + clustr.deep)
                println ("(" + l + ") cent   = " + cent)
            } // if
            //println (s"sizes = $sizes")
            //for (c <- 0 until k if sizes(c) == 0) assert (false, s"2!!!empty cluster c = $c; l = $l")
        }} // for
        if (DEBUG) println (s"clustr = ${clustr.deep}")
        clustr
    } // clusterHartigan
    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reassign each vector/point to the cluster with the closest centroid.
     *  This one follows a simplified version of the Hartigan-Wong algorithm.
     *  Indicate done, if no points changed clusters (for stopping rule).
     *  @see http://www.tqmp.org/RegularArticles/vol09-1/p015/p015.pdf
     */
    private def reassign2 (): Boolean =
    {
        var done = true                                          // done indicates no changes
        val iv   = raniv.igen                                    // randomly order the points
        for (i <- iv) {
            val c0     = clustr(i)                               // record point v's current cluster (c0)
//          println (s"reaassign2: try point $i in cluster $c0")
            if (sizes(c0) > 1) {                                 // make sure not to empty a cluster
                calcCentroids ()                                 // re-calculate the centroids
                val c1 = closestByR2 (i)                         // find closest cluster to point v
                sizes(c0) -= 1                                   // decrement the size of cluster c0
                clustr(i)  = c1                                  // reassign point v to cluster c1
                sizes(c1) += 1                                   // increment size of cluster c1
                if (c1 != c0) done = false                       // changed clusters => not done
            } // if
        } // for
        for (c <- 0 until k if sizes(c) == 0) assert (false, s"reassign2 (): empty cluster c = $c; sizes = $sizes")
        done                                                     // return whether there were no changes
    } // reassign2

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the closest cluster to point 'u'.
     */
    private def closest (u: VectorD): Int =
    {
        var cmin = 0; var min = distance (u, cent(0))
        for (c <- 1 until k) {
            val d = distance (u, cent(c))
            if (d < min) { cmin = c; min = d }
        } // for
        cmin
    } // closest

   //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the closest cluster to point 'u' according to the R2 value 
     *  described in the Hartigan-Wong algorithm.
     */
    private def closestByR2 (i: Int): Int =
    {
        val u    = x(i)
        var cmin = 0; var min = distance (u, cent(0))
        var r2   = 0.0
        val cc   = clustr(i)
        for (c <- 1 until k) {
            r2 = if (c == cc) (sizes(c) * distance (u, cent(c))) / (sizes(c) - 1)
                 else         (sizes(c) * distance (u, cent(c))) / (sizes(c) + 1)
            if (r2 < min) { cmin = c; min = r2 }
        } // for
        cmin
    } // closestByR2

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the closest 2 centroids to the point `u`, respectively.
     *  @param u  the vector/point
     */
    private def closest2 (u: VectorD): Array [Int] = 
    {
        // (0 until k).sortBy { case (l, r) => distance (u, centr(l)) < distance (u, centr(r)) }.take (2)
        var cmin1 = 0; var min1 = distance (u, cent(0))
        var cmin2 = 1; var min2 = distance (u, cent(1))
        if (cmin2 < cmin1) {
            val (ctemp, temp) = (cmin1, min1)
            cmin1 = cmin2
            min1  = min2
            cmin2 = ctemp
            min2  = temp
        } // if
        for (c <- 2 until k) {
            val d = distance (u, cent(c))
            if (d < min1) {
                min2  = min1
                cmin2 = cmin1
                min1  = d
                cmin1 = c
            } else if (d < min2) {
                min2  = d
                cmin2 = c
            } // if
        } // for
        Array (cmin1, cmin2)
    } // closest2

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a set of points/vectors, put them in clusters, returning the cluster
     *  assignment vector.  A basic goal is to minimize the sum of the distances
     *  between points within each cluster.
     */
    def cluster (): Array [Int] =
    {
        if (clustered) return clustr
        clustered = true
        initCentroids ()
        algo match {
            case LLOYD         => clusterLloyd ()
            case HARTIGAN      => clusterHartigan ()
        } // match
    } // cluster

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a new point/vector y, determine which cluster it belongs to.
     *  @param y  the vector to classify
     */
    def classify (y: VectorD): Int = 0 // TODO fix

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the sum of squared errors (distance sqaured from centroid `c`
     *  for all points) for centroid `c`.
     */
    def sse (c: Int): Double =
    {
        var sum = 0.0
        for (i <- x.range1) {
            val cli = clustr(i)
            if (cli == c) sum += distance (x(i), cent(cli))
        } // for
        sum
    } // sse

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the sum of squared errors (distance sqaured from centroid for all points)
     */
    def sse (): Double = sse (x)

} // KMeansPlusPlusClusterer

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** Companion object for the `KMeansPlusPlusClusterer` class.
 */
object KMeansPlusPlusClusterer
{

    private var streams = VectorI.range (0, 1000)

    def permuteStreams (stream: Int = 0)
    {
        import scalation.random.PermutedVecI
        streams = PermutedVecI (streams, stream).igen
    } // permuteStreams

    def run (x: MatrixD, k: Int, algo: Algorithm = LLOYD): (KMeansPlusPlusClusterer, Array [Int]) =
    {
//        println (s"run(stream = $streams)")
        var iters  = 0
        val check  = 3
        var seen   = 0
        var ssemin = Double.PositiveInfinity
        var kmppmin: KMeansPlusPlusClusterer = null
        for (s <- streams) {
            iters += 1
            val kmpp = new KMeansPlusPlusClusterer (x, k, algo, s = s)
            val cls = kmpp.cluster ()
            val sse = kmpp.sse ()
            //println (s"sse = $sse")
            if (sse == ssemin) {
                if (seen==check) return (kmppmin, kmppmin.clustr) // { println (s"iters = $iters"); return (kmppmin, kmppmin.clustr) }
                else             seen += 1
            } // if
            if (sse < ssemin) { ssemin = sse; kmppmin = kmpp; seen = 0 }
        } // for
//        { println (s"iters = $iters"); return (kmppmin, kmppmin.clustr) }
        return (kmppmin, kmppmin.clustr)
    } // run

} // KMeansPlusPlusClusterer object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansPlusPlusClustererTester` trait includes a `test` function to aid
 *  in the testing of the `KMeansPlusPlusClusterer` class. 
 */
trait KMeansPlusPlusClustererTester
{
    import scalation.stat.Statistic
    import scalation.plot.Plot

    def checkEmpty (cls: Array [Int], k: Int)
    {
        for (c <- 0 until k if !cls.contains(c)) assert(false, s"empty cluster found c = $c")
    } // checkEmpty

    def test (v: MatrixD, k: Int, opt: Double = -1, plot: Boolean = false, nstreams: Int = 1000)
    {
        banner (s"Testing KMeansPlusPlusCluster")
        println (s"v.dim1 = ${v.dim1}")
        println (s"v.dim2 = ${v.dim2}")
        println (s"     k = $k")
        if (plot) new Plot (v.col(0), v.col(1))
        for (algo <- Algorithm.values) {
            banner (s"test (algo = $algo)")
            val statSSE = new Statistic ("sse")
            var ok      = 0
            for (s <- 0 until nstreams) {
                val cl  = new KMeansPlusPlusClusterer (v, k, algo, s)
                val cls = cl.cluster ()
                val sse = cl.sse ()
                checkEmpty (cls, k)
                // println (s"stream $s, sse = $sse")
                statSSE.tally (sse)
                if ((opt != -1) && (sse <= opt)) ok += 1
            } // for
            println (Statistic.labels)
            println (statSSE)
            println (s"min sse = ${statSSE.min}")
            if (opt != -1) println (s"optimal = $ok / $nstreams")            
        } // for 
    } // test

    def test2 (v: MatrixD, k: Int, opt: Double = -1, plot: Boolean = false, nstreams: Int = 1000)
    {
        banner (s"Testing KMeansPlusPlusCluster object")
        println (s"v.dim1 = ${v.dim1}")
        println (s"v.dim2 = ${v.dim2}")
        println (s"     k = $k")
        if (plot) new Plot (v.col(0), v.col(1))
        for (algo <- Algorithm.values) {
            banner (s"test2 (algo = $algo)")
            val statSSE = new Statistic ("sse")
            var ok      = 0
            for (s <- 0 until nstreams) {
                KMeansPlusPlusClusterer.permuteStreams (s)
                val (cl, cls)  = KMeansPlusPlusClusterer.run (v, k, algo)
                val sse = cl.sse ()
                checkEmpty (cls, k)
                // println (s"stream $s, sse = $sse")
                statSSE.tally (sse)
                if ((opt != -1) && (sse <= opt)) ok += 1
            } // for
            println (Statistic.labels)
            println (statSSE)
            println (s"min sse = ${statSSE.min}")
            if (opt != -1) println (s"optimal = $ok / $nstreams")            
        } // for 
    } // test2
    

} // KMeansPlusPlusClutererTester


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansPlusPlusClustererTest` object is used to test the `KMeansPlusPlusCluterer`
 *  class.
 *  > run-main scalation.analytics.clusterer.KMeansPlusPlusClustererTest
 */
object KMeansPlusPlusClustererTest
       extends App with KMeansPlusPlusClustererTester
{
    val v  = new MatrixD ((6, 2), 1.0, 2.0,
                                  2.0, 1.0,
                                  5.0, 4.0,
                                  4.0, 5.0,
                                  9.0, 8.0,
                                  8.0, 9.0)
    val k   = 3
    val opt = 3.0
    test (v, k, opt)
    test2 (v, k, opt)

} // KMeansPlusplusclustererTest


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansPlusPlusClustererTest2` object is used to test the `KMeansPlusPlusCluterer`
 *  class.
 *  > run-main scalation.analytics.clusterer.KMeansPlusPlusClustererTest2
 */
object KMeansPlusPlusClustererTest2
       extends App with KMeansPlusPlusClustererTester
{
    import scalation.random.{Normal, Bernoulli}
    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 1.0)
    val dist2 = Normal (8.0, 1.0)
    val v     = new MatrixD (50, 2)
    for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen)
    val k   = 4
    val opt = 76         // rounded up

    test (v, k, opt)
    test2 (v, k, opt)

} // KMeansPlusplusclustererTest2

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansPlusPlusClustererTest3` object is used to test the `KMeansPlusPlusCluterer`
 *  class.
 *  > run-main scalation.analytics.clusterer.KMeansPlusPlusClustererTest3
 */
object KMeansPlusPlusClustererTest3
       extends App with KMeansPlusPlusClustererTester
{
    import scalation.random.{Normal, Bernoulli}
    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 1.0)
    val dist2 = Normal (8.0, 1.0)
    val v     = new MatrixD (50, 3)
    for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen)

    val k   = 8
    val opt = 106         // rounded up

    test (v, k, opt)
    test2 (v, k, opt)

} // KMeansPlusplusclustererTest3

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansPlusPlusClustererTest4` object is used to test the `KMeansPlusPlusCluterer`
 *  class.
 *  > run-main scalation.analytics.clusterer.KMeansPlusPlusClustererTest4
 */
object KMeansPlusPlusClustererTest4
       extends App with KMeansPlusPlusClustererTester
{
    import scalation.random.{Normal, Bernoulli}
    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 1.0)
    val dist2 = Normal (8.0, 1.0)
    val v     = new MatrixD (100, 4)
    for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen)

    val k   = 16
    val opt = 290         // rounded up

    test (v, k, opt)
    test2 (v, k, opt)

} // KMeansPlusplusclustererTest4

