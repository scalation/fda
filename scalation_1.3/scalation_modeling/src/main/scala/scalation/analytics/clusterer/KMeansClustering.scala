
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.3
 *  @date    Tue May 29 14:45:32 EDT 2012
 *  @see     LICENSE (MIT style license file).
 */

package scalation.analytics.clusterer

import scala.collection.mutable.Set
import scala.util.control.Breaks.{breakable, break}

import scalation.math.double_exp
import scalation.linalgebra.{MatrixD, VectorD, VectorI}
import scalation.random.{Randi, Uniform, RandomVecD, RandomVecI}
import scalation.util.{banner, Error}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClustering` class cluster several vectors/points using k-means
 *  clustering.  Either (1) randomly assign points to 'k' clusters or (2) randomly
 *  pick 'k' points as initial centroids (technique (1) to work better and is the
 *  primary technique).  Iteratively, reassign each point to the cluster containing
 *  the closest centroid.  Stop when there are no changes to the clusters.
 *-----------------------------------------------------------------------------
 *  @param x        the vectors/points to be clustered stored as rows of a matrix
 *  @param k        the number of clusters to make
 *  @param s        the random number stream (to vary the clusters made)
 *  @param primary  true indicates use the primary technique for initiating the clustering
 *  @param remote   whether to take a maximally remote or a randomly selected point
 *  @param post     whether to perform post processing by randomly swapping points to reduce error
 */
class KMeansClustering (x: MatrixD, k: Int, s: Int = 0, primary: Boolean = true, remote: Boolean = true, post: Boolean = true)
      extends Clusterer with Error
{
    if (k >= x.dim1) flaw ("constructor", "k must be less than the number of vectors")

    protected val DEBUG    = false                               // debug flag
    protected val MAX_ITER = 200                                 // the maximum number of iterations
    protected val cent     = new MatrixD (k, x.dim2)             // the k centroids of clusters
    protected val sizes    = new VectorI (k)                     // the cluster sizes
    protected val myDist   = new VectorD (x.dim1)                // distances from centroids
    protected val clustr   = Array.ofDim [Int] (x.dim1)          // assignment of vectors to clusters
    protected val dist     = new VectorD (x.dim1)                // distance to closest centroid

    var tc1: Double = 0.0
    var tc2: Double = 0.0

    dist.set (Double.PositiveInfinity)

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
    /** Randomly assign each vector/point 'x(i)' to a random cluster.
     *  Primary technique for initiating the clustering.
     */
    def assign ()
    {
        val ran = new Randi (0, k-1, s)                          // for random integers: 0, ..., k-1
        for (i <- x.range1) {
            clustr(i) = ran.igen                                 // randomly assign x(i) to a cluster
            sizes(clustr(i)) += 1                                // increment size of this cluster
        } // for
    } // assign

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Randomly pick vectors/points to serve as the initial 'k' centroids (cent).
     *  Secondary technique for initiating the clustering.
     */
    def pickCentroids ()
    {
        val rvi = RandomVecI (k, x.dim1-1, 0, stream = s).igen   // random vector of integers
        for (i <- 0 until k) cent(i) = x(rvi(i))                 // set the centroids
    } // pickCentroids

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Fix any empty clusters by taking a point from the largest cluster.

     */
    private def fixEmptyClusters (useDistance: Boolean = remote)
    {
        if (DEBUG) {
            banner (s"fixemptyclusters()")
            println (s"remote = $remote")
            println (s"Initial clustering = ${clustr.deep}")
        } // if
        for (c <- 0 until k if ! (clustr contains c)) {                  // for each empty cluster
            if (DEBUG) println (s"Cluster c=$c is empty!")
            val biggest = sizes.argmax ()                                // biggest cluster
            val indices = clustr.indices.filter (clustr(_) == biggest)   // indices of elements in biggest cluster            
            if (DEBUG) {
                println (s"Current cluster sizes = $sizes")
                println (s"Biggest cluster = $biggest")
                println (s"Biggest cluster indices = $indices")
            } // if

            var i       = 0                                              // element index to reassign
            if (useDistance) {
                i           = clustr.indexOf (biggest)                   // first element in biggest cluster
                var max     = distance (x(i), cent(biggest))             // first distance in biggest cluster
                for (ii <- indices) {                                    // find furthest in biggest cluster
                    val dist = distance (x(ii), cent(biggest))
                    if (dist > max) { max  = dist; i = ii }
                } // for
            } else {
                val ran = new Randi (0, indices.size-1)                  // random integer generator
                i       = indices(ran.igen)                              // randomly pick one point from biggest cluster
            } // if
            sizes(clustr(i)) -= 1                                        // decrement size of previous cluster
            clustr(i)         = c                                        // reassign vector x(i) to cluster c
            sizes(c)         += 1                                        // increment size of cluster c
            if (DEBUG) println (s"New clustering = ${clustr.deep}")
        } // for
    } // fixEmptyClusters

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Check for empty clusters and throw an execption if found.
     */
    private def emptyClusters ()
    {
        for (c <- 0 until k if ! (clustr contains c)) throw new Exception (s"Empty cluster c = $c")
    } // emptyClusters

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reassign each vector/point to the cluster with the closest centroid.
     *  Indicate done, if no points changed clusters (for stopping rule).
     */
    def reassign (): Boolean =
    {
        var done = true                                          // done indicates no changes
        for (i <- x.range1) {
            val v = x(i)                                         // let v be the ith vector
            for (c <- 0 until k) {
                val newDist = distance (v, cent(c))              // calc distance to centroid c
                if (newDist < dist(i)) {                         // is it closer than old distance
                    dist(i)           = newDist                  // make it the new distance
                    sizes(clustr(i)) -= 1                        // decrement size of previous cluster
                    clustr(i)         = c                        // reassign vector x(i) to cluster c
                    sizes(c)         += 1                        // increment size of cluster c
                    done              = false                    // changed clusters => not done
                } // if
            } // for
        } // for
        done                                                     // return whether there were no changes
    } // reassign

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Calculate the centroids based on current assignment of points to clusters.
     */
    def calcCentroids ()
    {
        val cx = new MatrixD (k, x.dim2)                         // to hold sum of vectors for each cluster
        val cs = new VectorD (k)                                 // to hold number of vectors in each cluster
        for (i <- x.range1) {
            val ci  = clustr(i)                                  // x(i) currently assigned to cluster ci
            cx(ci)  = cx(ci) + x(i)                              // add the next vector in cluster
            cs(ci) += 1.0                                        // add 1 to number in cluster
        } // for
        for (c <- 0 until k) cent(c) = cx(c) / cs(c)             // divide to get averages/means
    } // calcCentroids

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Iteratively recompute clusters until the assignment of points does not
     *  change, returning the final cluster assignment vector.
     */
    def cluster (): Array [Int] =
    {
        import scalation.util.gauge 

        tc1 = gauge {
        if (primary) {
            assign ()                                            // randomly assign points to clusters
            fixEmptyClusters (false)                                  // swap points into empty clusters
            calcCentroids ()                                     // calculate the initial centroids
        } else {
            pickCentroids ()                                     // alt., pick points for initial centroids
            fixEmptyClusters (false)                                  // swap points into empty clusters
        } // if

        if (DEBUG) {
            println ("(" + 0 + ") clustr = " + clustr.deep)
            println ("(" + 0 + ") cent   = " + cent)
        } // if

        breakable { for (l <- 1 to MAX_ITER) {
            if (reassign ()) break                               // reassign points to clusters (no change => break)
            fixEmptyClusters ()                                  // check for empty clusters            
            calcCentroids ()                                     // re-calculate the centroids
            if (DEBUG) {
                println ("(" + l + ") clustr = " + clustr.deep)
                println ("(" + l + ") cent   = " + cent)
            } // if
        }} // for
            emptyClusters ()                                         // should not have any empty clusters
        } // time

        tc2 = gauge {
            if (post) randomSwaps ()                                 // swap points to improve sse
        } // time
        clustr                                                   // return the cluster assignment vector
    } // cluster

    @inline
    private def swapPoints (i: Int, j: Int)
    {
        val temp = clustr(i)
        clustr(i) = clustr(j)
        clustr(j) = temp
        calcCentroids ()
    } // swapPoints

    private def randomSwaps ()
    {
        if (DEBUG) banner (s"randomSwaps ()")
        for (i <- 0 until x.dim1-1; j <- i+1 until x.dim1 if clustr(i) != clustr(j)) {
            val ssei1 = sse (clustr(i))
            val ssej1 = sse (clustr(j))
            val sum1  = ssei1 + ssej1
            swapPoints (i, j)
            val ssei2 = sse (clustr(i))
            val ssej2 = sse (clustr(j))
            val sum2  = ssei2 + ssej2            
            if (DEBUG) {
                println (s"ssei1 = $ssei1; ssej1 = $ssej1; sum1 = $sum1")
                println (s"ssei2 = $ssei2; ssej2 = $ssej2; sum2 = $sum2")
            } // if
            if (sum2 > sum1) {
                swapPoints (i, j)
                if (DEBUG) println (s"swapping back")
            } // if 
        }  // for
    } // randomSwaps

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a new point/vector 'y', determine which cluster it belongs to,
     *  i.e., the cluster whose centroid it is closest to.
     *  @param y  the vector to classify
     */
    def classify (y: VectorD): Int =
    {
        var dist = distance (y, cent(0))                         // calc distance to centroid 0
        var clus = 0                                             // assign y to cluster 0
        for (c <- 1 until k) {
            val newDist = distance (y, cent(c))                  // calc distance to centroid c
            if (newDist < dist) {                                // is it closer than old distance
                dist = newDist                                   // make it the new distance
                clus = c                                         // assign y to cluster c
            } // if
        } // for
        clus                                                     // return cluster y belongs to
    } // classify

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the sum of squared errors (distance sqaured from centroid for all points)
     */
    def sse (): Double =
    {
        //x.range1.view.map (i => distance (x(i), cent(clustr(i)))).sum
        var sum = 0.0
        for (i <- x.range1) sum += distance (x(i), cent(clustr(i)))
        sum
    } // sse

    def sse (c: Int): Double =
    {
//        x.range1.view
//                .filter (clustr(_) == c)
//                .map (i => distance (x(i), cent(clustr(i)))).sum
        var sum = 0.0
        for (i <- x.range1) {
            val cli = clustr(i)
            if (cli == c) sum += distance (x(i), cent(cli))
        } // for
        sum
    } // sse

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Check to see if the sum of squared errors is optimum.
     *  @param opt  the known (from human/oracle) optimum
     */
    def checkOpt (opt: Double): Boolean = sse () <= opt + 1 // =~ opt @TODO FIX 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// LETS MOVE THE RANDOM MATRIX AND GAP CALCULATIONS INTO A SUBCLASS

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the sum of pairwise distances for all points in cluster 'r'
     *  @param r  the cluster of interest
     */

    def d(r: Int): Double =
    {
	var sum = 0.0
	for (i <- 0 until x.dim1 if clustr(i)==r;
	     j <- 0 until x.dim1 if clustr(j)==r ) sum += distance( x(i) , x(j) )
	sum
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the pooled within cluster sum of squares around the cluster
    *   means. 
    */

    def w(): Double =
    {
	var sum = 0.0
	for( r <- 0 until k ) sum += d(r) / (2 * clustr.count(_ == r) )
	sum
    }

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the pooled within-cluster sum of squares around the cluster means
     *   of a simulated distribution. 
     */

    def wStar(): Double =
    {
	var max     = 0.0
	var min     = Double.MaxValue
	var uniDist = Array.ofDim[Double](x.dim1,x.dim2) 
	for( j <- 0 until x.dim2)
	{
		for( i <- 0 until x.dim1 )
		{
			if( x.apply(i,j) > max ) max = x.apply(i,j)
			if( x.apply(i,j) < min ) min = x.apply(i,j)
		}
		val rand = new Uniform(min,max,s)
		for( i <- 0 until x.dim1 ) uniDist(i)(j) = rand.gen 
	}//for
	var uniDMatrix = new MatrixD(x.dim1, x.dim2, uniDist)
	val cl = new KMeansClustering(uniDMatrix, k, s, primary)
	cl.w()
    }

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /**Return a tuple (gapK, sk) containing the estimated gap statistic for this
     * clustering in position one and the sk in position two.
     * @param B  the number of uniform distributions used to calculate the gap.
     */

    def gap(B: Int = 10): (Double,Double) =
    {
	import math.log
	import math.pow
	var gapK = 0.0
	var lBar = 0.0
	var sdk  = 0.0
	var sk   = 0.0
	var wStars = Array.ofDim[Double](B)
	for( b <- 0 until B )
	{
		wStars(b)  = wStar()
		gapK      += log( wStars(b) ) - log ( w() )
		lBar      += log( wStars(b) )
	}//for
	lBar /= B
	for( b <- 0 until B) sdk += pow( log( wStars(b) ) - lBar , 2 )
	sdk /= B
	sdk = pow( sdk, 1/2 )
	sk = sdk * pow( 1 + 1/B , 1/2)          //is this correct? The paper notation confused me...
	(gapK / B,  sk)
    }

    private def getMinMax (): (VectorD, VectorD) =
    {
        (VectorD (for (j <- x.range2) yield x.col(j).min()),
         VectorD (for (j <- x.range2) yield x.col(j).max()))
    } // getMinMax

    private def genRand (): MatrixD =
    {
        val xrand = new MatrixD (x.dim2, x.dim1)
        val (min, max) = getMinMax ()
        for (i <- xrand.range1) xrand(i) = RandomVecD (xrand.dim2, max(i), min(i), stream = i).gen 
        xrand.t
    } // genRand

    def simulate (): KMeansClustering =
    {
        val xrand = genRand ()
        println (s"random matrix = $xrand")
        val cl = new KMeansClustering (xrand, k, s)
        val clusters = cl.cluster ()
        println (s"sim clusters = ${clusters.deep}")
        cl 
    } // simulate

} // KMeansClustering class


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClusteringTest` object is used to test the `KMeansClustering` class.
 *  > run-main scalation.analytics.clusterer.KMeansClusteringTest
 */
trait KMeansClusteringTester
{
    import scalation.stat.Statistic
    def test (v: MatrixD, k: Int, primary: Boolean, remote: Boolean, post: Boolean, opt: Double = -1)
    {
        banner (s"test (primary = $primary, remote = $remote, post = $post)")
        val statSSE = new Statistic ()
        val statTC1 = new Statistic ()
        val statTC2 = new Statistic ()
        var ok = 0
        for (s <- 0 until 1000) {                         // test with different random streams
            //banner ("KMeansClustering for stream s = " + s)
            val cl  = new KMeansClustering (v, k, s, primary = primary, remote = remote, post = post)
            cl.cluster ()
            val sse = cl.sse ()            
            //println ("--- final cluster   = " + cl.cluster ().deep)
            //println ("--- final sse       = " + sse)            
            statSSE.tally (sse)
            statTC1.tally (cl.tc1)
            statTC2.tally (cl.tc2)
            if ((opt != -1) && (cl.checkOpt (opt))) ok += 1

            /*
            var sum = 0.0
            for (c <- 0 until k) {
                val csse = cl.sse (c)
                println (s"sse($c) = $csse")
                sum += csse
            } // for
            assert (sum == cl.sse ())
            println (s"sum sse(c) = $sum")
            //println (s"calc sse(c) = ${cl.sse ()}")
             */

        } // for
        if (opt != -1) println (s"ok = $ok")
        println (Statistic.labels)
        println (statSSE)
        println (statTC1)
        println (statTC2)

    } // test
} // KMeansCluteringTester

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClusteringTest` object is used to test the `KMeansClustering` class.
 *  > run-main scalation.analytics.clusterer.KMeansClusteringTest
 */
object KMeansClusteringTest extends App with KMeansClusteringTester
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

    val tf = Array (true, false)
    for (primary <- tf; remote <- tf; post <- tf) {
        test (v, k, primary, remote, post, 3)
    } // for

} // KMeansClusteringTest object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClusteringTest2` object is used to test the `KMeansClustering` class.
 *  > run-main scalation.analytics.clusterer.KMeansClusteringTest
 */
object KMeansClusteringTest2 extends App with KMeansClusteringTester
{

    val v = new MatrixD ((8, 2),  1.0,  1.0,
                                  1.0,  3.0,
                                  5.0, 18.0,
                                  5.0, 20.0,
                                  9.0, 10.0,
                                  9.0, 12.0,
                                  15.0, 30.0,
                                  15.0, 32.0)

    val k = 4

    println ("v = " + v)
    println ("k = " + k)
    println ("----------------------------------------------------")

    val tf = Array (true, false)
    for (primary <- tf; remote <- tf; post <- tf) {
        test (v, k, primary, remote, post, 8)
    } // for

} // KMeansClusteringTest2 object

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClusteringTest2` object is used to test the `KMeansClustering` class.
 *  > run-main scalation.analytics.clusterer.KMeansClusteringTest3
 */
object KMeansClusteringTest3 extends App with KMeansClusteringTester
{

    import scalation.random.{Normal, Bernoulli}

    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 1.0)
    val dist2 = Normal (8.0, 1.0)
    val v    = new MatrixD (50, 2)
    val k    = 4

    for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen)

    import scalation.plot.Plot    
    new Plot (v.col(0), v.col(1))    

    println ("v = " + v)
    println ("k = " + k)
    println ("----------------------------------------------------")

    val tf = Array (true, false)
    for (primary <- tf; remote <- tf; post <- tf) {
        test (v, k, primary, remote, post, 76.6)
    } // for

} // KMeansClusteringTest3 object

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClusteringTest4` object is used to test the `KMeansClustering` class.
 *  > run-main scalation.analytics.clusterer.KMeansClusteringTest4
 */
object KMeansClusteringTest4 extends App with KMeansClusteringTester
{

    import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer
    import scalation.random.{Normal, Bernoulli}

    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 1.0)
    val dist2 = Normal (8.0, 1.0)
    val v    = new MatrixD (50, 2)
    val k    = 4

    for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
                                        if (coin.gen == 0) dist1.gen else dist2.gen)

    import scalation.plot.Plot    
    new Plot (v.col(0), v.col(1))

    val cl = new KMeansPlusPlusClusterer (k)
    
    println ("v = " + v)
    println ("k = " + k)
    println ("----------------------------------------------------")

    val tf = Array (true, false)
    for (primary <- tf; remote <- tf; post <- tf) {
        test (v, k, primary, remote, post, 76.6)
    } // for

} // KMeansClusteringTest4 object
