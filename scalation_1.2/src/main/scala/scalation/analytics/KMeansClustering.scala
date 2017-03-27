
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.2
 *  @date    Tue May 29 14:45:32 EDT 2012
 *  @see     LICENSE (MIT style license file).
 */

package scalation.analytics

import scala.collection.mutable.Set
import scala.util.control.Breaks.{breakable, break}

import scalation.linalgebra.{MatrixD, VectorD}
import scalation.random.{Randi, Uniform}
import scalation.util.{banner, Error}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClustering` class cluster several vectors/points using k-means
 *  clustering.  Either (1) randomly assign points to 'k' clusters or (2) randomly
 *  pick 'k' points as initial centroids (technique (1) to work better and is the
 *  primary technique).  Iteratively, reassign each point to the cluster containing
 *  the closest centroid.  Stop when there are no changes to the clusters.
 *  @param x        the vectors/points to be clustered stored as rows of a matrix
 *  @param k        the number of clusters to make
 *  @param s        the random number stream (to vary the clusters made)
 *  @param primary  true indicates use the primary technique for initiating the clustering
 */
class KMeansClustering (x: MatrixD, k: Int, s: Int = 0, primary: Boolean = true)
      extends Clusterer with Error
{
    if (k >= x.dim1) flaw ("constructor", "k must be less than the number of vectors")

    private val DEBUG    = true                          // debug flag
    private val MAX_ITER = 200                           // the maximum number of iterations
    private val cent     = new MatrixD (k, x.dim2)       // the k centroids of clusters
    private val clustr   = Array.ofDim [Int] (x.dim1)    // assignment of vectors to clusters
    private val dist     = new VectorD (x.dim1)          // distance to closest centroid
    dist.set (Double.PositiveInfinity)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute a distance metric between vectors/points u and v.
     *  @param u  the first vector/point
     *  @param v  the second vector/point
     */
    def distance (u: VectorD, v: VectorD): Double =
    {
        (u - v).normSq       // squared Euclidean norm used for efficiency, may use other norms
    } // distance

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Assign each vector/point to a random cluster.  Primary technique for
     *  initiating the clustering.
     */
    def assign ()
    {
        val ran = new Randi (0, k - 1, s)              // for random integers: 0, ..., k-1
        for (i <- x.range1) clustr(i) = ran.igen       // randomly assign to a cluster
    } // assign

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reassign each vector/point to the cluster with the closest centroid.
     *  Indicate done, if no points changed clusters (for stopping rule).
     */
    def reassign (): Boolean =
    {
        var done = true                                // done indicates no changes
        for (i <- x.range1) {
            val v = x(i)                               // let v be the ith vector
            for (c <- 0 until k) {
                val newDist = distance (v, cent(c))    // calc distance to centroid c
                if (newDist < dist(i)) {               // is it closer than old distance
                    dist(i)  = newDist                 // make it the new distance
                    clustr(i) = c                      // assign vector i to cluster c
                    done = false                       // changed clusters => not done
                } // if
            } // for
        } // for
        done                                           // return whether there were no changes
    } // reassign

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Randomly pick vectors/points to serve as the initial k centroids (cent).
     *  Secondary technique for initiating the clustering.
     */
    def pickCentroids ()
    {
        val ran  = new Randi (0, x.dim1 - 1, s)        // for random integers: 0, ..., x.dim1-1
        val iSet = Set [Int] ()                        // set of integers already generated
        for (c <- 0 until k) {
            var i = ran.igen                           // generate a random integer
            while (iSet contains i) i = ran.igen       // do not allow repeats
            iSet   += i                                // add to set of generated integers
            cent(c) = x(i)                             // let centroid c be the ith vector
        } // for
    } // pickCentroids

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Calculate the centroids based on current assignment of points to clusters.
     */
    def calcCentroids ()
    {
        val cx = new MatrixD (k, x.dim2)               // to hold sum of vectors for each cluster
        val cs = new VectorD (k)                       // to hold number of vectors in each cluster
        for (i <- x.range1) {
            val ci  = clustr(i)                        // x(i) currently assigned to cluster ci
            cx(ci)  = cx(ci) + x(i)                    // add the next vector in cluster
            cs(ci) += 1.0                              // add 1 to number in cluster
        } // for
        for (c <- 0 until k) cent(c) = cx(c) / cs(c)   // divide to get averages/means
    } // calcCentroids

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Iteratively recompute clusters until the assignment of points does not
     *  change, returning the final cluster assignment vector.
     */
    def cluster (): Array [Int] =
    {
        if (primary) {
            assign ()                                  // randomly assign points to clusters
            calcCentroids ()                           // calculate the initial centroids
        } else {
            pickCentroids ()                           // alt., pick points for initial centroids 
        } // if
        println ("(" + 0 + ") clustr = " + clustr.deep)
        println ("(" + 0 + ") cent   = " + cent)

        breakable { for (l <- 1 to MAX_ITER) {
            if (reassign ()) break                     // reassign points to clusters (no change => break)
            calcCentroids ()                           // re-calculate the centroids
            if (DEBUG) {
                println ("(" + l + ") clustr = " + clustr.deep)
                println ("(" + l + ") cent   = " + cent)
            } // if
        }} // for   
        clustr                                         // return the cluster assignment vector
    } // cluster

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a new point/vector 'y', determine which cluster it belongs to,
     *  i.e., the cluster whose centroid it is closest to.
     *  @param y  the vector to classify
     */
    def classify (y: VectorD): Int =
    {
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
    } // classify

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

	
} // KMeansClustering class


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClusteringTest` object is used to test the `KMeansClustering` class.
 *  > run-main scalation.analytics.KMeansClusteringTest
 */
object KMeansClusteringTest extends App
{
    val v = new MatrixD ((6, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 5.0, 4.0,
                                 4.0, 5.0,
                                 9.0, 8.0,
                                 8.0, 9.0)
    val y = VectorD (10.0, 10.0)
    println ("v = " + v)
    println ("y = " + y)
    println ("----------------------------------------------------")

    for (s <- 0 to 4) {                         // test with different random streams
        banner ("KMeansClustering for stream s = " + s)
        val cl = new KMeansClustering (v, 3, s)                 
        println ("--- final cluster = " + cl.cluster ().deep + "\n")
        println ("--- classify " + y + " = " + cl.classify (y) + "\n")
    } // for

} // KMeansClusteringTest object
