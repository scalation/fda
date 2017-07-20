

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.3
 *  @date    Sun Mar 12 16:49:17 EDT 2017
 *  @see     LICENSE (MIT style license file).
 *
 *  @see www.jstor.org/stable/3695642?seq=1#page_scan_tab_contents
 */

package scalation.analytics.clusterer

import scala.util.control.Breaks.{breakable, break}
import scala.collection.mutable.{ArrayBuffer, Set}
import scala.math.min

import scalation.linalgebra.{MatrixD, SparseMatrixD}
import scalation.random.RandomVecSample
import scalation.util.SortingI

import Algorithm._

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `TightClusterer` class uses tight clustering to eliminate points that
 *  do not not fit well in any cluster.
 *  @param x     the vectors/points to be clustered stored as rows of a matrix
 *  @param k0    the number of clusters to make
 *  @param kmin  the minimum number of clusters to make
 *  @param s     the random number stream (to vary the clusters made)
 */
//class TightClusterer (x: MatrixD, k0: Int, kmin: Int, s: Int = 0)
class TightClusterer (x: MatrixD, k0: Int, kmin: Int, s: Int = 0, alpha: Double = 0.7, beta: Double = 0.9, q: Int = 7, b: Int = 10, ratio: Double = 0.7, levels:Int = 3)
//      extends Clusterer
{
        // NOTE q = top.can
	//I think seq.num is levels?
	
    private val DEBUG = false                                  // debug flag
    //private val ratio = 0.7                                   // subsampling ratio
    //private val alpha = 0.2                                   // how far below 1 to set threshold
    private val thres = 1 - alpha                             // membership threshold for high scores
    //private val beta  = 0.9                                   // similarity threshold
    //private val b     = 10                                    // number of times to resample
    //private val q     = 7                                     // number of candidates for each k
    private val n     = x.dim1                                // size of whole sample/population
    private val avail = Array.fill(x.dim1)(true)     	      // the not yet tightly clustered data points

    //private val levels   = 3                                          // number of levels to try
    private val clusters = new ArrayBuffer [Set [Int]] ()
    private val topClubs = Array.ofDim [ArrayBuffer [Set [Int]]] (levels)

    private val mda = Array.ofDim [Double] (n, n)
    private val da  = Array.ofDim [Double] (n, n)
    private val ya  = Array.ofDim [Double] (x.dim1, x.dim2)

    private val md = new MatrixD (n, n, mda)                           // mean comembership matrix
    private val d  = new MatrixD (n, n, da)                           // comembership matrix for current sample


/*    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a new reandom subsample.
     */
    def createSubsample (): (MatrixD, Array [Int]) =
    {
	val nn	  = avail.count(_ == true)		      // the number of available rows (i.e. - rows which haven't been tight clustered yet...)
   	val ns    = (nn * ratio).toInt			      // size of a random subsample
	//println(s"From subsamp ns: $ns")
    	val sr    = 0 until ns                                // sample range
	val strm  = (System.currentTimeMillis % 1000).toInt
    	val rsg   = RandomVecSample (nn, ns, strm)	      // random sample generator
	
        val indexMap = rsg.igen ().toArray                    // select e.g. 5th, 3rd, 7th  // FIX - why toArray
	//print(s"indexMap: ${indexMap.deep}")
	val subsamp  = new MatrixD(indexMap.length,x.dim2)    // a matrix to hold the specified vectors from the positions specified by indexMap 
	val arrayMap = avail.zipWithIndex.map{case (e,i) =>
	    	       			  if(e) i else -1}.filterNot(_ == -1) 	// the indices of the rows specified in indexMap e.g. 5th => index 7, 3rd => 3, 7th => 9
        //val subsamp  = x.selectRows (arrayMap)                    	 	// generate random subsample
	//println(s"arrayMap: ${arrayMap.deep}") 
	for( i <- subsamp.range1 ) {
	    //println(s"i: $i")
	    //println(s"indexMap(i): ${indexMap(i)}")
	    //println(s"arrayMap(indexMap(i)): ${arrayMap(indexMap(i))}")
	    subsamp(i) = x(arrayMap(indexMap(i)))    	// fill the subsamp with the rows from x specified by arrayMap e.g. x(5), x(7), x(9)
	}
        //println (s"subsamp = $subsamp")
        (subsamp, indexMap) 
    } // createSubsample
 */

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a new reandom subsample.
     */
    def createSubsample (subsamp: MatrixD): Array [Int] =
    { 
	val nn	  = avail.count(_ == true)		      // the number of available rows (i.e. - rows which haven't been tight clustered yet...)
   	val ns    = (nn * ratio).toInt			      // size of a random subsample
	//println(s"From subsamp ns: $ns")
    	val sr    = 0 until ns                                // sample range
	val strm  = (System.currentTimeMillis % 1000).toInt
    	val rsg   = RandomVecSample (nn, ns, strm)	      // random sample generator
	
        val indexMap = rsg.igen ().toArray                    // select e.g. 5th, 3rd, 7th  // FIX - why toArray
	//print(s"indexMap: ${indexMap.deep}")
	// val subsamp  = new MatrixD(indexMap.length,x.dim2)    // a matrix to hold the specified vectors from the positions specified by indexMap 
	val arrayMap = avail.zipWithIndex.map{case (e,i) =>
	    	       			  if(e) i else -1}.filterNot(_ == -1) 	// the indices of the rows specified in indexMap e.g. 5th => index 7, 3rd => 3, 7th => 9
        //val subsamp  = x.selectRows (arrayMap)                    	 	// generate random subsample
	//println(s"arrayMap: ${arrayMap.deep}") 
	for( i <- subsamp.range1 ) {
	    //println(s"i: $i")
	    //println(s"indexMap(i): ${indexMap(i)}")
	    //println(s"arrayMap(indexMap(i)): ${arrayMap(indexMap(i))}")
	    subsamp(i) = x(arrayMap(indexMap(i)))    	// fill the subsamp with the rows from x specified by arrayMap e.g. x(5), x(7), x(9)
	}
        //println (s"subsamp = $subsamp")
        indexMap 
    } // createSubsample

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Computet the mean comembership matrix by averaging results from several subsamples.
     */
    def computeMeanComembership (k: Int, y: MatrixD): MatrixD =
    {
	val unclustered = avail.count(_ == true)
    	val nn	  = (unclustered * ratio).toInt
	val clustr2 = Array.ofDim[Int](n)			   // to hold the future clustering of our data classified by centroids of some subset sample clustering
                                                                   // val d = new MatrixD (n, n)                           // comembership matrix for current sample
        md.clear ()
	
        for (l <- 0 until b) {
            d.clear ()                                             // clear the comembreship matrix
            y.clear ()
            val imap = createSubsample (y)                // create a new subsample
	    val kmc	  = new KMeansPPClusterer(y,k,s=(s+l)%1000)
            val clustr 	  = kmc.cluster ()                       // get the clusters
	    val cents 	  = kmc.centroids()
            //println (s"clustr = ${clustr.deep}, cents: $cents")
	    for (i <- x.range1 ) clustr2(i) = if( avail(i) ) kmc.classify2(x(i), cents) else -1
	    if (DEBUG) println(s"clustr2: ${clustr2.deep}")

            for (i <- x.range1; j <- x.range1 if (clustr2(i) == clustr2(j) && clustr2(i) >= 0)) {
	    	//println(s"i: $i")
		//println(s"j: $j")
	        d(i, j) = 1.0
	    }
            //println (s"d = $d")
            md += d 
        } // for
        md /= b // ratio * b                                       // compute mean
        md                                                    // return result
    } // computeMeanComembership

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Form candidate clusters by collecting points with high average comembership
     *  scores together in clusters (clubs).
     *  @param md  the mean comembership matrix
     */
    def formCandidateClusters (md: MatrixD): ArrayBuffer [Set [Int]] =
    {
	//println(s"From formCandidateClusters, n: $n")
	// I don't think we should be using the available list here...
        //val avail = Array.fill (n)(true)                      // whether a point is available
        val clubs = new ArrayBuffer [Set [Int]] ()            // list of clubs
        for (i <- 0 until md.dim1 if avail(i)) {
            val club = Set (i)                                // put i in a club
            //avail(i) = false                                  // make i unavailable
            for (j <- i until md.dim1 if ( avail(i) && md(i,j) >= thres )) { club += j}//; avail(j) = false }
            if( club.size > 1 ) clubs += club
        } // for
        clubs
    } // formCandidateClusters

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Order the clubs (candidate clusters) by size, returning the rank order
     *  (largest first).
     *  @param clubs  the candidate clusters
     */
    def orderBySize (clubs: ArrayBuffer [Set [Int]]): Array [Int] =
    {
        val sizes = clubs.map (_.size).toArray                // record sizes of clubs
        new SortingI (sizes).iselsort2 ()                     // indirectly sort by size
    } // orderBySize
/*
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Select candidates for tight clusters in the K-means algorithm for a given
     *  number of clusters 'k'.  This corresponds to Algorithm A in the paper/URL.
     *  @param k  the number of clusters
     */
    def selectCandidateClusters (k: Int): (ArrayBuffer [Set [Int]], Array [Int]) =
    {
        println ("AAAAAAAAAAAAAA")
        val md    = computeMeanComembership (k)               // mean comembership
        println ("BBBBBBBBBBBBBB")        
        val clubs = formCandidateClusters (md)                // form candidate clusters (clubs)
        println ("CCCCCCCCCCCCCC")        
        val order = orderBySize (clubs)                       // determine rank order by club size
        println ("DDDDDDDDDDDDDD")
        if (DEBUG) {
            println (s"mean = $md")
            println (s"clubs = $clubs")
            println (s"order = ${order.deep}")
        } // if
        (clubs, order)
    } // selectCandidateClusters
 */
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Select candidates for tight clusters in the K-means algorithm for a given
     *  number of clusters 'k'.  This corresponds to Algorithm A in the paper/URL.
     *  @param k  the number of clusters
     */
    def selectCandidateClusters (k: Int, y: MatrixD): (ArrayBuffer [Set [Int]], Array [Int]) =
    {
        val md    = computeMeanComembership (k,y)               // mean comembership
        val clubs = formCandidateClusters (md)                // form candidate clusters (clubs)
        val order = orderBySize (clubs)                       // determine rank order by club size
        if (DEBUG) {
            println (s"mean = $md")
            println (s"clubs = $clubs")
            println (s"order = ${order.deep}")
        } // if
        (clubs, order)
    } // selectCandidateClusters


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Pick the top q clubs based on club size.
     *  @param clubs  all the clubs (candidate clusters)
     *  @param order  the rank order (by club size) of all the clubs
     */
    def pickTopQ (clubs: ArrayBuffer [Set [Int]], order: Array [Int]):  ArrayBuffer [Set [Int]] =
    {
        val ml = ArrayBuffer [Set [Int]] ()
        for (i <- 0 until min (q, clubs.size)) ml += clubs(order (i))
        ml
    } // pickTopQ

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the similarity of two clubs as the ratio of the size of their
     *  intersection to their union.
     *  @param c1  the first club
     *  @param c2  the second club
     */
    def sim (c1: Set [Int], c2: Set [Int]): Double = (c1 & c1).size / (c1 union c2).size

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find a the first tight and stable cluster from the top candidate clubs.
     *  To be stable, a club must have a similar club at the next level (next k value).
     *  @param topClubs  the top clubs for each level to be search for stable clusters
     */
    def findStable (topClubs: Array [ArrayBuffer [Set [Int]]]): (Int, Set [Int]) =
    {
        for (lev <- 0 until topClubs.length-1) {
            for (c1 <- topClubs (lev); c2 <- topClubs (lev+1)) {
                if (sim (c1, c2) >= beta) return (lev+1, c2)        // found a stable cluster
            } // for	     	    	  	       		  
        } // for
        return (-1, null)                                         // none found
    } // findStable

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a set of points/vectors, put them in clusters, returning the cluster
     *  assignment vector.  A basic goal is to minimize the sum of the distances
     *  between points within each cluster.
     */
    def cluster (): ArrayBuffer [Set [Int]] =
    {
	// var done    = false
        breakable { for (kc <- k0 to kmin by -1) {                            // iteratively decrement kc (k current value)
	                                                                      //println(s"kc : $kc")
            val nn = (avail.count( _ == true) * ratio).toInt
            if( nn == 1 ) break
//	    if( avail.count( _ == true ) == 1 ) break // done = true
// if( !done ) break
            val y = new MatrixD(nn, x.dim2, ya)
            for (k <- kc until kc + levels) {
                y.clear()
                val (clubs, order) = selectCandidateClusters (k, y)
                topClubs(k-kc) = pickTopQ (clubs, order)
            } // for
            if (DEBUG) println (s"topClubs = ${topClubs.deep}")
            val (lev, stable) = findStable (topClubs)             // next stable cluster
            if (DEBUG) println (s"(lev, stable) = ($lev, $stable)")
            if (lev >= 0) {
                clusters      += stable                           // add to stable clusters
                //topClubs(lev) -= stable                           // remove from top clubs (WHY?)
		for( i <- topClubs ) i.clear()
		for( i <- stable ) avail(i) = false
	        if( avail.count(_ == true ) == 0 ) break // done = true
            } else {
                if (DEBUG) println (s"no stable cluster found for kc = $kc: $stable")
            } // if
        }} // for // breakable
        if (DEBUG) println (s"clusters = $clusters")
        clusters
    } // cluster

} // TightClusterer class

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `TightClustererTest` is used to test the `TightClusterer` class.
 *  > run-main scalation.analytics.clusterer.TightClustererTest
 */
object TightClustererTest extends App
{
    val v = new MatrixD ((7, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 5.0, 4.0,
                                 4.0, 5.0,
                                 9.0, 8.0,
                                 8.0, 9.0,
				 19.0, 32.0)

   val (k0, kmin) = (5,1)
   for (s <- 0 until 5) {
       println(s"\n\n\n//::::::::::::::::::::::::::\nTight Cluster test for s = $s\n//::::::::::::::::::::::::::\n\n\n")
       val tcl = new TightClusterer (v, k0, kmin, s)
       val clust = tcl.cluster ()
       assert(!clust.flatten.contains(6))
   } // for

} // TightClustererTest object

