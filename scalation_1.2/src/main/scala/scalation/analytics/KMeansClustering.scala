
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.2
 *  @date    Tue May 29 14:45:32 EDT 2012
 *  @see     LICENSE (MIT style license file).
 */

package scalation.analytics

import scala.collection.mutable.Set
import scala.util.control.Breaks.{breakable, break}
import scalation.math.double_exp
import scala.math.sqrt
import scalation.linalgebra.{MatrixD, VectorD, VectorI}
import scalation.random.{Randi, Uniform}
import scalation.util.{banner, Error}
import com.opencsv.CSVWriter
import java.io.{FileWriter,File}
import scala.math.exp

object KMeansClustering
{
    private var streams = VectorI.range (0, 10)
    
    def permuteStreams (stream: Int = 0)
    {
        import scalation.random.PermutedVecI
        streams = PermutedVecI (streams, stream).igen
    } // permuteStreams

}

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
	fixCluster()
    } // assign

    //::
    /**/
    def cSizes(): VectorI =
    {
	val cSizes = new VectorI(k) 
    	for (i <- x.range1) cSizes(clustr(i)) += 1
	cSizes
    }
    
    //::
    /** Returns true if we made any reassignment, false if we did nothing. 
     */
    def fixCluster(): Boolean = 
    {
	var fixed = false
	var sizes = cSizes()
	var min = sizes.min()
	var minPos = sizes.argmin()
	var max = sizes.argmax()
	while( min == 0 )
	{
	   fixed = true
	   breakable{
		for(i <- x.range1 if clustr(i) == max) {
	   	   clustr(i) = minPos
		   break
	   	}
	   } // breakable
	   sizes = cSizes()
	   min = sizes.min()
	   max = sizes.argmax()
	   minPos = sizes.argmin()
	}
	fixed
    }

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reassign each vector/point to the cluster with the closest centroid.
     *  Indicate done, if no points changed clusters (for stopping rule).
     */
    def reassign (): Boolean =
    {
        var done = true                                // done indicates no changes
	var cs = cSizes()
        for (i <- x.range1) {
	    //if( cs( clustr(i) ) >= 2 ){
	    	val v = x(i)				   // let v be the ith vecto
            	breakable{for (c <- 0 until k) {
                    val newDist = distance (v, cent(c))    // calc distance to centroid c
		    println(s"old dist: ${dist(i)}, new dist: $newDist")
                    if (newDist < dist(i)) {               // is it closer than old distance
                       dist(i)  = newDist                 // make it the new distance
                       clustr(i) = c                      // assign vector i to cluster c
                       done = false                       // changed clusters => not done
		       break
                    } // if
            	}} // for
	    //} // if
        } // for
	if( fixCluster() ) done = false 
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
	for (i <- x.range1) dist(i) = distance( x(i) , cent(clustr(i)) )
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
        println (s"(0) clustr = ${clustr.deep}")
        println (s"(0) cent   = $cent")

        breakable { for (l <- 1 to MAX_ITER) {
            if (reassign ()) break                     // reassign points to clusters (no change => break)
            calcCentroids ()                           // re-calculate the centroids
            if (DEBUG) {
                println (s"($l) clustr = ${clustr.deep}")
                println (s"($l) cent   = $cent")
            } // if
        } // for
	}
	
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

    def getCentroids(): MatrixD =
    {
	cent.copy
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the sum of pairwise distances for all points in cluster 'r'
     *  @param r  the cluster of interest
     */

    def d(r: Int): Double =
    {
	var sum = 0.0
	for (i <- 0 until x.dim1 if clustr(i)==r;j <- 0 until x.dim1 if clustr(j)==r )  {
	    sum += distance( x(i) , x(j) )						//aggregate the pairwise distances
	}
	
	sum
    } // d

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the pooled within cluster sum of squares around the cluster
    *   means. 
    */
    def w(): Double =
    {
	var sum = 0.0
	for( r <- 0 until k ){
	     if(clustr.count(_ == r) != 0) sum += d(r) / (2 * clustr.count(_ == r) )	
	} // for
	sum
    } // w

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the pooled within-cluster sum of squares around the cluster means
     *   of a simulated distribution using uniform distribution. 
     */
    def wStar(): Double =
    {
	var max = Double.MinValue
	var min = Double.MaxValue
	
	var uniDist = new MatrixD(x.dim1, x.dim2)
	for( j <- 0 until x.dim2 ){	
		for( i <- 0 until x.dim1 ){
			if( x(i,j) > max ) max = x(i,j)
			if( x(i,j) < min ) min = x(i,j)
		} // for
		//println(s"max: $max, min: $min")
		val rand2 = new Uniform(min,max,(System.currentTimeMillis()%1000).toInt)
		for( i <- 0 until x.dim1 ) uniDist(i, j) = rand2.gen
		max = Double.MinValue;
		min = Double.MaxValue;
	} // for
	
	//println("uniDist: " + uniDist)
	
	val cl = new KMeansClustering(uniDist, k, s, primary)
	cl.cluster()
	cl.w()
    } // wStar

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the pooled within-cluster sum of squares around the cluster means
     *   of a simulated distribution using SVD. 
     */
    def wStar2(): Double =
    {
	import scalation.linalgebra.SVD
	//println(s"x: ${x}")
	var max = Double.MinValue
	var min = Double.MaxValue
	
	val svd = new SVD(x)
	val (u,d,vt) = svd.factor()
	val v = vt.t
	
	var xTransform = x * v
	
	var zTransform = new MatrixD(xTransform.dim1, xTransform.dim2)
	
	for( j <- 0 until xTransform.dim2 ){	
		for( i <- 0 until xTransform.dim1 ){
			if( xTransform(i,j) > max ) max = xTransform(i,j)
			if( xTransform(i,j) < min ) min = xTransform(i,j)
		} // for
		//println(s"max: $max, min: $min")
		val rand2 = new Uniform(min,max,(System.currentTimeMillis()%1000).toInt)
		for( i <- 0 until zTransform.dim1 ) zTransform(i, j) = rand2.gen
		max = Double.MinValue;
		min = Double.MaxValue;
	} // for
	
	//println("uniDist: " + uniDist)

	var z = zTransform * vt
	
	val cl = new KMeansClustering(z, k, s, primary)
	cl.cluster()
	cl.w()
    } // wStar

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /**Return a tuple (gapK, sk) containing the estimated gap statistic for this
     * clustering in position one and the sk in position two.
     * @param B  the number of uniform distributions used to calculate the gap.
     */

    def gap(B: Int = 10): (Double,Double,Double,Double) =
    {
	import math.log
	import math.pow
	var gapK   = 0.0								// the four parts of the tuple
	var lBar   = 0.0
	var sdk    = 0.0
	var sk     = 0.0
	var wStars = Array.ofDim[Double](B)					// an array to hold w* values
	
	for (b <- 0 until B) {
	       	println(s"Monte Carlo simulation iteration b")
		//wStars(b)  = wStar()						// repeat w* b many times for the monte carlo simulation
		wStars(b)  = wStar2() 
		lBar      += log( wStars(b) )					// aggregate the results for later computation
	} // for
	lBar /= B								// this gives us the estimated expected value
	gapK = lBar-log(w()) //log(wVals(k-1)) //log ( wVals(k-1) )		// this gives us the diff btwn the est. exp. value and the observed (i.e. - gap)
	for ( b <- 0 until B) sdk += pow( log( wStars(b) ) - lBar , 2 )		// finding the standard deviation

	sdk /= B     	      	     	       		     	      		
	sdk = pow( sdk, 0.5 )
	sk = pow((1 + 1.0/B) , 0.5)          //is this correct? The paper notation confused me...
	sk *= sdk
	(log(w()),lBar, gapK, sk)
	//(w(),lBar, gapK, sk)
    } // gap

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
    def sse (): Double =
    {
        //x.range1.view.map (i => distance (x(i), cent(clustr(i)))).sum
        var sum = 0.0
        for (i <- x.range1) sum += distance (x(i), cent(clustr(i)))
        sum
    } // sse

} // KMeansClustering class



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClusteringTest` object is used to test the `KMeansClustering` class.
 *  > run-main scalation.analytics.KMeansClusteringTest
 */
object KMeansClusteringTest extends App
{
    val v = new MatrixD ((6, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 4.0, 5.0,
                                 5.0, 4.0,
                                 8.0, 9.0,
                                 9.0, 8.0)
    val y = VectorD (10.0, 10.0)
    println ("v = " + v)
    println ("y = " + y)
    println ("----------------------------------------------------")
    var count = 0
    for (s <- 0 until 100) {                         // test with different random streams
        banner ("KMeansClustering for stream s = " + s)
        val cl = new KMeansClustering (v, 3, s)                 
        println ("--- final cluster = " + cl.cluster ().deep + "\n")
        println ("--- classify " + y + " = " + cl.classify (y) + "\n")
	println (s"w: ${cl.w()}")
	if (cl.w() =~ 3.0 ) count += 1
    } // for
    println(s"Success rate: $count/100.0")
} // KMeansClusteringTest object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `KMeansClusteringTest2` object is used to test the wStar() and Gap()
 *  statistic methods of the `KMeansClustering` class.
 *  > run-main scalation.analytics.KMeansClusteringTest2
 */
object KMeansClusteringTest2 extends App
{
    val v = new MatrixD ((6, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 5.0, 4.0,
                                 4.0, 5.0,
                                 9.0, 8.0,
                                 8.0, 9.0)
    println ("v = " + v)
    println ("----------------------------------------------------")

    for (s <- 0 to 4) {                         // test with different random streams
        banner ("KMeansClustering for stream s = " + s)
        val cl = new KMeansClustering (v, 3, s)                 
	println("Gap: " + cl.gap())	
    } // for


} // KMeansClusteringTest object

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The 'GapStatisticTest' object is used to test the estimate of a gap statistic
 *  using the fully itterated process.
 *  > run-main scalation.analytics.GapStatisticTest
 */

object KMeansGapTest extends App
{
	
	val FILE_NAME  = "chicks.csv"									//or whatever CSV file contains your data
	val file       = new File(FILE_NAME.substring(0,FILE_NAME.indexOf(".")) + "Output.csv")		//for recording the results of the 10 iterations
	val file2      = new File(FILE_NAME.substring(0,FILE_NAME.indexOf(".")) + "GapOut.txt")		//for recording the result of the optimal estimated k
	
	val LENGTH     = 5										//the # of features the data is measured on
	val SIZE       = 578										//the # of data points
	
	val data       = io.Source.fromFile(FILE_NAME);							//the source of the data

	var dataArray  = Array.ofDim[Double](SIZE,LENGTH)						//to create a matrix of the data
	var i          = 0
	for(line <-data.getLines.drop(1)){								//skip the 1st line of the csv file which is usually a header
		 val cols = line.split(",").map(_.trim)							//tokenize each line of the csv file
		 for(j <- 1 to LENGTH-1){								//skip the first column which is the row name
		       cols(j) = cols(j) replaceAll ("\"","")						//can only use floating point numbers
		       cols(j) = cols(j) replaceAll ("c","")						//" "
  		       if(cols(j) != "NA") dataArray(i)(j-1)=cols(j).toDouble				//" " 
		       else dataArray(i)(j-1)=0
		 } // for
		 //println("")
		 i=i+1
	} // 
	val dataMatrix = new MatrixD(dataArray)
	var clus = new KMeansClustering(dataMatrix,3,(System.currentTimeMillis % 1000).toInt)
	print(clus.cluster())
	
}//KMeansGapTest

object ClusGapNewTest extends App
{
	import math.log
	val maxClusters = 10
	import scalation.plot.FPlot
	var kv = VectorD.range(1,maxClusters+1)
	val clusters = Array.ofDim[Array[Int]](maxClusters)
	//val centroids = Array.ofDim[MatrixD](maxClusters)
	println(s"kv: $kv")
	val sse = Array.ofDim[Double](maxClusters)
	//var sse = new VectorD(maxClusters)
	//var usse = new VectorD(maxClusters)
	//val gaps = Array.ofDim[(Double,Double,Double,Double)](maxClusters)
	/*val v = new MatrixD ((6, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 5.0, 4.0,
                                 4.0, 5.0,
                                 9.0, 8.0,
                                 8.0, 9.0)*/

import scalation.random.{Normal, Bernoulli}
    val coin  = Bernoulli ()
    val dist1 = Normal (2.0, 1.0)
    val dist2 = Normal (8.0, 1.0)
    val v     = new MatrixD (50, 2)
    //for (i <- v.range1) v(i) = VectorD (if (coin.gen == 0) dist1.gen else dist2.gen,
          //                              if (coin.gen == 0) dist1.gen else dist2.gen)
					
    for (i <- 0 until 15) v(i) = VectorD (dist1.gen,dist1.gen)
    for (i <- 0 until 15) v(i) = VectorD (dist1.gen,dist2.gen)
    for (i <- 0 until 15) v(i) = VectorD (dist2.gen,dist1.gen)
    for (i <- 0 until 15) v(i) = VectorD (dist2.gen,dist2.gen)
    
	//var wVals = Array.ofDim[Double](maxClusters)
	//for(i <- 1 to maxClusters){//for k = 1...10
	for( i <- 1 to maxClusters){
	       	    println(s"\n\nNew Clustering...\n\n")
	      	    val s = (System.currentTimeMillis % 1000).toInt
	            val cl = new KMeansClustering(v, 4, s)
		    //println(s"new random seed: $s")
		    clusters(i-1) = cl.cluster()
		    //centroids(i-1) = cl.getCentroids()
		    //wVals(i-1) = cl.w()
		    sse(i-1) = cl.sse()
		    //println(s"w() results for $i many clusters: ${log(wVals(i-1))}")
		    //gaps(i-1)  = cl.gap()
		    //println(s"gap results for $i many clusters: ${gaps(i-1)}")
		    //sse(i-1)   = gaps(i-1)._1
		    //usse(i-1)  = gaps(i-1)._2
		    //(w(),lBar, gapK, sk)
		   
	} // for

/*	for (i <- 0 until maxClusters) println(s"${gaps(i)}")

	    breakable{for(i <- 1 until maxClusters){
	        val (lwki, elwki, gapi, ski) = gaps(i)
		val (lwki1, elwki1, gapi1, ski1) = gaps(i+1)
		println(s" gapi: ${gapi}, gapi1: ${gapi1}")
		if(gapi >= gapi1){
	      		 println(s"${i+1} is optimal k")
			 
			 break
	        } // if 
	      
	} // for
	} // breakable

	println(s"usse: $usse")

	println("Clusters: ")*/
	
	//println(s"sse: $sse")
	var percentage = 0.0
	//for(i <- sse.indices ) if (sse(i) <= 76) percentage += 1.0
	//println(s"Accuracy: ${percentage/1000}")
	//for(cluster <- clusters) println(cluster.deep)

	//println("Centroids")
	//for(i <- centroids.indices ; j <- centroids(i).range1 )  println(centroids(i)(j))
	//val plot = new FPlot(kv,sse, kv, (x: Double)=>usse(x.toInt - 1),"Miller fixes Michael's mistakes")
}

object wTester extends App
{
	import math.log
	val s = (System.currentTimeMillis % 1000).toInt
	val v = new MatrixD ((6, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 5.0, 4.0,
                                 4.0, 5.0,
                                 9.0, 8.0,
                                 8.0, 9.0)
	for(i <- 1 to 5) {
	    val cl = new KMeansClustering(v, i, s)
	    cl.cluster()
	    println(s"w(): ${log(cl.w())}")
	}
}

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
        
        
            val statSSE = new Statistic ("sse")
            var ok      = 0
            for (s <- 0 until nstreams) {
                val cl  = new KMeansClustering (v, k, s)
                cl.cluster ()
                val sse = cl.sse ()
                // println (s"stream $s, sse = $sse")
                statSSE.tally (sse)
                if ((opt != -1) && (sse <= opt)) ok += 1
            } // for
            println (Statistic.labels)
            println (statSSE)
            println (s"min sse = ${statSSE.min}")
            if (opt != -1) println (s"optimal = $ok / $nstreams")            
    } // test

    def test2 (v: MatrixD, k: Int, opt: Double = -1, plot: Boolean = false, nstreams: Int = 1000)
    {
        banner (s"Testing KMeansPlusPlusCluster object")
        println (s"v.dim1 = ${v.dim1}")
        println (s"v.dim2 = ${v.dim2}")
        println (s"     k = $k")
        if (plot) new Plot (v.col(0), v.col(1))
        
            
            val statSSE = new Statistic ("sse")
            var ok      = 0
            for (s <- 0 until nstreams) {
                KMeansClustering.permuteStreams (s)
                val cl = new KMeansClustering(v, k, s)
		cl.cluster()
                val sse = cl.sse ()
                // println (s"stream $s, sse = $sse")
                statSSE.tally (sse)
                if ((opt != -1) && (sse <= opt)) ok += 1
            } // for
            println (Statistic.labels)
            println (statSSE)
            println (s"min sse = ${statSSE.min}")
            if (opt != -1) println (s"optimal = $ok / $nstreams")            
        
    } // test2
    

} // KMeansPlusPlusClutererTester



object newTester extends App with KMeansPlusPlusClustererTester
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

}