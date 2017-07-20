package apps

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest6` is used to test the `Smoothing_F` class.
 *  > run-main apps.GeneAnalysis
 */
object GeneAnalysis extends App
{
    import scalation.plot.{PlotM,FPlot}
    import scalation.analytics.clusterer.{GapStatistic, KMeansPPClusterer, TightClusterer}
    import scalation.linalgebra.{VectorD,MatrixD}
    import scalation.util.banner
    import scala.collection.mutable.{ArrayBuffer,Set}
    import scalation.relalgebra.Relation
    import scalation.analytics.fda.{Smoothing_F,SmoothingMethod}
    import scalation.math.FunctionS2S
    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* SET PARMETERS */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    val inFileName  = "Drosophila.csv"		// the name of the file with the data points
    val outFileName = "DrosophilaClusters"  	// the name of the file to write output to, ex. name: "outFileName_LC_COEFFS.csv"
    val gap    	    = 1				// to set the number of knots in the smoother (gap==1 => n many knots, gap==2 => n/2, etc.)
    val ord   	    = 4				// the order for the B-Splines (order 4 => degree 3 polynomial B-Splines)
    val rowSum	    = 100			// the filter threshold for important data (i.e. - minimum value for the sum of a row of a data)
    val kMin 	    = 1				// the minimum number of clusters to find
    val kMax 	    = 6				// the maximum number of clusters to find
    val kEst        = 6				// an estimated number of clusters to make when performing loose clustering
    val useSVD	      	    = true 		// use singlular value decomposition during the GapStatistic process?
    val clusterRawData      = true		// should we cluster the raw data (i.e. - unsmoothed) ? 
    val clusterSmoothedData = true		// should we cluster the smoothed data? 
    val clusterCoefficients = true		// should we cluster the coefficients
    val clusterLoose	    = true		// should we cluster the data without resorting to tight clustering?
    val clusterTight  	    = true		// should we tight cluster the data? 
    val clusterGap	    = false		// should we cluster the data without resorting to tight clustering but using the Gap Statistic to find opt k?
    val plots = true
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* LOAD, FILTER AND PREPARE DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    banner ("Loading observed data...")
    
    val taggedDataBuf  = ArrayBuffer.empty[       String]	// an array buffer for the names of the genes 
    val cleanedDataBuf = ArrayBuffer.empty[Array[Double]]	// an array buffer for the doubles parsed from the strings read in from the CSV 
    val source 	       = io.Source.fromFile(inFileName)		// the source of the data 
    for( line <- source.getLines.drop(1) ){
    	 val cols  = line.split(",").map(_.trim)
	 val clean = (for ( i <- 1 until cols.length ) yield cols(i).toDouble )
	 if ( clean.sum > rowSum ) {
	    cleanedDataBuf += clean.toArray
 	    taggedDataBuf  += cols(0)
	 } // if
    }
    
    val nzdata = new MatrixD(cleanedDataBuf.toArray)
    
    //var nzdata = MatrixD (for (i <- data.range1 if data(i).sum > rowSum) yield data(i), false) // filter out the unimportant data

    //println (s" -    data.min = ${data.min()}; data.max = ${data.max()}")
    //println (s" -   data.dim1 = ${data.dim1}")
    //println (s" - nzdata.dim1 = ${nzdata.dim1}")
    
    val t      = VectorD.range (0, nzdata.dim2)   
    if (plots) new PlotM (t, nzdata, null, s"OBSERVED", lines = true)
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM SMOOTHING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    banner ("Smoothing...")
    
    val smooth = new MatrixD (nzdata.dim1, nzdata.dim2)
    val coeffs = new MatrixD (nzdata.dim1, nzdata.dim2/gap + ord - 2)    // matrix for smoothing spline coefficients
    val funcs  = Array.ofDim [FunctionS2S] (nzdata.dim1)
    for (i <- nzdata.range1) {
        val y     = nzdata(i)
        val moo   = new Smoothing_F (y, t, null, ord, SmoothingMethod.RIDGE, -1, gap) // smoother
        coeffs(i) = moo.train()
        smooth(i) = moo.predict (t)                     // predict for all time points
        funcs(i)  = moo.predict
        //new Plot (t, y, z, s"Smoothing_F ord = $ord; λ = ${moo.getLambda}", lines = true)
    } // for
    
    val t0 = 0.0
    val tn = nzdata.dim2.asInstanceOf[Double]
    if (plots) new FPlot (t0 to tn by 0.1, funcs.toSeq, s"SmoothedData", lines = true)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM CLUSTERING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    /*	Take a matrix of data, cluster it into a specified number of clusters, and plot out the clusters.
     *  @param  data  - the data to cluster
     *	@param  k     - the number of intended clusters
     *	@param  label - a label describing the data/clustering circumstances
     *	@param  orig  - the original data 
     */
    def clusterKMPP (data: MatrixD, k: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): Array[Int] =
    {
        banner (s"K-Means++ Clustering: $label")
        println (s"k = $k")
        val cl    = new KMeansPPClusterer (data, k)
        val cls   = cl.cluster ()
        val sse   = cl.sse ()
        val rs    = cl.rSquared (data)
        println (s"SSE = $sse; rSquared = $rs")
	for (c <- 0 until k) {
            val fclust = for (i <- 0 until cls.size if cls(i) == c) yield forig(i)
            new FPlot (t0 to tn by 0.1, fclust.toSeq, s"KM++: $label; c = $c; n = ${fclust.size}", true)
        } // for
        cls
    } // clusterKMPP

    /*	Take a matrix of data, cluster it into a specified number of clusters using the gap statistic to estimate opt k, and plot out the clusters.
     *  @param  data  - the data to cluster
     *	@param  kMax     - the number of intended clusters
     *	@param  label - a label describing the data/clustering circumstances
     *	@param  orig  - the original data 
     */
    def clusterGap (data: MatrixD, kMax: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): Array[Int] =
    {
        banner (s"Gap Clustering: $label")
	val (cl,cls,k) = GapStatistic.kMeansPP (data,kMax,useSVD = useSVD, plot = false)
        val sse   = cl.sse ()
        val rs    = cl.rSquared (data)
        println (s"SSE = $sse; rSquared = $rs")
	for (c <- 0 until k) {
            val fclust = for (i <- 0 until cls.size if cls(i) == c) yield forig(i)
            new FPlot (t0 to tn by 0.1, fclust.toSeq, s"Gap: $label; c = $c; n = ${fclust.size}", true)
        } // for
        cls
    } // clusterGap

    def clusterTight (data: MatrixD, kMin: Int, kMax: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): ArrayBuffer[Set[Int]] =
    {
        banner (s"Tight Clustering: $label")
        println (s"kMin = $kMin; kMax = $kMax")
	val cl    = new TightClusterer (data, kMax, kMin)
        val cls   = cl.cluster ()
        for (set <- cls) {
            val fclust = for (i <- set) yield forig(i)
            new FPlot (t0 to tn by 0.1, fclust.toSeq, s"TC: $label; n = ${fclust.size}", true)
        } // for
        cls
    } // clusterTight

    def writeOutCluster(cls: Array[Int], label: String) =
    {
    	val vecs = ArrayBuffer.empty[Vector[Any]]
    	for( i <- 0 until cls.size ){
	     vecs += Vector[Any](taggedDataBuf(i), cls(i)) ++ cleanedDataBuf(i).toVector
	}
	val clustRel = Relation( s"genes_${label}",
		    	     	 Seq("gene_name",
				     "cluster_assignment",
				     "em0-2hr"  ,"em2-4hr"  ,"em4-6hr"  ,"em6-8hr"  ,"em8-10hr" ,"em10-12hr",
				     "em12-14hr","em14-16hr","em16-18hr","em18-20hr","em20-22hr","em22-24hr"
				    ),
				 vecs.toSeq,
				 0,
				 "SDDDDDDDDDDDD"
			       ) 
	clustRel.writeCSV(s"${outFileName}_${label}.csv")
    } // writeOutCluster

    def writeOutTightCluster(cls: ArrayBuffer[Set[Int]], label: String) =
    {
    	val vecs = ArrayBuffer.empty[Vector[Any]]
	val seen = ArrayBuffer.empty[Int]
	var c = -1
    	for (set <- cls) { 
	    	c += 1
		for( i <- set ) {
		     vecs += Vector[Any](taggedDataBuf(i) , c) ++ cleanedDataBuf(i).toVector
		     seen += i
		} // for
        } // for 
	for( i <- cleanedDataBuf.toArray.indices if !seen.contains(i) ) vecs += Vector[Any](taggedDataBuf(i), -1) ++ cleanedDataBuf(i).toVector  
	val clustRel = Relation( s"genes_${label}",
		    	     	 Seq("gene_name",
				     "cluster_assignment",
				     "em0-2hr"  ,"em2-4hr"  ,"em4-6hr"  ,"em6-8hr"  ,"em8-10hr" ,"em10-12hr",
				     "em12-14hr","em14-16hr","em16-18hr","em18-20hr","em20-22hr","em22-24hr"),
			         vecs.toSeq,
			     	 0,
			     	 "SDDDDDDDDDDDD")
	clustRel.writeCSV(s"${outFileName}_${label}.csv") 
    } // writeOutTightCluster
    
    val clusts = Array(clusterRawData,clusterSmoothedData,clusterCoefficients)
    val dats   = Array(nzdata,smooth,coeffs)
    val labls  = Array("OBSERVED","SMOOTHED","COEFFS")
  
    for ( i <- 0 to 2 ){
    	if ( clusts(i) ) {
	   if ( clusterGap ){
	      println(s"Calculating gap statistic for ${labls(i)} data...")
	      val cls = clusterGap (dats(i), kMax, labls(i), smooth, funcs)
	      writeOutCluster(cls,s"GSC_${labls(i)}")
	      println(s"done.")
	   }
	   if ( clusterLoose ){
	   	println(s" loose clustering ${labls(i)} data...")
	   	val cls2 = clusterKMPP (dats(i), kEst, labls(i), smooth, funcs)
	   	writeOutCluster(cls2, s"LC_${labls(i)}")
	  	println(s"Done.")
	   }
	   if ( clusterTight ){
	      	println(s"Tight clustering ${labls(i)} data.")
	   	val cls3 = clusterTight (dats(i), kMin, kMax, labls(i), smooth, funcs)
	   	writeOutTightCluster(cls3, s"TC_${labls(i)}")
	   	println(s"Done.")
	   } // if 
	} // if clusts(i)
    } // for i

} // Smoothing_FTest6


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `GeneAnalysis2` is used to test the `Smoothing_F` class.
 *  > run-main apps.GeneAnalysis2
 */
object GeneAnalysis2 extends App
{
    import scalation.plot.{PlotM,FPlot}
    import scalation.analytics.clusterer.{GapStatistic, KMeansPPClusterer, TightClusterer}
    import scalation.linalgebra.{VectorD,MatrixD}
    import scalation.util.banner
    import scala.collection.mutable.{ArrayBuffer,Set}
    import scalation.relalgebra.Relation
    import scalation.analytics.fda.{Smoothing_F,SmoothingMethod}
    import scalation.math.FunctionS2S

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* SET PARMETERS */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

/*
    val inFileName  = if (args(0) == null) "Drosophila.csv"     else args(0) // the name of the file with the data points
    val outFileName = if (args(1) == null) "DrosophilaClusters" else args(1) // the name of the file to write output to, ex. name: "outFileName_LC_COEFFS.csv"
*/
    val inFileName    = args(0)
    val outFileName   = args(1)
    val rowSum        = args(2).toInt                   // the filter threshold for important data (i.e. - minimum value for the sum of a row of a data)
    val kMax 	      = args(3).toInt
    val useSVD	      = if (args(4).toInt>0) true else false
    val plots	      = if (args(5).toInt>0) true else false
    val kEst	      = args(6).toInt
    val alpha	      = args(7).toDouble
    val beta	      = args(8).toDouble
    val q	      = args(9).toInt
    //val t             = args(10)
    
    val clusterRawData       = if (args(11).toInt>0) true else false
    val clusterSmoothedData  = if (args(12).toInt>0) true else false
    val clusterCoefficients  = if (args(13).toInt>0) true else false 
    val clusterLoose         = if (args(14).toInt>0) true else false 
    val clusterTight         = if (args(15).toInt>0) true else false 
    val clusterGap           = if (args(16).toInt>0) true else false 
    val gap    	    = 1				// to set the number of knots in the smoother (gap==1 => n many knots, gap==2 => n/2, etc.)
    val ord   	    = 4				// the order for the B-Splines (order 4 => degree 3 polynomial B-Splines)
    
/*  val rowSum	    = if (args(2) == null) 100  // the filter threshold for important data (i.e. - minimum value for the sum of a row of a data)
    		      else args(2).toInt
*/

	
    val kMin 	    = 1				// the minimum number of clusters to find
/*    
    val kMax 	    = if (args(3) == null) 6
    		      else args(3).toInt
						// the maximum number of clusters to find
    val kEst        = 6				// an estimated number of clusters to make when performing loose clustering
    val useSVD	      	    = true 		// use singlular value decomposition during the GapStatistic process?
    val clusterRawData      = true		// should we cluster the raw data (i.e. - unsmoothed) ? 
    val clusterSmoothedData = true		// should we cluster the smoothed data? 
    val clusterCoefficients = true		// should we cluster the coefficients
    val clusterLoose	    = true		// should we cluster the data without resorting to tight clustering?
    val clusterTight  	    = true		// should we tight cluster the data? 
    val clusterGap	    = false		// should we cluster the data without resorting to tight clustering but using the Gap Statistic to find opt k?
*/

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* LOAD, FILTER AND PREPARE DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    banner ("Loading observed data...")
    
    val taggedDataBuf  = ArrayBuffer.empty[       String]	// an array buffer for the names of the genes 
    val cleanedDataBuf = ArrayBuffer.empty[Array[Double]]	// an array buffer for the doubles parsed from the strings read in from the CSV 
    val source 	       = io.Source.fromFile(inFileName)		// the source of the data
    var m=0
    for( line <- source.getLines.drop(1) ){
    	 val cols  = line.split(",").map(_.trim)
	 val clean = (for ( i <- 1 until cols.length ) yield cols(i).toDouble )
	 if ( clean.sum > rowSum ) {
	    cleanedDataBuf += clean.toArray
	    /*val tag = if ( !(cols(0).equals("") || cols(0).equals("val"))) cols(0) else s"$m"
 	    taggedDataBuf  += tag
	    */
	    taggedDataBuf += s"$m"
	    m+=1
	 } // if
    }
    
    val nzdata = new MatrixD(cleanedDataBuf.toArray)
    
    //var nzdata = MatrixD (for (i <- data.range1 if data(i).sum > rowSum) yield data(i), false) // filter out the unimportant data

    //println (s" -    data.min = ${data.min()}; data.max = ${data.max()}")
    //println (s" -   data.dim1 = ${data.dim1}")
    //println (s" - nzdata.dim1 = ${nzdata.dim1}")
    
    val t      = VectorD.range (0, nzdata.dim2)   
    if (plots) new PlotM (t, nzdata, null, s"OBSERVED", lines = true)
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM SMOOTHING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    banner ("Smoothing...")
    
    val smooth = new MatrixD (nzdata.dim1, nzdata.dim2)
    val coeffs = new MatrixD (nzdata.dim1, nzdata.dim2/gap + ord - 2)    // matrix for smoothing spline coefficients
    val funcs  = Array.ofDim [FunctionS2S] (nzdata.dim1)
    for (i <- nzdata.range1) {
        val y     = nzdata(i)
        val moo   = new Smoothing_F (y, t, null, ord, SmoothingMethod.RIDGE, -1, gap) // smoother
        coeffs(i) = moo.train()
        smooth(i) = moo.predict (t)                     // predict for all time points
        funcs(i)  = moo.predict
        //new Plot (t, y, z, s"Smoothing_F ord = $ord; λ = ${moo.getLambda}", lines = true)
    } // for
    
    val t0 = 0.0
    val tn = nzdata.dim2.asInstanceOf[Double]
    if (plots) new FPlot (t0 to tn by 0.1, funcs.toSeq, s"SmoothedData", lines = true)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM CLUSTERING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    /*	Take a matrix of data, cluster it into a specified number of clusters, and plot out the clusters.
     *  @param  data  - the data to cluster
     *	@param  k     - the number of intended clusters
     *	@param  label - a label describing the data/clustering circumstances
     *	@param  orig  - the original data 
     */
    def clusterKMPP (data: MatrixD, k: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): Array[Int] =
    {
        banner (s"K-Means++ Clustering: $label")
        println (s"k = $k")
        val cl    = new KMeansPPClusterer (data, k)
        val cls   = cl.cluster ()
        val sse   = cl.sse ()
        val rs    = cl.rSquared (data)
        println (s"SSE = $sse; rSquared = $rs")
	//for (c <- 0 until k) {
            //val fclust = for (i <- 0 until cls.size if cls(i) == c) yield forig(i)
            // new FPlot (t0 to tn by 0.1, fclust.toSeq, s"KM++: $label; c = $c; n = ${fclust.size}", true)
        //} // for
        cls
    } // clusterKMPP

    /*	Take a matrix of data, cluster it into a specified number of clusters using the gap statistic to estimate opt k, and plot out the clusters.
     *  @param  data  - the data to cluster
     *	@param  kMax     - the number of intended clusters
     *	@param  label - a label describing the data/clustering circumstances
     *	@param  orig  - the original data 
     */
    def clusterGap (data: MatrixD, kMax: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): Array[Int] =
    {
        banner (s"Gap Clustering: $label")
	val (cl,cls,k) = GapStatistic.kMeansPP (data,kMax,useSVD = useSVD, plot = false)
        val sse   = cl.sse ()
        val rs    = cl.rSquared (data)
        println (s"SSE = $sse; rSquared = $rs")
	for (c <- 0 until k) {
            val fclust = for (i <- 0 until cls.size if cls(i) == c) yield forig(i)
            // new FPlot (t0 to tn by 0.1, fclust.toSeq, s"Gap: $label; c = $c; n = ${fclust.size}", true)
        } // for
        cls
    } // clusterGap

    def clusterTight (data: MatrixD, kMin: Int, kMax: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): ArrayBuffer[Set[Int]] =
    {
        banner (s"Tight Clustering: $label")
        println (s"kMin = $kMin; kMax = $kMax")
	val cl    = new TightClusterer (data, kMax, kMin)
        val cls   = cl.cluster ()
        for (set <- cls) {
            val fclust = for (i <- set) yield forig(i)
            // new FPlot (t0 to tn by 0.1, fclust.toSeq, s"TC: $label; n = ${fclust.size}", true)
        } // for
        cls
    } // clusterTight

    def writeOutCluster(cls: Array[Int], label: String) =
    {
    	val vecs = ArrayBuffer.empty[Vector[Any]]
    	for( i <- 0 until cls.size ){
	     vecs += Vector[Any](taggedDataBuf(i), cls(i)) ++ cleanedDataBuf(i).toVector
	}
	val clustRel = Relation( s"genes_${label}",
		    	     	 Seq("gene_name",
				     "cluster_assignment",
				     "em0-2hr"  ,"em2-4hr"  ,"em4-6hr"  ,"em6-8hr"  ,"em8-10hr" ,"em10-12hr",
				     "em12-14hr","em14-16hr","em16-18hr","em18-20hr","em20-22hr","em22-24hr"
				    ),
				 vecs.toSeq,
				 0,
				 "SDDDDDDDDDDDD"
			       ) 
	clustRel.writeCSV(s"${outFileName}_${label}.csv")
    } // writeOutCluster

    def writeOutTightCluster(cls: ArrayBuffer[Set[Int]], label: String) =
    {
    	val vecs = ArrayBuffer.empty[Vector[Any]]
	val seen = ArrayBuffer.empty[Int]
	var c = -1
    	for (set <- cls) { 
	    	c += 1
		for( i <- set ) {
		     vecs += Vector[Any](taggedDataBuf(i) , c) ++ cleanedDataBuf(i).toVector
		     seen += i
		} // for
        } // for 
	for( i <- cleanedDataBuf.toArray.indices if !seen.contains(i) ) vecs += Vector[Any](taggedDataBuf(i), -1) ++ cleanedDataBuf(i).toVector  
	val clustRel = Relation( s"genes_${label}",
		    	     	 Seq("gene_name",
				     "cluster_assignment",
				     "em0-2hr"  ,"em2-4hr"  ,"em4-6hr"  ,"em6-8hr"  ,"em8-10hr" ,"em10-12hr",
				     "em12-14hr","em14-16hr","em16-18hr","em18-20hr","em20-22hr","em22-24hr"),
			         vecs.toSeq,
			     	 0,
			     	 "SDDDDDDDDDDDD")
	clustRel.writeCSV(s"${outFileName}_${label}.csv") 
    } // writeOutTightCluster
    
    val clusts = Array(clusterRawData,clusterSmoothedData,clusterCoefficients)
    val dats   = Array(nzdata,smooth,coeffs)
    val labls  = Array("OBSERVED","SMOOTHED","COEFFS")
  
    for ( i <- 0 to 2 ){
    	if ( clusts(i) ) {
	   if ( clusterGap ){
	      println(s"Calculating gap statistic for ${labls(i)} data...")
	      val cls = clusterGap (dats(i), kMax, labls(i), smooth, funcs)
	      writeOutCluster(cls,s"GSC_${labls(i)}")
	      println(s"done.")
	   }
	   if ( clusterLoose ){
	   	println(s" loose clustering ${labls(i)} data...")
	   	val cls2 = clusterKMPP (dats(i), kEst, labls(i), smooth, funcs)
	   	writeOutCluster(cls2, s"LC_${labls(i)}")
	  	println(s"Done.")
	   }
	   if ( clusterTight ){
	      	println(s"Tight clustering ${labls(i)} data.")
	   	val cls3 = clusterTight (dats(i), kMin, kMax, labls(i), smooth, funcs)
	   	writeOutTightCluster(cls3, s"TC_${labls(i)}")
	   	println(s"Done.")
	   } // if 
	} // if clusts(i)
    } // for i


    println(s"GeneAnalysis application finished. Close any graph to terminate the application.")
} // GeneAnalysisApp2

