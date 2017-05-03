
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell, John Miller, Dong-Yu Yu
 *  @version 1.3
 *  @date    Fri Mar 24 20:24:40 2017
 *  @see     LICENSE (MIT style license file).
 *
 *  @see open.uct.ac.za/bitstream/item/16664/thesis_sci_2015_essomba_rene_franck.pdf?sequence=1
 *  @see www.jstatsoft.org/article/view/v051i04/v51i04.pdf
 *  @see Functional Data Analysis, Second Edition, Chapter 4
 *  @see http://link.springer.com/book/10.1007%2Fb98888
 */

package scalation.analytics.fda

import scalation.analytics.RidgeRegression
import scalation.linalgebra.{MatrixD, VectoD, VectorD}
import scalation.math.FunctionS2S
import scalation.plot.{FPlot, Plot}
import scalation.util.Error

object SmoothingMethod extends Enumeration
{
    type SmoothingMethod = Value
    val ROUGHNESS, RIDGE, WLS, OLS = Value
} // SmoothingMethod

import SmoothingMethod._

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_F` class fits a time-dependent data vector 'y' to B-Splines.
 *  <p>
 *      y(t(i)) = x(t(i)) + ε(t(i))
 *      x(t) = cΦ(t)
 *  <p>
 *  where 'x' is the signal, 'ε' is the noise, 'c' is a coefficient vector and
 *  'Φ(t)' is a vector of basis functions. 
 *-----------------------------------------------------------------------------
 *  @param y       the (raw) data points/vector
 *  @param t       the data time points/vector
 *  @param τ       the time points/vector for the knots
 *  @param ord     the order (degree+1) of B-Splines (2, 3, 4, 5 or 6)
 *  @param method  the smoothing method
 *  @param lambda  the regularization parameter (>= 0 or -1 to use GCV) 
 */
class Smoothing_F (y: VectorD, t: VectorD, private var τ: VectorD = null, ord: Int = 4, method: SmoothingMethod = ROUGHNESS, lambda: Double = 1E-5, gap: Int = 1)
      extends Error
{
    private val DEBUG = false                   // debug flag
    private val GAP   = gap  // 1                      // gap between time points and knots
    private val m     = t.dim                  // number of data time points
    if (τ == null) τ  = makeKnots
    private val n     = τ.dim                  // number of time points for the knots
    private val bs    = new B_Spline (τ, ord)  // use B-Spline basis functions
    private val ns    = bs.size ()
    private val Φ     = bs.phi ()(t)           // matrix: jth spline at time ti
    private val Σ     = bs.penalty ()(t)       // penalty matrix
    private val I     = MatrixD.eye(ns)         // identity matrix
    private val W     = MatrixD.eye(m) // calcCov (y).inverse    //MatrixD.eye(m)         // weight matrix

    private var λopt  = lambda                 // regularization parameter    
    private var c: VectoD = null               // coefficient vector
    private var e: VectoD = null               // residual/error vector    
    private var sse   = 0.0                    // sum of squared error

    def size = ns

    def summary
    {
        val cov    = calcCov(y)
        val covInv = cov.inverse
        println (s"calcVar     = $cov")
        println (s"1/calcVar   = $cov")
        println (s"cov * covInc = ${cov * covInv}")
        println (s"ns          = $ns")
        println (s"size of λΣ  = ${λopt * Σ.norm}")
        println (s"size of Φ'Φ = ${(Φ.t * Φ).norm}")
        println (s"ratio       = ${(λopt * Σ.norm) / (Φ.t * Φ).norm}")
    } // summary

    def getLambda = λopt

    if (y.dim != m) flaw ("constructor", "require # data points == # data time points")
    if (n > m)      flaw ("constructor", "require # knot points <= # data time points")

    private def makeKnots: VectorD = VectorD.range (0, m/GAP) / ((m-1)/GAP) * t(t.dim-1)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the partial fit for the model using one of several training 
     *  methods. If you multitply the training data by this `MatrixD`, then you
     *  get the symmetric "hat" matrix. If you multitply this `MatrixD` by the
     *  response vector, then you will get the vector estimated model
     *  coefficients.
     *  @param λ  the regularization parameter
     */
    private def pfit (λ: Double): MatrixD = method match {
        case ROUGHNESS => ((Φ.t * W * Φ) + (Σ * λ)).inverse * Φ.t * W
        case RIDGE     => ((Φ.t * W * Φ) + (I * λ)).inverse * Φ.t * W
        case WLS       => (Φ.t * W * Φ).inverse * Φ.t * W
        case OLS       => (Φ.t * Φ).inverse * Φ.t
    } // pfit


    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Calculate the correlation matrix for the basis functions.
     */
    def calcCov (yy: VectorD, k: Int = 1): MatrixD =
    {
        import scalation.stat.vectorD2StatVector
        val avar = yy.variance
        val cov  = MatrixD.eye (yy.dim) * avar
        val acov = yy.acov (k)
        for (i <- 0 until yy.dim - 1) {
            cov(i, i+1) = acov
            cov(i+1, i) = acov
        } // for
        cov
    } // calcCov
    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Computes the "hat" matrix.
     *  @param λ  the regularization parameter
     */
    private def H (λ: Double) : MatrixD = Φ * pfit (λ)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the model, i.e., determine the optimal coeifficient 'c' for the
     *  basis functions by finding optimal Lamdba to minimize gcv.
     */
    def train (): VectoD =
    {
        if (λopt == -1) useGCV ()
        c   = pfit (λopt) * y      //
        e   = y - Φ * c
        sse = e dot e
        c
    } // train

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the degrees of freedom based on the trace of hat matrix.
     *  @param λ  the regularization parameter
     */
    private def df (λ: Double): Double = H (λ).trace

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the Generalized Cross Validation (GCV) score.
     *  @param λ  the regularization parameter
     */
    private def gcv (λ: Double): Double =
    {
        val obj = (m / (m - df (λ))) * (sse / (m - df (λ)))
        //println (s"λ = $λ; ns = $ns; m = $m; df = ${df(λ)}; sse = $sse; obj = $obj")
        // (ns / (ns - df (λ))) * (sse / (ns - df (λ)))
        // (m / (m - df (λ))) * (sse / (m - df (λ)))
        obj 
    } // gcv

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Use the Generalized Cross Validation (GCV) method to find the optimal
     *  λ value. It finds the λ that minimizes the GCV score. 
     */
    private def useGCV ()
    {
        import scalation.minima.GoldenSectionLS
        def f (l: Double): Double =
        {
            c   = pfit (l) * y
            e   = y - Φ * c
            sse = e dot e
            gcv (l)
        } // f
        val gs   = new GoldenSectionLS (f)
        val step = 10
        λopt     = gs.search (step)
        //if (DEBUG) println (s"λopt = $λopt")
    } // useGCV

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Predict the y-value at time point 'tt'.
     *  @param tt  the given time point
     */
    def predict (tt: Double): Double =
    {
        var sum = 0.0
        for (j <- bs.range ()) sum += c(j) * bs (ord) (j, tt)
        sum
    } // predict

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Predict the y-values at all time points in vector 'tv'.
     *  @param tv  the given vector of time points
     */
    def predict (tv: VectorD): VectorD = tv.map (predict (_))

    def plotBasis (tt: VectorD = t) = bs.plot () (tt)

} // Smoothing_F class

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest
 */
object Smoothing_FTest extends App
{
    import scalation.random.Normal
    import math.pow

    val normal1 = Normal (0, 100.0)                                   // normal random variate generator
    val normal2 = Normal (0, 5000.0)                                   // normal random variate generator
    val normal3 = Normal (0, 100.0)                                   // normal random variate generator    

    val t      = VectorD.range (0, 100) / 17.00                      // time points
    val tt     = VectorD.range (0, 1000) / 170.00                    // time points    
/*
    val t      = VectorD.range (0, 50) / (17.0 / 2.0)                      // time points
    val tt     = VectorD.range (0, 500) / (170.00 / 2.0)                    // time points    
 */  
    val y      = t.map ((x: Double) => pow(x-4, 5) + 5.0 * pow(x-4, 4) - 20.0 * pow(x-4, 2) + 4.0 * (x-4))
    val z      = y.map ((e: Double) => e + (if (e < 2) normal1.gen else if (e < 3) normal2.gen else normal3.gen))
    val mMin   = 4 // 2                                              // minimum order to try
    val mMax   = 4                                                   // maximum order to try
    val method = SmoothingMethod.ROUGHNESS // RIDGE                  // smoothing method
    val lambda = -1 // 0.07318687795319875

    new Plot (t, y, z, s"TRUE DATA vs. WITH NOISE", lines = true)

    for (ord <- mMin to mMax) {
        val τ   = null                                               // let `Smoothing_F` nake the knots
        val moo = new Smoothing_F (z, t, τ, ord, method, lambda)     // smoother (use GCV)
        moo.plotBasis (tt)
        val c   = moo.train ()                                       // train -> set coefficients
        val x   = moo.predict (t)                                    // predict for all time points
        new Plot (t, z, x, s"B-Spline Fit: ord = $ord; λ = ${moo.getLambda}", lines = true)
        new Plot (t, y, x, s"TRUTH vs. SMOOTHED", lines = true)
        moo.summary
    } // for

    {
        import scalation.analytics.PolyRegression
        val prg   = new PolyRegression (t, z, mMax+1)
        prg.train ()
        val yp = prg.predictExpand (t)
        new Plot (t, y, yp, s"TRUTH vs. POLYNOMIAL", lines = true)
    }

} // Smoothing_FTest object



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest2` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest2
 */
object Smoothing_FTest2 extends App
{
    import scalation.analytics.clusterer.{GapStatistic, KMeansPPClusterer}
    import java.time.LocalDateTime

    val time = LocalDateTime.now()
                            .toString()
                            .replace(":", "-")
                            .replace(".", "-")

    val DATA_FILE = s"$time-data.csv"
    val SMOO_FILE = s"$time-smoothed.csv"
    val SSES_FILE = s"$time-sse.csv"

    val ord  = 4                            // b-spline order
    val kMax = 30                           // maximum number of clusters

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* LOAD DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    println ("Loading observed data...")
    val data = MatrixD ("../data/gene_expression.csv")  // observed data

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PREPARE DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    println ("Preparing data...")
    val x = data.slice (0, 1000)
    val z = new MatrixD (x.dim1, x.dim2)    // matrix for smoothed sample
    val t = VectorD.range (0, x.dim2)       // vector for time points
    val c = new MatrixD (x.dim1, x.dim2)    // matrix for smoothing spline coefficients

    x.write (DATA_FILE)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM SMOOTHING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    println ("Smoothing...")
    for (i <- x.range1) {
        val y   = x(i)
        val moo = new Smoothing_F (y, t, null, ord, lambda = -1) // smoother
        c(i)    = moo.train ()                                   // train -> set coefficients
        z(i)    = moo.predict (t)                                // predict for all time points
        //new Plot (t, y, z, s"Smoothing_F ord = $ord; λ = ${moo.getLambda}", lines = true)
    } // for

    z.write (SMOO_FILE)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM KMEANS++ CLUSTERING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    val sseObs = new VectorD (kMax)
    val sseSmo = new VectorD (kMax)
    val rsObs  = new VectorD (kMax)
    val rsSmo  = new VectorD (kMax)
    val kVals  = VectorD.range (1, kMax + 1)

    println ("Clustering... ")
    for (k <- 1 to kMax) {
        print (s"k = $k; ")
        //val (cl1, cls1) = KMeansPPClusterer (x, k)
        val cl1  = new KMeansPPClusterer (x, k)
        val cls1 = cl1.cluster ()
        val sse1 = cl1.sse ()
        val rs1  = cl1.rSquared (x)
        sseObs (k-1) = sse1
        print (s"observed SSE = $sse1; rs = $rs1; ")
        //val (cl2, cls2) = KMeansPPClusterer (z, k)
        val cl2  = new KMeansPPClusterer (z, k)
        val cls2 = cl2.cluster ()
        val sse2 = cl2.sse ()
        val rs2  = cl2.rSquared (z)
        sseSmo (k-1) = sse2        
        print (s"smoothed SSE = $sse2; rs = $rs2")
        println ()
    } // for

    val sses = MatrixD.++^ (kVals, sseObs) :^+ sseSmo
    sses.write (SSES_FILE)

    new Plot (kVals, sseObs, sseSmo, "k-means++ SSEs: Observed vs. Smoothed", true)
    new Plot (kVals, sseObs.map(math.log _), sseSmo.map(math.log _), "k-means++ log SSEs: Observed vs. Smoothed", true)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM TIGHT CLUSTERING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    // TODO implement

} // Smoothing_FTest2


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest3` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest3
 */
object Smoothing_FTest3 extends App
{
    import scalation.plot.PlotM
    import scalation.analytics.clusterer.{GapStatistic, KMeansPPClusterer, TightClusterer}
    import scalation.linalgebra.VectorD
    import scalation.util.banner

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a new random subsample.
     */
    def createSubsample (x: MatrixD, ns: Int, s: Int = 0): (MatrixD, Array [Int]) =
    {
        import scalation.random.RandomVecSample
        val rsg      = RandomVecSample (x.dim1, ns, s) 
        val indexMap = rsg.igen ().array                    // select e.g., 5, 3, 7  // FIX - why toArray
        val subsamp  = x.selectRows (indexMap)              // generate random subsample
        //println (s"subsamp = $subsamp")
        (subsamp, indexMap) 
    } // createSubsample

    import java.time.LocalDateTime

    val time = LocalDateTime.now()
                            .toString()
                            .replace(":", "-")
                            .replace(".", "-")

    val DATA_FILE = s"$time-data-observed.csv"
    val SMOO_FILE = s"$time-data-smoothed.csv"
    val COEF_FILE = s"$time-data-smoothed-coef.csv"    

    val ord  = 4                            // b-spline order
    val kMin = 5                            // maximum number of clusters
    val kMax = 6                            // maximum number of clusters

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* LOAD DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    banner ("Loading observed data...")
    val data   = MatrixD ("../data/gene_expression.csv")  // observed data
    val _0     = new VectorD (data.dim2)
    var nzdata = MatrixD (for (i <- data.range1 if data(i).sum > 100) yield data(i), false)

    println (s" -    data.min = ${data.min()}; data.max = ${data.max()}")
    println (s" -   data.dim1 = ${data.dim1}")
    println (s" - nzdata.dim1 = ${nzdata.dim1}")

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PREPARE DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    banner ("Preparing data...")
    val x = nzdata //.slice (0, 1000)

    //val s = (System.currentTimeMillis % 1000).toInt
    //val (sample, map) = createSubsample (x, 1000, s)
    //val t = VectorD.range (0, sample.dim2)
    val t      = VectorD.range (0, x.dim2)
    val sample = x

    // sample.write (DATA_FILE)    

    new PlotM (t, sample, null, s"OBSERVED", lines = true)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM SMOOTHING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    banner ("Smoothing...")
    val gap    = 4
    val smooth = new MatrixD (sample.dim1, sample.dim2)
    val coeffs = new MatrixD (sample.dim1, sample.dim2/gap + ord - 2)    // matrix for smoothing spline coefficients
    val funcs  = Array.ofDim [FunctionS2S] (sample.dim1)
    for (i <- sample.range1) {
        val y     = sample(i)
        val moo   = new Smoothing_F (y, t, null, ord, SmoothingMethod.RIDGE, -1, gap) // smoother
        coeffs(i) = moo.train()
        smooth(i) = moo.predict (t)                     // predict for all time points
        funcs(i)  = moo.predict
        //new Plot (t, y, z, s"Smoothing_F ord = $ord; λ = ${moo.getLambda}", lines = true)
    } // for

    new PlotM (t, smooth, null, s"SMOOTHED", lines = true)
    new FPlot (0.0 to 11 by 0.1, funcs.toSeq, lines = true)
    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM CLUSTERING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    def clusterKMPP (data: MatrixD, k: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): KMeansPPClusterer =
    {
        banner (s"K-Means++ Clustering: $label")
        println (s"k = $k")
        val cl    = new KMeansPPClusterer (data, k)
        val cls   = cl.cluster ()
        val sse   = cl.sse ()
        val rs    = cl.rSquared (data)
        val clust = Array.ofDim [MatrixD] (k)        
        println (s"SSE = $sse; rSquared = $rs")
        for (c <- 0 until k) {
            if (orig == null) clust(c) = MatrixD (for (i <- 0 until cls.size if cls(i) == c) yield data(i), false)
            else              clust(c) = MatrixD (for (i <- 0 until cls.size if cls(i) == c) yield orig(i), false)
            new PlotM (t, clust(c), _title = s"KM++: $label; c = $c; n = ${clust(c).dim1}", lines = true)
            if (forig != null) {
                val fclust = for (i <- 0 until cls.size if cls(i) == c) yield forig(i)
                new FPlot (0.0 to 11 by 0.1, fclust.toSeq, s"KM++: $label; c = $c; n = ${clust(c).dim1}", true)
            } // if
        } // for
        cl
    } // clusterKMPP

    def clusterTight (data: MatrixD, kMin: Int, kMax: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): TightClusterer =
    {
        banner (s"Tight Clustering: $label")
        println (s"kMin = $kMin; kMax = $kMax")
        val cl    = new TightClusterer (data, kMax, kMin)
        val cls   = cl.cluster ()
        val sse   = 0.0 // cl.sse ()
        val rs    = 0.0 // cl.rSquared (data)
        println (s"SSE = $sse; rSquared = $rs")
        for (set <- cls) {
            val clust = MatrixD (set.map(i => if (orig == null) data(i) else orig(i)).toSeq, false)
            new PlotM (t, clust, _title = s"TC: $label; n = ${clust.dim1}", lines = true)
            if (forig != null) {
                val fclust = for (i <- set) yield forig(i)
                new FPlot (0.0 to 11 by 0.1, fclust.toSeq, s"TC: $label; n = ${clust.dim1}", true)
            } // if
        } // for
        cl
    } // clusterTight

//    val ocl = clusterKMPP (sample, 6, "OBSERVED")
//    val scl = clusterKMPP (smooth, 6, "SMOOTHED")

//    println(s"new SSE = ${scl.sse (sample)}; rSquared = ${scl.rSquared (sample)}")

    val ccl = clusterKMPP (coeffs, kMax, "  COEFFS", smooth, funcs)

//    println(s"new SSE = ${ccl.sse (sample)}; rSquared = ${ccl.rSquared (sample)}")

//    clusterTight (sample, 1, 6, "OBSERVED")
//    clusterTight (smooth, 1, 6, "SMOOTHED")
    clusterTight (coeffs, kMin, kMax, "  COEFFS", smooth, funcs)        

} // Smoothing_FTest3


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest4` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest
 */
object Smoothing_FTest4 extends App
{
    import scalation.random.Normal
    import math.pow

    val s   = (System.currentTimeMillis % 1000).toInt
    val e   = Normal (0, 1.0, stream = s)
    val rho = 0.3
    val y   = VectorD.range (0, 50)
    y(0)    = e.gen; for (i <- 1 until y.dim) y(i) = rho * y(i-1) + e.gen
    val t   = VectorD.range (0, 50)                          // time points
    val τ   = null                                           // let `Smoothing_F` nake the knots
    val moo = new Smoothing_F (y, t, τ, 4, lambda = -1)      // smoother (use GCV)
    val c   = moo.train ()                                   // train -> set coefficients
    val z   = moo.predict (t)                                // predict for all time points
    moo.summary

    println (s"y = $y")
    println (s"var = ${y.variance}")

    new Plot (t, y, z, lines = true)

} // Smoothing_FTest4 object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest5` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest5
 */
object Smoothing_FTest5 extends App
{
    import scalation.plot.PlotM
    import scalation.analytics.clusterer.{GapStatistic, KMeansPPClusterer, TightClusterer}
    import scalation.linalgebra.VectorD
    import scalation.util.banner

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* LOAD DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    banner ("Loading observed data...")
    val data   = MatrixD ("/Users/mepcotterell/Desktop/polyester/polyester/chr22_small_timecourse2.csv")  // observed data
    val _0     = new VectorD (data.dim2)
    var nzdata = MatrixD (for (i <- data.range1 if data(i).sum > 100) yield data(i), false)

    println (s" -    data.min = ${data.min()}; data.max = ${data.max()}")
    println (s" -   data.dim1 = ${data.dim1}")
    println (s" - nzdata.dim1 = ${nzdata.dim1}")

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PREPARE DATASET */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    val t      = VectorD.range (0, nzdata.dim2)
    new PlotM (t, nzdata, null, s"OBSERVED", lines = true)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM SMOOTHING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    banner ("Smoothing...")
    val ord    = 4
    val gap    = 2
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

    //new PlotM (t, smooth, null, s"SMOOTHED", lines = true)
    new FPlot (0.0 to 11 by 0.1, funcs.toSeq, lines = true)

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /* PERFORM CLUSTERING */
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    val kMin = 1
    val kMax = 6

    def clusterKMPP (data: MatrixD, k: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): KMeansPPClusterer =
    {
        banner (s"K-Means++ Clustering: $label")
        println (s"k = $k")
        val cl    = new KMeansPPClusterer (data, k)
        val cls   = cl.cluster ()
        val sse   = cl.sse ()
        val rs    = cl.rSquared (data)
        val clust = Array.ofDim [MatrixD] (k)        
        println (s"SSE = $sse; rSquared = $rs")
        for (c <- 0 until k) {
            if (orig == null) clust(c) = MatrixD (for (i <- 0 until cls.size if cls(i) == c) yield data(i), false)
            else              clust(c) = MatrixD (for (i <- 0 until cls.size if cls(i) == c) yield orig(i), false)
            //new PlotM (t, clust(c), _title = s"KM++: $label; c = $c; n = ${clust(c).dim1}", lines = true)
            if (forig != null) {
                val fclust = for (i <- 0 until cls.size if cls(i) == c) yield forig(i)
                new FPlot (0.0 to 11 by 0.1, fclust.toSeq, s"KM++: $label; c = $c; n = ${clust(c).dim1}", true)
            } // if
        } // for
        cl
    } // clusterKMPP

    def clusterTight (data: MatrixD, kMin: Int, kMax: Int, label: String, orig: MatrixD = null, forig: Seq[FunctionS2S] = null): TightClusterer =
    {
        banner (s"Tight Clustering: $label")
        println (s"kMin = $kMin; kMax = $kMax")
        val cl    = new TightClusterer (data, kMax, kMin)
        val cls   = cl.cluster ()
        val sse   = 0.0 // cl.sse ()
        val rs    = 0.0 // cl.rSquared (data)
        println (s"SSE = $sse; rSquared = $rs")
        for (set <- cls) {
            val clust = MatrixD (set.map(i => if (orig == null) data(i) else orig(i)).toSeq, false)
            //new PlotM (t, clust, _title = s"TC: $label; n = ${clust.dim1}", lines = true)
            if (forig != null) {
                val fclust = for (i <- set) yield forig(i)
                new FPlot (0.0 to 11 by 0.1, fclust.toSeq, s"TC: $label; n = ${clust.dim1}", true)
            } // if
        } // for
        cl
    } // clusterTight

//    val ocl = clusterKMPP (sample, 6, "OBSERVED")
//    val scl = clusterKMPP (smooth, 6, "SMOOTHED")

//    println(s"new SSE = ${scl.sse (sample)}; rSquared = ${scl.rSquared (sample)}")

      val ccl = clusterKMPP (coeffs, kMax, "  COEFFS", smooth, funcs)

//    println(s"new SSE = ${ccl.sse (sample)}; rSquared = ${ccl.rSquared (sample)}")

//    clusterTight (sample, 1, 6, "OBSERVED")
//    clusterTight (smooth, 1, 6, "SMOOTHED")
    clusterTight (coeffs, kMin, kMax, "  COEFFS", smooth, funcs)        

} // Smoothing_FTest5

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Smoothing_FTest6` is used to test the `Smoothing_F` class.
 *  > run-main scalation.analytics.fda.Smoothing_FTest6
 */
object Smoothing_FTest6 extends App
{
    import scalation.plot.PlotM
    import scalation.analytics.clusterer.{GapStatistic, KMeansPPClusterer, TightClusterer}
    import scalation.linalgebra.VectorD
    import scalation.util.banner
    import scala.collection.mutable.{ArrayBuffer,Set}
    import scalation.relalgebra.Relation
    
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
    new PlotM (t, nzdata, null, s"OBSERVED", lines = true)
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
    new FPlot (t0 to tn by 0.1, funcs.toSeq, s"SmoothedData", lines = true)

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
