
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell
 *  @version 1.2
 *  @date    Thu Jan 25 16:05:58 EST 2017
 *  @see     LICENSE (MIT style license file).
 */

package apps.analytics.fda

import scala.io.Source.fromFile
import scalation.analytics.{BASE_DIR, KMeansClustering}
import scalation.analytics.fda.Smoothing_F
import scalation.linalgebra.{MatrixD, VectorD}
import scalation.util.banner

object GeneIsoform extends App
{

    // LOAD DATASET
    banner ("Loading Dataset")

    val file   = BASE_DIR + "fda/gene_isoform.csv"
    val sp     = ','
    val lskip  = 1 // skip first line (headers)
    val cskip  = 1 // skip first column (ids)
    val lines  = fromFile (file).getLines.toArray
    val (m, n) = (lines.length-lskip, lines(0).split(sp).length-cskip)
    var data   = new MatrixD (m, n)
    //val fits   = new MatrixD (m, n+1) // +1 for intercept
    val fits   = new MatrixD (m, n-1)

    for (i <- data.range1) {
        val line = lines(i+lskip).split (sp)
        for (j <- data.range2) {
            data(i, j) = line(j+cskip).toDouble
        } // for
    } // for

    // intercept model (may not need it)
    // data = VectorD.one (m) +^: data

    println(data)

    // FIT EACH ISOFORM TO SMOOTHING SPLINE
    banner ("Fitting Models with Optimal λ-values")

    for (i <- data.range1) {
        val ord = 4
        val y   = data(i) 
        val t   = VectorD.range (1, y.dim+1)
        // val τ   = VectorD.range (1, 3) / 3.0
        val moo = new Smoothing_F (y, t, null, ord)
        val c   = moo.train()
        println (s"c = $c")
        fits(i) = c(1 until c.dim)
    } // for

    println(fits)
   
    // CLUSTER USIGN K-MEANS
    banner ("Clustering with K-Means (k=5)")
    val k  = 5 // # clusters (will use gap estimate later)
    val cl = new KMeansClustering (fits, k)
    println ("--- final cluster = " + cl.cluster ().deep)

} // GeneIsoform

