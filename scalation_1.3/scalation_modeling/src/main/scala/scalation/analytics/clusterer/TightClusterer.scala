
package scalation.analytics.clusterer

import scalation.linalgebra.MatrixD
import scalation.random.RandomVecSample

class TightClusterer (x: MatrixD, k: Int, s: Int = 0)
{
    val b  = 10                                           // times to resample
    val ns = (x.dim1 * 0.7).toInt                         // size of random sample
    val rs = RandomVecSample (x.dim1, ns, s)              // random sample generator
    val d  = new MatrixD (ns, ns)                         // comembership matrix
    val md = new MatrixD (ns, ns)                         // mean comembership matrix

    for (l <- 0 until b) {
        val y = x.selectRows (rs.igen ().toArray)         // generate random subsample  // FIX - why toArray
        println (s"y = $y")                               
        val kmc    = new KMeansClustering (y, k, s)       // apply clustering to the sample
        val clustr = kmc.cluster ()                       // get the clusters
        println (s"clustr = ${clustr.deep}")

        for (i <- y.range1; j <- y.range1) d(i, j) = if (clustr(i) == clustr(j)) 1.0 else 0.0
        println (s"d = $d")
        md += d
    } // for

    md /= b                                               // compute mean
    println (s"mean = $md")

} // TightClusterer class


// > run-main scalation.analytics.clusterer.TightClustererTest

object TightClustererTest extends App
{
    val v = new MatrixD ((6, 2), 1.0, 2.0,
                                 2.0, 1.0,
                                 5.0, 4.0,
                                 4.0, 5.0,
                                 9.0, 8.0,
                                 8.0, 9.0)

   new TightClusterer (v, 3)

} // TightClustererTest object

