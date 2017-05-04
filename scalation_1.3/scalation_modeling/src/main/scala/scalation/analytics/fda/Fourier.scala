
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller, Michael Cotterell
 *  @version 1.3
 *  @date    Tue May 2 21:45:58 EDT 2017
 *  @see     LICENSE (MIT style license file).
 *
 *  @see https://en.wikipedia.org/wiki/Fourier_series
 */

package scalation.analytics.fda

import math.{cos, Pi, sin}

import scalation.linalgebra.{MatrixD, VectorD}
import scalation.math.double_exp
import scalation.plot.FPlot
import scalation.util.Error

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Fourier` class provides Fourier basis functions. Such basis functions
 *  are useful are useful for fitting periodic data in Functional Data Analysis.
 *  @see https://en.wikipedia.org/wiki/Fourier_series
 *-----------------------------------------------------------------------------
 *  @param k     the number of terms in the series (k >= 1)
 *  @param w     the fundamental frequency parameter
 */
class Fourier (k: Int = 4, w: Double = 2.0 * Pi)
      extends Error
{
    private val DEBUG = true                                // debug flag
    private val ns    = 2 * k + 1

    def apply (j: Int, t: Double): Double = if (j == 0)          1.0
                                            else if (j % 2 == 0) sin ((j+1)/2 * w * t)
                                            else                 cos ((j+1)/2 * w * t)

    override def toString =
    {
        var s = s"Fourier(k = $k, w = $w) \n"
        for (j <- range) if (j == 0)          s +=  "   1.0 \n"
                         else if (j % 2 == 0) s += s" + sin (${(j+1)/2} * w * t) \n"
                         else                 s += s" + cos (${(j+1)/2} * w * t)"
        s
    } // toString
    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    /** Computes a matrix of the basis functions evaluated at the given time
     *  points.
     *  @param k  the order of the basis expansion (# terms = 2k+1)
     *  @param t  the time parameter
     */
    def phi (t: VectorD): MatrixD =
    {
        val nt = t.dim
        val Φ  = new MatrixD (nt, ns)
        for (i <- Φ.range1; j <- Φ.range2) Φ(i, j) = this (j, t(i))
        Φ
    } // phi

    def range = 0 until 2 * k + 1

} // Fourier

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `FourierTest` object is used to test the `Fourier` class.
 *  > run-main scalation.analytics.fda.FourierTest
 */
object FourierTest extends App
{
    import scala.language.postfixOps

    // census data (1790--1990)
    val t = VectorD (1790.0 to 1991 by 10 toSeq)
    println (t)
    val y = VectorD (  3.9000,
                       5.3000,
                       7.2000,
                       9.6000,
                      12.9000,
                      17.1000,
                      23.1000,
                      31.4000,
                      38.6000,
                      50.2000,
                      62.9000,
                      76.0000,
                      92.0000,
                     105.7000,
                     122.8000,
                     131.7000,
                     150.7000,
                     179.0000,
                     205.0000,
                     226.5000,
                     248.7000)

    import scalation.plot.Plot
    new Plot (t, y, lines = true)

    val k    = 8 
    val L    = y.max () - y.min ()
    val w    = 2.0 * Pi / L
    val four = new Fourier (k, w) // 0.020282546792471)
    val Φ    = four.phi (t)
    val I    = MatrixD.eye (t.dim)
    val λ    = 1.0E-10
    val c    = ((Φ.t * Φ) + (I * λ)).inverse * Φ.t * y

    def x (tt: Double): Double =
    {
        var sum = 0.0
        for (j <- four.range) sum += c(j) * four (j, tt)
        sum
    } // predict

    val z    = VectorD (for (tt <- 0 until t.dim) yield x (tt))
    val e    = y - z
    val sse  = e dot e

    println (s"Φ = $Φ")
    println (s"c = $c")
    println (four)

    import scalation.stat.vectorD2StatVector

    for (i <- 1 until t.dim) println (y.acorr(i) / (2.0 * Pi))

    new FPlot (1789.0 to 1991 by 1, Seq(x), lines = true)

} // FourierTest
