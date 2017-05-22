//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.3
 *  @date    Sat Apr 29 13:22:47 EDT 2017
 *  @see     LICENSE (MIT style license file).
 *
 *  @see introcs.cs.princeton.edu/java/92symbolic/Polynomial.java.html
 */


package scalation.calculus

import scala.math.abs

import scalation.linalgebra.VectorD

import scala.collection.immutable.Vector


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `PiecewisePoly` (pwp) class provides simple operations on piecewise polynomials.
 *  There are as many pairs of booleans in the knot vector as there are intervals defining the pwp:
 *  	The first member in each pair indicates whether the left endpoint of the corresponding interval
 *  	is an inclusive (True) end point or an exclusive (False) end point,
 *	and similarly for the second member and right end point. 
 *  @param knot   the 'knot' vector describing the intervals on which the pwp is defined on
 *  @param ends   the list of ordered pairs describing the inclusion/exclusion rules for endpoints on each interval of the pwp
 *  @param polys  the polynomials defined on each interval of the knot vector
 */

// TODO : knot should probably be a vector of Interval type
case class PiecewisePoly (knot: VectorD, ends: Vector[(Boolean,Boolean)], polys: Vector[Poly])
{
	
    if ( knot.size-1 != ends.length  ) println(s"Error in constructor: number of end points must match number of knot intervals.")
    if ( knot.size-1 != polys.length ) println(s"Error in constructor: number of knot intervals must match number of polynmoials.")
    if ( ends.length != polys.length ) println(s"Error in constructor: number of end points must match number of polynomials.")
    if ( endPointsInvalid ) 	       println(s"Error in constructor: bad end points description.")

    assert( knot.size-1 == ends.length && knot.size-1 == polys.length && ends.length == polys.length && !endPointsInvalid )
    
    private val DEBUG = true                         // debug flag
    val intervals:Int = polys.length                // number of intervals the piecewise polynomial is defined on

    if (DEBUG) for ( i <- 1 to intervals ) println (s"PiecwisePoly defined on ${knot(i-1)} to ${knot(i)} as ${polys(i-1)}")

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Apply/evaluate the polynomial at 'x'.
     *  @param x  the value of the variable
     */
    def apply (x: Double): Double =
    {
	val interval = getApplyInterval(x)-1
	if ( interval < 0 ) println(s"Error at apply: tried to evaluate for an undefined value.")
	assert(interval >= 0)
	if ( DEBUG ) println(s"Piecewise function for evaluating piecewise function at $x: ${polys(interval)}")
        return polys(interval)(x) 
    } // apply

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Add 'this' piecewise polynomial and the 'q' piecewise polynomial.
     *  @param q  the other polynomial
     */

    //TODO maybe the new piecewise poly should be defined on q.knot INTERSECT this.knot
    def + (q: PiecewisePoly): PiecewisePoly =
    {
	assert( knot.equals(q.knot) && ends.equals(q.ends) )
       	val newFuncs = for(i <- 0 until intervals) yield polys(i)+q.polys(i)
	PiecewisePoly(knot,ends,newFuncs.toVector)
    } // +

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Subtract 'this' piecewise polynomial and the 'q' piecewise polynomial.
     *  @param q  the other polynomial
     */
    def - (q: PiecewisePoly): PiecewisePoly =
    {
	assert( knot.equals(q.knot) && ends.equals(q.ends) )
       	val newFuncs = for(i <- 0 until intervals) yield polys(i)-q.polys(i)
	PiecewisePoly(knot,ends,newFuncs.toVector)
    } // -

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Multiply 'this' piecewise polynomial and the 'q' piecewise polynomial.
     *  @param q  the other polynomial
     */
    def * (q: PiecewisePoly): PiecewisePoly =
    {
	assert( knot.equals(q.knot) && ends.equals(q.ends) )
       	val newFuncs = for(i <- 0 until intervals) yield polys(i)*q.polys(i)
	PiecewisePoly(knot,ends,newFuncs.toVector)
    } // *

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Multiply 'this' piecewise polynomial and the polynomial 'q' NOT a piecewise polynomial.
     *  @param q  the other polynomial
     */
    def * (q: Poly): PiecewisePoly =
    {
       	val newFuncs = for(i <- 0 until intervals) yield q*polys(i)
	PiecewisePoly(knot,ends,newFuncs.toVector)
    } // **

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Take the derivative of 'this' piecewise polynomial,
     *  returning the result as a piecewise polynomial.
     */
    def derivative: PiecewisePoly =
    {
	val polys = for( i<- 0 until intervals ) yield this.polys(i).derivative
	PiecewisePoly (knot,ends,polys.toVector)
    }

    def Ⅾ : PiecewisePoly =
    {
        val polys = for( i<- 0 to intervals ) yield this.polys(i).derivative
	PiecewisePoly (knot,ends,polys.toVector)
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Integrate 'this' polynomial, returning the result as a polynomial.
     *  Note, the arbitrary constant 'c' for the indefinite integral is set to 1.
     */
    def integrate: PiecewisePoly =
    {
	val polys = for( i<- 0 until intervals ) yield this.polys(i).integrate
	PiecewisePoly (knot,ends,polys.toVector)
    }

    def ∫ : PiecewisePoly =
    {
	val polys = for( i<- 0 to intervals ) yield this.polys(i).integrate
	PiecewisePoly (knot,ends,polys.toVector)
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Integrate 'this' polynomial on the interval 'on', returning its value as
     *  a double.
     *  @param on  the interval of integration
     */
    def integrate (on: Interval): Double =
    {
	var (start,finish) = (0,intervals+1)
	val (fKnot,lKnot)  = (knot(0),knot(intervals))
	val (begin,end)    = (on._1,on._2)
	if    ( begin <  fKnot || end >  lKnot || begin >  end ) println(s"Invalid integration interval: ${on}")
	assert( begin >= fKnot && end <= lKnot && begin <= end )
	for ( i <- 1 to intervals ){
	    if ( start  == 0       && knot(i) >= begin ) start  = i
	    if ( finish == intervals+1  && knot(i) >= end   ) finish = i
	} // for
	var sum = 0.0
	for ( i <- start to finish ){
	    val a = if (knot(i-1) < begin) begin else knot(i-1)
	    val b = if (knot(i) > end) end else knot(i)
	    if (DEBUG) println(s"Integrating ${polys(i-1)} on ($a,$b)")
	    sum += polys(i-1).integrate((a,b))
	}
	sum
    } // integrate

    def ∫ (on: Interval): Double = integrate(on)
    
    def getApplyInterval(x: Double) : Int =
    {
	var ret = 0
	if ( x < knot(0) || x > knot(knot.size-1) ) {
	   println("invalid apply interval")
	} // if 
    	for( i <- 1 to intervals) {
	     if ( x == knot(i-1) && ends(i)._1 ) ret = i
	     else if ( x > knot(i-1) && x < knot(i) ) {
	     	ret = i
	     } // if
	     else if ( x == knot(i) && ends(i)._2 ) ret = i
	} // for
	ret
    }

    /*TODO more sophisticated endpoints checking.
    	   only endpoint checking now makes sure that the function is defined over the entire not vector.
	   maybe we should be more flexible and just make sure that the general functional requirements are ensured?
	   	 i.e. - f(x1)=f(x2) <-> x1=x2
    */
    
    def endPointsInvalid : Boolean =
    {
	if  ( ! (ends(0)._1 && ends(ends.length-1)._2 ) )return true 
	for ( i <- 1 to intervals-1 ) if( !(ends(i-1)._2 ^ ends(i)._1 && ends(i)._2 ^ ends(i+1)._1) ) return true
	false
    } // endPointsInvalid

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Convert the polynomial to a string.
     */
    override def toString: String =
    {
	var ret = "\n"
        for( i <- 0 to intervals-1 ){
	     val l = if (ends(i)._1) "[" else "("
	     val r = if (ends(i)._2) "]" else ")"
	     ret += s"${l}${knot(i)},${knot(i+1)}${r} : ${polys(i)}\n"
	} 
	ret += "\n"
	ret
    } // toString

} // PiecewisePoly class



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `PolyTest` object is used to test the `Poly` class.
 *  > run-main scalation.calculus.PiecewisePolyTest
 */
object PiecewisePolyTest extends App
{
    import scalation.plot.Plot
    import scalation.calculus
    val mpoly = Poly(1,1,1)
    val funs  = Array.ofDim[Poly](4)
    val funs2 = Array.ofDim[Poly](4)
    val knots = VectorD(0.0,1.0,2.0,3.0,4.0)
    val p = (true,false)
    val ends = Vector(p,p,p,(true,true))
    println("//::::::::::::::::::::::::::::\nCreating Polynomials\n//::::::::::::::::::::::::::::")    
    for ( i <- 0 until 4 ) {
    	funs(i)  = Poly(i+1,i+1,i+1,i+1)
    	funs2(i) = Poly(i+2,i+2,i+2,i+2)
    } // for
    println("\n\n")
    println("//::::::::::::::::::::::::::::\nCreating Piecewise Polynomials: pwp1 & pwp2\n//::::::::::::::::::::::::::::")
    println("\n\n")
    val pwp1 = PiecewisePoly(knots,ends,funs.toVector)
        println("\n\n")
    val pwp2 = PiecewisePoly(knots,ends,funs2.toVector)
        println("\n\n")
    println("//::::::::::::::::::::::::::::\nDifferentiating Piecewise Polynomials: dpl1 & dpl2\n//::::::::::::::::::::::::::::")
    println("\n\n")
    val dpl1 = pwp1.derivative                         // its derivative
        println("\n\n")
    val dpl2 = pwp2.derivative
        println("\n\n")
    println("//::::::::::::::::::::::::::::\nIntegrating (Indefinite) Piecewise Polynomials: ipl1 & ipl2\n//::::::::::::::::::::::::::::")
    println("\n\n")
    val ipl1 = pwp1.integrate                          // one of its indefinite integrals
    println("\n\n")
    val ipl2 = pwp2.integrate                          // one of its indefinite integrals
    println("\n\n")
    println("//::::::::::::::::::::::::::::\nIntegrating (Definite) Piecewise Polynomials: jpl1 & jpl2\n//::::::::::::::::::::::::::::")
        println("\n\n")
    val jpl1 = pwp1.integrate ((0.0, 2.0))             // one of its definite integrals
    println("\n\n")
    val jpl2 = pwp2.integrate ((0.0, 2.0))             // one of its definite integrals
    println("\n\n")
    println("//::::::::::::::::::::::::::::\nSum Of Piecewise Polynomials and Their Derivatives: spl1 & spl2\n//::::::::::::::::::::::::::::")
    println("\n\n")
    val spl1 = pwp1 + dpl1                              // sum of polynomials and its dervivate
        println("\n\n")
    val spl2 = pwp2 + dpl2				// sum of polynomials and its dervivate
        println("\n\n")
    println("//::::::::::::::::::::::::::::\nDifference of Piecewise Polynomials and Their Derivatives: mpl1 & mpl2\n//::::::::::::::::::::::::::::")
    println("\n\n")
    val mpl1 = pwp1 - dpl1                              // difference of polynomial and its dervivate
    println("\n\n")
    val mpl2 = pwp2 - dpl2                              // difference of polynomial and its dervivate
    println("\n\n")
    println("//::::::::::::::::::::::::::::\nProduct of Piecewise Polynomials and Their Derivatives: tpl1 & tpl2\n//::::::::::::::::::::::::::::")
    println("\n\n")
    val tpl1 = pwp1 * dpl1                              // product of polynomial and its dervivate
    println("\n\n")
    val tpl2 = pwp2 * dpl2                              // product of polynomial and its dervivate
    println("//::::::::::::::::::::::::::::\nProduct of Piecewise Polynomials and Non-Piecewise Polynomials: zpl1 & zpl2\n//::::::::::::::::::::::::::::")
    println("\n\n")
    val zpl1 = pwp1 * mpoly                              // product of polynomial and its dervivate
    println("\n\n")
    val zpl2 = pwp2 * mpoly                             // product of polynomial and its dervivate
    println("\n\n")
    println (s"pwp1      = $pwp1")
    println (s"dpl1     = $dpl1")
    println (s"ipl1     = $ipl1")
    println (s"spl1     = $spl1")
    println (s"mpl1     = $mpl1")
    println (s"tpl1     = $tpl1")
    println (s"zpl1     = $zpl1")
  
    println (s"pwp1 (2)  = ${pwp1 (2)}")
    println (s"dpl1 (2) = ${dpl1 (2)}")
    println (s"ipl1 (2) = ${ipl1 (2)}")
    println (s"jpl1     = $jpl1")
    println (s"spl1 (2) = ${spl1 (2)}")
    println (s"mpl1 (2) = ${mpl1 (2)}")
    println (s"tpl1 (2) = ${tpl1 (2)}")
    println (s"zpl1 (2) = ${zpl1 (2)}")
    
    println (s"pwp2      = $pwp2")
    println (s"dpl2     = $dpl2")
    println (s"ipl2     = $ipl2")
    println (s"spl2     = $spl2")
    println (s"mpl2     = $mpl2")
    println (s"tpl2     = $tpl2")
    println (s"zpl2     = $zpl2")
    println (s"pwp2 (2)  = ${pwp2 (2)}")
    println (s"dpl2 (2) = ${dpl2 (2)}")
    println (s"ipl2 (2) = ${ipl2 (2)}")
    println (s"jpl2     = $jpl2")
    println (s"spl2 (2) = ${spl2 (2)}")
    println (s"mpl2 (2) = ${mpl2 (2)}")
    println (s"tpl2 (2) = ${tpl2 (2)}")
    println (s"zpl2 (2) = ${zpl2 (2)}")

/*    val x = VectorD.range (0, 20) / 5.0
    val y = x.map (pl (_))
    val z = x.map (dpl (_))

    new Plot (x, y, z)
*/

} // PolyTest object

