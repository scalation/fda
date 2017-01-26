
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author   Vishnu Gowda Harish, John Miller
 *  @version  1.0
 *  @date     Mon March 28 5:10:20 EDT 2016
 *  @see      LICENSE (MIT style license file).
 */

package scalation.linalgebra

import scala.collection.{breakOut, Traversable}
import scala.collection.mutable.{HashSet, IndexedSeq}
import scala.math.{abs => ABS, max => MAX, min => MIN, sqrt}

import scalation.math.long_exp

import scalation.util.{ReArray, SortingL}
import scalation.util.SortingL.{iqsort, qsort2}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `TripletL` class holds a run length encoded triplet, e.g., (1l, 2, 3)
 *  Triplets are used for storing compressed vectors and matrices.
 *  @param value     the value itself
 *  @param count     the number of times the value is repeated in a run
 *  @param startPos  the starting position of the value
 */
class TripletL (var value: Long, var count: Int, var startPos: Int)
{
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Update 'this' triplet.
     *  @param x  the tuple value to which 'this' triplet should be updated
     */
    def update (x: (Long, Int, Int))
    {
        value = x._1; count = x._2; startPos = x._3
    } // update

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Override 'equals' to determine whether 'this' triplet equals triplet 'y'.
     *  @param y  the triplet to compare with this
     */
    override def equals (y: Any): Boolean =
    {
        val yy = y.asInstanceOf [TripletL]
        value =~ yy.value && startPos == yy.startPos && count == yy.count
    } // equals

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Must also override 'hashCode' for 'this' vector to be compatible with equals.
     */
    override def hashCode (): Int = ABS ((value.hashCode()) + (31 * startPos) + (13 * count))

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Converts 'this' triplet to a string.
     */
    override def toString: String = "(" + value + ", " + count + ", " + startPos + ")"

} // TripletL class


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `RleVectorL` object is the companion object for the `RleVectorL` class.
 */
object RleVectorL 
{
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a `RleVectorL` (i.e., Run Length Encoding compression) from `VectorL`.
     *  @param y  the `VectorL` to be compressed
     */
    def apply (y: VectoL): RleVectorL = 
    {
        if (y.size == 0) return new RleVectorL (0)
        val x = y.asInstanceOf [VectorL]
        var c1 = new ReArray [TripletL] ()
        var idx = 0
        var count = 1
        var z = 0  
        for (i <- 1 until x.dim) {
            if (x(i) != x(i - 1)) {
                c1(idx) = new TripletL (x(i - 1), count,(i) - count); idx += 1
                count = 1; z += 1
            } else {
                count += 1
            } // if
        } // for
        c1(idx) = new TripletL (x(x.dim - 1), count, x.dim - count)
        new RleVectorL (x.dim, c1)
    } //apply 
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a `RleVectorL` from a sequence of `Long`.
     *  @param xs  the sequence of the `Long`
     */
    def apply (xs: Seq [Long]): RleVectorL =
    { 
        if (xs.length == 0 )  return new RleVectorL (0)
        var c1 = new ReArray [TripletL] ()
        var idx = 0
        var count = 1
        var z = 0  
        for (i <- 1 until xs.size) {
            if (xs(i) != xs(i - 1)) {
                c1(idx) = new TripletL (xs(i - 1), count, (i) - count); idx += 1
                count = 1; z += 1
            } else { count += 1 } // if
        } // for
       c1(idx) = new TripletL (xs(xs.size - 1), count, xs.size - count)
       new RleVectorL (xs.size, c1)
    } // apply
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a `RleVectorL` from one or more values (repeated values `Long`*).
     *  @param x   the first `Long number
     *  @param xs  the rest of the `Long` numbers
     */
    def apply (x: Long, xs: Long*): RleVectorL =
    {
        var c1 = new ReArray [TripletL] ()
        var idx = 0
        var count = 1
        var z = 0
        if (xs(1) != x) {
            c1(idx) = new TripletL (x, 1, 0); idx += 1
            count = 1; z += 1
        } else { count += 1 } // if
        for (i <- 2 until xs.size) {      
            if (xs(i) != xs(i - 1)) {
                c1(idx) = new TripletL (xs(i - 1), count, i - count); idx += 1
                count = 1; z += 1
            } else { count += 1 } // if
        } // for
        c1(idx) = new TripletL (xs(xs.size - 1), count, xs.size - count)
        new RleVectorL (xs.size + 1, c1)   
    } // apply

} // RleVectorL object


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `RleVectorL` class stores and operates on compressed Numeric Vectors 
 *  of base type `Long`. 
 *  @param dim  the dimension/size of the vector (uncompressed data)
 *  @param v    the 1D array used to store vector elements
 */   
class RleVectorL (val dim: Int, protected var v: ReArray [TripletL] = null)  
      extends VectoL
{      
    if (v == null) {
        v = new ReArray [TripletL] () 
        v(0) = new TripletL (0l, dim, 0)
    } // if
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Number of elements in the compressed vector `RleVectorL` as a double.
     */
    def cnd: Double = v.length.toDouble
           
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Range for the compressed storage array.
     */
    def crange: Range = 0 until v.length
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the size (number of elements) of compressed vector `RleVectorL`.
     */
    def csize: Int = v.length
     
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Construct a vector and assign values from vector 'u'.
     *  @param u  the other vector
     */
    def this (u: VectoL) { this (u.dim); for (i <- range) this(i) = u(i) }
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Clone (make a copy) a vector from this vector.
     */
    def copy: RleVectorL = new RleVectorL (this)
  
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Get the uncompressed vector's element at the 'i'-th index position. 
     *  @param i  the given index
     */
    def apply (i: Int): Long = v(getTriplet(i)).value

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Get 'this' vector's elements within the given range (vector slicing).
     *  @param r  the given range
     */
    def apply (r: Range): VectorL = slice (r.start, r.end)
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Get 'this' vector's entire data as an array.
     */
    def apply (): IndexedSeq [Long] = 
    {
        (for (i <- v; j <- 0 until i.count) yield i.value) (breakOut)                                            
    } // apply
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Set 'this' vector's element at the 'i'-th index position. 
     *  @param i  the given index
     *  @param x  the value to assign
     */
    def update (i: Int, x: Long)
    {      
        var tripPos = getTriplet (i)      
        var (value, count, sPos) = (v(tripPos).value, v(tripPos).count, v(tripPos).startPos) 
           
        if (x =~ value) return 
      
        if (count == 1) mergeBoth (tripPos, x)                 
        else if (i == sPos && count > 1) mergeLeft (tripPos, x)
        else if (i == sPos + count - 1 && count > 1) mergeRight (tripPos, x)            
        else split (i, tripPos, x)      
    } // update
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Merge the 'tripPos' triplet with the right triplet. 
     *  @param tripPos the position of the triplet to be updated
     *  @param value   the new value  
     */
    private def mergeRight (tripPos: Int, value: Long) = 
    {       
        if (tripPos + 1 >= csize) { 
            v(tripPos).count -= 1
            v(tripPos + 1) = new TripletL (value, 1, v(tripPos).count + v(tripPos).startPos)      
        } else if (value =~ v(tripPos + 1).value) {
            v(tripPos).count -= 1
            v(tripPos + 1).count += 1; v(tripPos + 1).startPos -= 1
        } else {
            v(tripPos).count -= 1; 
            v.shiftRight (tripPos + 1);
            v(tripPos + 1) = new TripletL (value, 1, v(tripPos).count + v(tripPos).startPos)
        } // if           
    } // mergeRight
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Merge the 'tripPos' triplet with the left triplet.  
     *  @param tripPos the position of the triplet to be updated
     *  @param value   the new value  
     */
    private def mergeLeft (tripPos: Int,value: Long) = 
    {
        if (tripPos == 0) {
            v.shiftRight (tripPos)
            v(tripPos) = new TripletL (value, 1, 0)
            v(tripPos + 1).count -= 1; v(tripPos + 1).startPos = 1
        } else if (value =~ v(tripPos - 1).value) {
            v(tripPos).count -= 1; v(tripPos).startPos += 1           
            v(tripPos - 1).count += 1           
        } else {
            v.shiftRight (tripPos)           
            v(tripPos) = new TripletL (value, 1, v(tripPos + 1).startPos)
            v(tripPos + 1).count -= 1; v(tripPos + 1).startPos += 1
        } // if
    } // mergeLeft
 
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Split the current triplet in to three triplets. 
     *  @param idx     the position of the update
     *  @param tripPos the position of the triplet to be updated
     *  @param value   the new value  
     */
    private def split (idx: Int, tripPos: Int, value: Long) = 
    { 
        v.shiftRight2 (tripPos)
        val oldCount = v(tripPos).count             
        v(tripPos).count = idx - v(tripPos).startPos
        v(tripPos + 1) = new TripletL (value, 1, v(tripPos).startPos + v(tripPos).count)
        v(tripPos + 2) = new TripletL (v(tripPos).value, oldCount - (v(tripPos).count + 1), 
                                       v(tripPos + 1).startPos + v(tripPos + 1).count)       
    } // split
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Merge the 'tripPos' triplet and its right triplet into its left triplet 
     *  @param tripPos the position of the triplet to be updated
     *  @param value   the new value  
     */
    private def mergeBoth (tripPos: Int, value: Long) = 
    {
        if (tripPos == 0 && csize > 1) {
            if (value =~ v(tripPos + 1).value) {
                v(1).startPos = 0; v(1).count += 1
                v.shiftLeft (tripPos)              
            } else v(0).value = value      
        } else if (tripPos == csize - 1 && csize > 1 ) {
            if (value =~ v(tripPos - 1).value) {
                v(tripPos - 1).count += 1
                v.shiftLeft (tripPos)
            } else v(tripPos).value = value 
        } else if (csize > 1) {
            val curr = v(tripPos)
            val prev = v(tripPos - 1)
            val next = v(tripPos + 1)
            if (value =~ prev.value && value =~ next.value) {
                prev.count = prev.count + 1 + next.count
                v.shiftLeft2 (tripPos) 
            } else if (value =~ prev.value) {
                prev.count += 1
                v.shiftLeft (tripPos)
            } else if (value =~ next.value) {
                next.count += 1; next.startPos -= 1
                v.shiftLeft (tripPos)
            } else curr.value = value 
        } else v(0).value = value 
    } // mergeBoth 
   
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Get the triplet corresponding to the element at the 'i'-th index position. 
     *  @param i the given index
     */
    def getTriplet (i: Int): Int = 
    {   
        if (i < 0 || i >= dim) throw new IndexOutOfBoundsException (i.toString)        
        var start = 0
        var end   = csize - 1
        var mid   = 0
        while (start != end) {
            mid = (start + end) / 2
            if (v(mid).startPos < i) start = mid + 1 else end = mid	        
        } // while                           
        if (i == 0 || i == v(end).startPos || i > v(end).startPos) end else end - 1                    
    } // getTriplet
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Set 'this' vector's elements over the given range (vector slicing).
     *  @param r  the given range
     *  @param x  the value to assign
     */
    def update (r: Range, x: Long) { for (i <- r) this(i) = x }
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Set 'this' vector's elements over the given range (vector slicing).
     *  @param r  the given range
     *  @param u  the vector to assign
     */
    def update (r: Range, u: VectoL) { for (i <- r) this(i) = u(i - r.start) }
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Set each value in 'this' vector to 'x'.
     *  @param x  the value to be assigned
     */
    def set (x: Long) 
    {      
        v = new ReArray [TripletL] (1) 
        v(0) = new TripletL (x, dim, 0)
    } // set
   
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Set the values in 'this' vector to the values in array 'u'.
     *  @param u  the array of values to be assigned
     */
    def set (u: Seq [Long]) { for (i <- range)  this(i) = u(i) } 
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Iterate over 'this' vector element by element.
     *  @param f  the function to apply
     */
    def foreach [U] (f: Long => U)
    {
        var i = 0    
        while (i < csize) { f (v(i).value); i += 1 }
    } // foreach
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Convert 'this' `RleVectorL` into a `VectorI`.
     */
    def toInt: VectorI =
    {
        val ar = v()
        VectorI (for (i <- 0 until v.length; j <-0 until ar(i).count) yield ar(i).value.toInt)
    } // toInt
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Convert 'this' `RleVectorL` into a `VectorL`.
     */
    def toLong: VectorL =
    {
        val ar = v()
        VectorL (for (i <- 0 until v.length ; j <- 0 until ar(i).count) yield ar(i).value.toLong)
    } // toLong

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Convert 'this' `RleVectorL` into a `VectorD`.
     */
    def toDouble: VectorD =
    {   
        val ar = v()
        VectorD (for (i <- 0 until v.length; j <- 0 until ar(i).count) yield ar(i).value.toDouble)
    } // toDouble
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Expand the size 'dim' of 'this' vector by 'more' elements.
     *  @param more  the number of new elements to add
     */
    def expand (more: Int = dim): VectoL =
    {     
        var c = toDense
        if (more < 1) c       // no change
        else new VectorL (dim + more, Array.concat (c().toArray, new Array [Long] (more)))
    } // expand
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a zero vector (all elements are zero) of length 'size'.
     *  @param size  the size of the new vector
     */
    def zero (size: Int): RleVectorL = new RleVectorL (size)
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a one vector (all elements are one) of length 'size'.
     *  @param size  the size of the new vector
     */
    def one (size: Int): RleVectorL = 
    {
        val trip = new TripletL (1l, size, 0)
        new RleVectorL (size, ReArray.fill (1) (trip))
    } // one

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create an compressed vector of the form (0, ... 1l, ... 0) where the 1l 
     *  is at position j.
     *  @param j     the position to place the 1l
     *  @param size  the size of the vector (upper bound = size - 1)
     */
    def oneAt (j: Int, size: Int = dim): RleVectorL =
    {   
        val c = new RleVectorL (size)
        c(j)  = 1l
        c
    } // oneAt
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create an compressed vector of the form (0, ... -1l, ... 0) where the -1l 
     *  is at position j.
     *  @param j     the position to place the -1l
     *  @param size  the size of the vector (upper bound = size - 1)
     */
    def _oneAt (j: Int, size: Int = dim): RleVectorL =
    {
        val c = new RleVectorL (size)
        c(j)  = -1l
        c
    } // _oneAt
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Filter the elements of 'this' vector based on the predicate 'p', returning
     *  a new vector.
     *  @param p  the predicate (`Boolean` function) to apply
     */
    override def filter (p: Long => Boolean): RleVectorL = 
    {    
        val c = new ReArray [Long] ()
        var idx = 0
        for (i <- crange if p (v(i).value)) 
            for (j <-  0 until v(i).count) { c(idx) = v(i).value; idx += 1 }      
        RleVectorL (c)    
    } // filter
   
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Filter the elements of 'this' vector based on the predicate 'p', returning
     *  the index positions.
     *  @param p  the predicate (`Boolean` function) to apply
     */
    def filterPos (p: Long => Boolean): IndexedSeq [Int] =
    { 
        val ar = v()
        val dd = (for (i <- crange if p (ar(i).value); j <- 0 until ar(i).count) yield ar(i).startPos + j) (breakOut)     
        dd.toArray
    } // filterPos
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Map the elements of 'this' vector by applying the mapping function 'f'.
     *  @param f  the function to apply
     */
    def map (f: Long => Long): RleVectorL = 
    {    
        val ar = v()
        RleVectorL (for (i <- crange; j <- 0 until ar(i).count) yield f (ar(i).value))     
    } // map
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Slice 'this' vector 'from' to 'end'.
     *  @param from  the start of the slice (included)
     *  @param till  the end of the slice (excluded)
     */
    override def slice (from: Int, till: Int): VectorL = new VectorL (till - from, 
                                                          this().toArray.slice (from, till))
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Select a subset of elements of 'this' vector corresponding to a 'basis'.
     *  @param basis  the set of index positions (e.g., 0, 2, 5)
     */
    def select (basis: Array [Int]): VectoL =
    {
        val c = new VectorL (basis.length)
        for (i <- c.range) c(i) = v(basis(i)).value
        c
    } // select
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Concatenate 'this' vector and vector' b' and return the uncompressed vector
     *  @param b  the vector to be concatenated
     */
    def ++ (b: VectoL): VectoL =
    {    
        val c = new VectorL (dim + b.dim)
        for (i <- c.range) c(i) = if (i < dim) this(i) else b(i - dim)
        c     
    } // ++
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Concatenate 'this' vector and scalar 's' and return the compressed vector
     *  @param s  the scalar to be concatenated
     */
    def ++ (s: Long): RleVectorL =
    {    
        val c = new RleVectorL (dim + 1, new ReArray [TripletL] (csize))  
        for (i<- crange) c.v(i) = v(i)
        if (v(csize - 1).value == s) c.v(csize - 1).count += 1     
        else c.v(csize) = new TripletL (s, 1, v(csize - 1).startPos + v(csize - 1).count)      
        c
    } // ++
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Add 'this' vector and vector 'b'.
     *  @param b  the vector to add (any kind)
     */
    def + (b: VectoL): VectoL =
    {     
        b match {
        case _: RleVectorL => this + b.asInstanceOf [RleVectorL]
        case _             => val c = new VectorL (dim)
                              var k = 0
                              for (i <- crange; j <- 0 until v(i).count) {        
                                  c(k) = v(i).value + b(k)
                                  k += 1
                              } // for   
                              c                              
        } // match       
    } // +
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Add 'this' vector and vector 'b'.
     *  @param b  the vector to add (same kind, more efficient)
     */
    def + (x: RleVectorL): VectoL =  
    {      
        var (i, j) = (0, 0) 
        var (a, b) = (v(0).count, x.v(0).count)  
        var v1 = new ReArray [TripletL] ()
     
        var sum = v(0).value + x.v(0).value
        var currCount = MIN (a, b)
        var startPos = 0
        var idx = 0
        var currTriple = new TripletL (sum, currCount, startPos)
        v1(idx) = currTriple; idx += 1 
        var oldSum = sum
        var totalCount = currCount
          
        while (totalCount < dim) {
            if (a == currCount) { i += 1; a = v(i).count }   else a -=  currCount
            if (b == currCount) { j += 1; b = x.v(j).count } else b -= currCount
            sum = v(i).value + x.v(j).value
            currCount = MIN (a, b)
            totalCount += currCount
            startPos = currTriple.startPos + currTriple.count
            if (sum == oldSum) {
                currTriple.count += currCount
            } else {
                currTriple = new TripletL (sum, currCount, startPos)
                v1(idx) = currTriple; idx += 1 
                oldSum = sum
            } // if    
        } // while
        new RleVectorL (dim, v1)    
    } // +
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Add 'this' vector and scalar 's' and return the compressed vector
     *  @param s  the scalar to add
     */
    def + (s: Long): RleVectorL =
    {
        val c = new RleVectorL (dim, new ReArray [TripletL] (csize))
        for (i <- crange)  c.v(i) = new TripletL (v(i).value + s, v(i).count, v(i).startPos) 
        c
    } // +
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** From 'this' vector add scalar 's._2' only at position 's._1' and return 
     *  the uncompressed vector
     *  @param s  the (scalar, position) to subtract
     */
    def + (s: Tuple2 [Int, Long]): VectoL = 
    {
        val c = new VectorL (dim)
        var k = 0
        for (i <- crange) {        
            for (j <- 0 until v(i).count) {
                c(k) = if (v(i).startPos + j == s._1)  v(i).value + s._2 else v(i).value
                k += 1
            } // for
        } // for    
       c    
    } // +
                  
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Add in-place 'this' vector and vector 'b' and return the compressed 
     *  vector.
     *  @param b  the vector to add
     */
    def += (b: VectoL): VectoL =
    { 
        for (i <- range) this(i) += b(i)
        this 
    } // +=
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Add in-place 'this' vector and vector 'b' and return the compressed 
     *  vector.
     *  @param x  the vector to add
     */
    def += (x: RleVectorL) : RleVectorL = 
    {      
        var (i, j) = (0, 0) 
        var (a, b) = (v(0).count, x.v(0).count)  
        var v1 = new ReArray [TripletL] ()
     
        var sum = v(0).value + x.v(0).value
        var currCount = MIN (a, b)
        var startPos = 0
        var idx = 0
        var currTriple = new TripletL (sum, currCount, startPos)
        v1(idx) = currTriple; idx += 1
        var oldSum = sum
        var totalCount = currCount
          
        while (totalCount < dim) {
            if (a == currCount) { i += 1; a = v(i).count }   else a -=  currCount
            if (b == currCount) { j += 1; b = x.v(j).count } else b -= currCount
            sum = v(i).value + x.v(j).value
            currCount = MIN (a, b)
            totalCount += currCount
            startPos = currTriple.startPos + currTriple.count
            if (sum == oldSum) {
                currTriple.count += currCount
            } else {
                currTriple = new TripletL (sum, currCount, startPos)
                v1(idx) = currTriple; idx += 1 
                oldSum = sum
            } // if    
        } // while
        v = v1
        this     
    } // +=
         
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Add in-place 'this' vector and scalar 's' and return the compressed vector
     *  @param s  the scalar to add
     */
    def += (s: Long): RleVectorL = 
    { 
      for (i <- crange) v(i).update (v(i).value + s, v(i).count, v(i).startPos)
      this 
    } // +=
         
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the negative of 'this' vector (unary minus).
     */
    def unary_- (): RleVectorL =
    {
        val c = new RleVectorL (dim, new ReArray [TripletL] (csize))  
        for (i <- crange) c.v(i) = new TripletL (-v(i).value, v(i).count, v(i).startPos)
        c
    } // unary_-
             
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** From 'this' vector subtract vector 'b'.
     *  @param b  the vector to subtract (any kind)
     */
    def - (b: VectoL): VectoL =
    {     
        b match {
        case _: RleVectorL => this - b.asInstanceOf [RleVectorL]
        case _             => val c = new VectorL (dim)
                              var k = 0
                              for (i <- crange; j <- 0 until v(i).count) {        
                                  c(k) = v(i).value - b(k)
                                  k += 1
                              } // for   
                              c                              
        } // match       
    } // -
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** From 'this' vector subtract vector 'b'.
     *  @param b  the vector to subtract (same kind, more efficient)
     */
    def - (x: RleVectorL): VectoL =  
    {      
        var (i, j) = (0, 0) 
        var (a, b) = (v(0).count, x.v(0).count)  
        var v1 = new ReArray [TripletL] ()
     
        var diff = v(0).value - x.v(0).value
        var currCount = MIN (a, b)
        var startPos = 0
        var idx = 0
        var currTriple = new TripletL (diff, currCount, startPos)
        v1(idx) = currTriple; idx += 1 
        var oldDiff = diff
        var totalCount = currCount
          
        while (totalCount < dim) {
            if (a == currCount) { i += 1; a = v(i).count }   else a -=  currCount
            if (b == currCount) { j += 1; b = x.v(j).count } else b -= currCount
            diff = v(i).value - x.v(j).value
            currCount = MIN (a, b)
            totalCount += currCount
            startPos = currTriple.startPos + currTriple.count
            if (diff == oldDiff) {
                currTriple.count += currCount
            } else {
                currTriple = new TripletL (diff, currCount, startPos)
                v1(idx) = currTriple; idx += 1 
                oldDiff = diff
            } // if    
        } // while
        new RleVectorL (dim, v1)    
    } // -
    
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** From 'this' vector subtract scalar 's'.
     *  @param s  the scalar to subtract
     */
    def - (s: Long): RleVectorL =
    {
        val c = new RleVectorL (dim, new ReArray [TripletL] (csize))
        for (i <- crange) c.v(i) = new TripletL (v(i).value - s, v(i).count, v(i).startPos) 
        c
    } // -
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** From 'this' vector subtract scalar 's._2' only at position 's._1' and return 
     *  the uncompressed vector
     *  @param s  the (scalar, position) to subtract
     */
    def - (s: Tuple2 [Int, Long]): VectoL = 
    {
      val c = new VectorL (dim)
      var k = 0
      for (i <- crange) {        
          for (j <- 0 until v(i).count) {
              c(k) = if (v(i).startPos + j == s._1) v(i).value - s._2 else v(i).value
              k += 1
          } // for
      }  // for    
      c    
    } // -
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Subtract in-place 'this' vector and vector 'b' and return the compressed 
     *  vector.
     *  @param b  the vector to add
     */
    def -= (b: VectoL): RleVectorL = { for (i <- range) this(i) -= b(i); this }
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Subtract in-place 'this' vector and scalar 's' and return the compressed vector
     *  @param s  the scalar to add
     */
    def -= (s: Long): RleVectorL = 
    { 
      for (i <- crange) v(i).update (v(i).value - s, v(i).count, v(i).startPos)
      this 
    } // -=
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Multiply 'this' vector by vector 'b'.
     *  @param b  the vector to multiply by (any kind)
     */
    def * (b: VectoL): VectoL =
    {     
        b match {
        case _: RleVectorL => this * b.asInstanceOf [RleVectorL]
        case _             => val c = new VectorL (dim)
                              var k = 0
                              for (i <- crange; j <- 0 until v(i).count) {        
                                  c(k) = v(i).value * b(k)
                                  k += 1
                              } // for   
                              c                              
        } // match       
    } // *
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Multiply 'this' vector by vector 'x'.
     *  @param x  the vector to multiply by (same kind, more efficient)
     */
    def * (x: RleVectorL): VectoL =  
    {      
        var (i, j) = (0, 0) 
        var (a, b) = (v(0).count, x.v(0).count)  
        var v1 = new ReArray [TripletL] ()
     
        var prod = v(0).value * x.v(0).value
        var currCount = MIN (a, b)
        var startPos = 0
        var idx = 0
        var currTriple = new TripletL (prod, currCount, startPos)
        v1(idx) = currTriple; idx += 1 
        var oldProd = prod
        var totalCount = currCount
          
        while (totalCount < dim) {
            if (a == currCount) { i += 1; a = v(i).count }   else a -=  currCount
            if (b == currCount) { j += 1; b = x.v(j).count } else b -= currCount
            prod = v(i).value * x.v(j).value
            currCount = MIN (a, b)
            totalCount += currCount
            startPos = currTriple.startPos + currTriple.count
            if (prod == oldProd) {
                currTriple.count += currCount
            } else {
                currTriple = new TripletL (prod, currCount, startPos)
                v1(idx) = currTriple; idx += 1 
                oldProd = prod
            } // if    
        } // while
        new RleVectorL (dim, v1)    
    } // *
         
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Multiply 'this' vector by scalar 's'.
     *  @param s  the scalar to multiply by
     */
    def * (s: Long): RleVectorL =
    {
        val c = new RleVectorL (dim, new ReArray [TripletL] (csize))
        for (i <- crange)  c.v(i) = new TripletL (v(i).value * s, v(i).count, v(i).startPos) 
        c
    } // *
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Multiply in-place 'this' vector and vector 'b'.
     *  @param b  the vector to add
     */
    def *= (b: VectoL): RleVectorL = { for (i <- range) this(i) *= b(i); this }
   
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Multiply in-place 'this' vector and scalar 's'.
     *  @param s  the scalar to add
     */
    def *= (s: Long): RleVectorL = 
    {  
        for (i <- crange) v(i).update (v(i).value * s, v(i).count, v(i).startPos)
        this 
    } // *=
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Divide 'this' vector by vector 'b' (element-by-element).
     *  @param b  the vector to divide by (any kind)
     */
    def / (b: VectoL): VectoL =
    {     
        b match {
        case _: RleVectorL => this / b.asInstanceOf [RleVectorL]
        case _             => val c = new VectorL (dim)
                              var k = 0
                              for (i <- crange; j <- 0 until v(i).count) {        
                                  c(k) = v(i).value / b(k)
                                  k += 1
                              } // for   
                              c                              
        } // match       
    } // /
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Divide 'this' vector by  vector 'x' (element-by-element).
     *  @param x  the vector to divide by (same kind, more efficient)
     */
    def / (x: RleVectorL): VectoL =  
    {      
        var (i, j) = (0, 0) 
        var (a, b) = (v(0).count, x.v(0).count)  
        var v1 = new ReArray [TripletL] ()
     
        var div = v(0).value / x.v(0).value
        var currCount = MIN (a, b)
        var startPos = 0
        var idx = 0
        var currTriple = new TripletL (div, currCount, startPos)
        v1(idx) = currTriple; idx += 1 
        var oldDiv = div
        var totalCount = currCount
          
        while (totalCount < dim) {
            if (a == currCount) { i += 1; a = v(i).count }   else a -=  currCount
            if (b == currCount) { j += 1; b = x.v(j).count } else b -= currCount
            div = v(i).value / x.v(j).value
            currCount = MIN (a, b)
            totalCount += currCount
            startPos = currTriple.startPos + currTriple.count
            if (div == oldDiv) {
                currTriple.count += currCount
            } else {
                currTriple = new TripletL (div, currCount, startPos)
                v1(idx) = currTriple; idx += 1 
                oldDiv = div
            } // if    
        } // while
        new RleVectorL (dim, v1)    
    } // /
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Divide 'this' vector by scalar 's'.
     *  @param s  the scalar to divide by
     */
    def / (s: Long): RleVectorL =
    {
        val c = new RleVectorL (dim, new ReArray [TripletL] (csize))
        for (i <- crange)  c.v(i) = new TripletL (v(i).value / s, v(i).count, v(i).startPos) 
        c
    } // /
       
     //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Divide in-place 'this' vector and vector 'b'.
     *  @param b  the vector to add
     */
    def /= (b: VectoL): RleVectorL = { for (i <- range) this(i) /= b(i); this }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Divide in-place 'this' vector and scalar 's'.
     *  @param s  the scalar to add
     */
    def /= (s: Long): RleVectorL = 
    {  
      for (i <- crange) v(i).update (v(i).value / s, v(i).count, v(i).startPos)
      this 
    } // /=
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the vector containing each element of 'this' vector raised to the
     *  s-th power.
     *  @param s  the scalar exponent
     */
    def ~^ (s: Long): RleVectorL =
    {
        val c = new RleVectorL (dim, new ReArray [TripletL] (csize))
        for (i <- crange)  c.v(i) = new TripletL (v(i).value ~^ s, v(i).count, v(i).startPos)
        c
    } // ~^
           
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Raise each element of 'this' vector to the 's'-th power.
     *  @param s  the scalar exponent
     */
    def ~^= (s: Long): RleVectorL = 
    {  
      for (i <- crange) v(i).update (v(i).value ~^ s, v(i).count, v(i).startPos); this
    } // ~^=
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the vector containing the reciprocal of each element of 'this' vector.
     */
    def recip: RleVectorL =
    {
        val c = new RleVectorL (dim, new ReArray [TripletL] (csize))
        for (i <- crange) c.v(i) = new TripletL (1l / v(i).value, v(i).count, v(i).startPos)
        c
    } // recip
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the vector that is the element-wise absolute value of 'this' vector.
     */
    def abs: RleVectorL =
    {    
        val c = new RleVectorL (dim, new ReArray [TripletL] (csize))
        var oldTriple = v(0)
        var oldValue = ABS (oldTriple.value)
        var idx = 0
        c.v(idx) = oldTriple; idx += 1
        for (i <- 1 until v.length) 
            {          
                 val currValue = ABS (v(i).value)
                 if (currValue == oldValue) {
                     oldTriple.update (oldTriple.value, oldTriple.count + v(i).count, oldTriple.startPos)
                 } else {
                     oldTriple = new TripletL (currValue, v(i).count, c.v(i-1).startPos + c.v(i-1).count)
                     oldValue = oldTriple.value
                     c.v(idx) = oldTriple; idx += 1
                 } // if
            } // for
        c
    } // abs
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Get the sum of the elements of 'this' vector.
     */
    def sum: Long = v.foldLeft (0l) ((s, a) => s + (a.value * a.count)) 
          
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Sum the absolute value of the elements of 'this' vector.
     */
    def sumAbs: Long = v.foldLeft (0l) ((s, a) => s + (ABS (a.value) * a.count)) 
        
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Sum the elements of 'this' vector skipping the 'i'-th element (Not Equal 'i').
     *  @param i  the index of the element to skip
     */
    def sumNE (i: Int): Long = sum - this(i)
     
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Sum the positive (> 0) elements of 'this' vector.
     */    
    def sumPos: Long = v.foldLeft (0l) ((s, a) => s + MAX ((a.value * a.count), 0l))
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Establish the rank order of the elements in 'self' vector, e.g.,
     *  (8.0, 2.0, 4.0, 6.0) is (3, 0, 1, 2).
     */
    def rank: VectorI = VectorI (iqsort (this().toArray))
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Cumulate the values of 'this' vector from left to right (e.g., create a
     *  CDF from a pmf).  Example: (4, 2, 3, 1) --> (4, 6, 9, 10)
     */
    def cumulate: VectoL =
    {        
        val c = new VectorL (dim)
        var sum = 0l
        var k = 0
        for (i <- crange) {        
            for (j <- 0 until v(i).count) {
                sum += v(i).value
                c(k) = sum
                k += 1
            } // for
        }  // for    
       c   
    } // cumulate
      
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Normalize 'this' vector so that it sums to one (like a probability vector).
     */
    def normalize: RleVectorL = this * (1l / sum)

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Normalize 'this' vector so its length is one (unit vector).
     */
    def normalizeU: RleVectorL = this * (1l / norm)

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Normalize 'this' vector to have a maximum of one.
     */
    def normalize1: RleVectorL = this * (1l / max ())
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the dot product (or inner product) of 'this' vector with vector 'b'.
     *  @param b  the other vector
     */
    def dot (b: VectoL): Long =
    {
     
      b match {
      case _: RleVectorL => this dot b.asInstanceOf [RleVectorL]
      case _             => var sum = 0l
                            for (i <- range) sum += this(i) * b(i)
                            sum 
      } // match     
    } // dot
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the dot product (or inner product) of 'this' compressed vector 
     *  with compressed vector 'b'.
     *  @param b  the other vector
     */
    def dot (x: RleVectorL): Long =
    {
        var (i, j)     =  (0, 0)
        var (a, b)     =  (v(0).count, x.v(0).count)
        var currCount  =  MIN (a, b)
        var sum        =  (v(0).value * x.v(0).value) * currCount
        var totalCount = currCount
        while (totalCount < dim) {
            if (a == currCount) { i += 1; a = v(i).count   } else a -= currCount
            if (b == currCount) { j += 1; b = x.v(j).count } else b -= currCount
            currCount = MIN (a, b)
            sum += (v(i).value * x.v(j).value) * currCount
            totalCount += currCount
        } // while
        sum
    } // dot
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the Euclidean norm (2-norm) squared of 'this' vector.
     */
    override def normSq: Long = 
    {
        var sum = 0l
        for (i <- crange) {
           val x = v(i).value
           sum += (x * x) * v(i).count 
        } // for
        sum
    } // normSq
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the Manhattan norm (1-norm) of 'this' vector.
     */  
    def norm1: Long = v.foldLeft (0l) ((s, a) => s + (ABS (a.value) * a.count))
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the maximum element in 'this' compressed vector.
     *  @param e  the ending index (exclusive) for the search
     */
    def max (e: Int = dim): Long =
    {
        if (e <= 0) flaw ("max", "the ending index e can't be negative")
        var x = v(0).value
        for (i <- 1 until csize if v(i).startPos < e && v(i).value > x) x = v(i).value
        x
    } // max
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Take the maximum of 'this' vector with vector 'b' (element-by element).
     *  I.e., the 'dim' of this Vector and 'b' should be same.
     *  @param b  the other vector
     */
    def max (b: VectoL): VectoL =
    {
        val c = new VectorL (dim)
        var k = 0
        for (i <- crange) {
            for (j <- 0 until v(i).count) {
                c(k) = if (b(k) > v(i).value) b(k) else v(i).value
                k += 1
            } // for
        } // for
        c
    } // max
      
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the minimum element in 'this' vector.
     *  @param e  the ending index (exclusive) for the search
     */
    def min (e: Int = dim): Long =
    {
        if (e <= 0) flaw ("min", "the ending index e can't be negative")
        var x = v(0).value
        for (i <- 1 until csize if v(i).startPos < e && v(i).value < x) x = v(i).value
        x
    } // min
      
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Take the minimum of 'this' vector with vector 'b' (element-by element).
     *  @param b  the other vector
     */
    def min (b: VectoL): VectoL =
    {
        val c = new VectorL (dim)
        var k = 0
        for (i <- crange) {
            for (j <- 0 until v(i).count) {
                c(k) = if (b(k) < v(i).value) b(k) else v(i).value
                k += 1
            } // for
        }  // for
        c
    } // min
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the argument maximum of 'this' vector (index of maximum element).
     *  @param e  the ending index (exclusive) for the search
     */
    def argmax (e: Int = dim): Int =
    {
        var j = 0
        var temp = v(0).value
        for (i <- 1 until csize if v(i).startPos + v(i).count - 1 < e && v(i).value > temp) {
            temp = v(i).value; j = v(i).startPos
        } // for
        j
    } // argmax
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find the argument minimum of 'this' vector (index of minimum element).
     *  @param e  the ending index (exclusive) for the search
     */
    def argmin (e: Int = dim): Int =
    {
        var j = 0
        var temp = v(0).value
        for (i <- 1 until csize if v(i).startPos + v(i).count - 1 < e && v(i).value < temp) {
            temp = v(i).value; j = v(i).startPos
        } // for
        j
    } // argmin
     
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the argument minimum of 'this' vector (-1 if its not negative).
     *  @param e  the ending index (exclusive) for the search
     */
    def argminNeg (e: Int = dim): Int =
    {
        var j = 0
        var temp = v(0).value
        for (i <- 1 until csize if v(i).startPos + v(i).count - 1 < e && v(i).value < temp) {
            temp = v(i).value; j = v(i).startPos
        } // for
        if (temp < 0l) j else -1
    } // argminNeg
  
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the argument maximum of 'this' vector (-1 if its not positive).
     *  @param e  the ending index (exclusive) for the search
     */
    def argmaxPos (e: Int = dim): Int =
    {
        var j = 0
        var temp = v(0).value
        for (i <- 1 until csize if v(i).startPos + v(i).count - 1 < e && v(i).value > temp) {
            temp = v(i).value; j = v(i).startPos
        } // for
        if (temp > 0l) j else -1
    } // argmaxPos
        
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the index of the first negative element in 'this' vector (-1 otherwise).
     *  @param e  the ending index (exclusive) for the search
     */
    def firstNeg (e: Int = dim): Int =
    {
        for (i <- crange if v(i).startPos + v(i).count - 1 < e && v(i).value < 0l) {
            return v(i).startPos
        } // for
        -1
    } // firstNeg
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the index of the first positive element in 'this' vector (-1 otherwise).
     *  @param e  the ending index (exclusive) for the search
     */
    def firstPos (e: Int = dim): Int =
    {
        for (i <- crange if v(i).startPos + v(i).count - 1 < e && v(i).value > 0l) {
            return v(i).startPos
        } // for
        -1
    } // firstPos
       
     //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the index of the first occurrence of element 'x' in 'this' vector,
     *  or -1 if not found.
     *  @param x  the given element
     *  @param e  the ending index (exclusive) for the search
     */
    def indexOf (x: Long, e: Int = dim): Int =
    {
        for (i <- crange if v(i).startPos + v(i).count - 1 < e &&  v(i).value == x) {
            return v(i).startPos
        } // for
        -1
    } // indexOf
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Find and return index of first element satisfying predicate 'p', or
     *  -1 if not found.
     *  @param p  the predicate to check
     */
    def indexWhere (p: (Long) => Boolean): Int = 
    {      
      for (i <- crange if p (v(i).value)) return v(i).startPos
     -1     
    } // indexWhere
    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count the number of strictly negative elements in 'this' vector.
     */
    def countNeg: Int =
    {
        var count = 0
        for (i <- crange if v(i).value < 0l) count += v(i).count
        count
    } // countNeg
    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count the number of strictly positive elements in 'this' vector.
     */
    def countPos: Int =
    {
        var count = 0
        for (i <- crange if v(i).value > 0l) count += v(i).count
        count
    } // countPos
        
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count the number of distinct elements in 'this' vector.
     */
    def distinct: RleVectorL =
    {     
        val u = new ReArray [TripletL] ()
        var idx = 0
        val s = HashSet [Long] ()  
        var sPos = 0
        for (t <- v if !(s contains t.value)) {
            s += t.value  
            u(idx) = new TripletL (t.value, 1, sPos); idx += 1
            sPos += 1
        } // for
        new RleVectorL (u.length, u)
    } // distinct
        
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count the number of distinct elements in 'this' vector.
     */
    def countinct: Int =
    {
        val s = HashSet [Long] ()  
        for (t <- v) s += t.value
        s.size
    } // countinct
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Determine whether 'x' is contained in 'this' vector.
     *  @param x  the element to be checked
     */
    def contains (x: Long): Boolean =
    {
        for (i <- crange if v(i).value == x) true
        false
    } // contains
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reverse the order of the elements in 'this' vector.
     */
    def reverse (): RleVectorL = 
    {      
        val nv = new ReArray [TripletL] (csize)
        for (i <-crange.reverse) nv(i) = v(i) 
        new RleVectorL (dim, nv)       
    } // reverse
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Determine whether 'this' vector is in sorted (ascending) order.
     */
    def isSorted: Boolean = (new SortingL (this().toArray)).isSorted 
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Sort 'this' vector in-place in ascending (non-decreasing) order.
     */
    def sort () 
     { 
         throw new NoSuchMethodException ("sort is not implemented yet")
     } // sort

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Sort 'this' vector in-place in descending (non-increasing) order.
     */
    def sort2 () 
    {
        throw new NoSuchMethodException ("sort2 is not implemented yet") 
    } // sort2 
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Swap elements 'i' and 'j' in 'this' vector. 
     *  @param i  the first element in the swap
     *  @param j  the second element in the swap
     */
    def swap (i: Int, j: Int)
    {
        val t = this(j); this(j) = this(i); this(i) = t
    } // swap
       
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Swap elements 'i' and 'j' in 'this' vector. (2nd implementation)
     *  @param i  the first element in the swap
     *  @param j  the second element in the swap
     */
    def swap2 (i: Int, j: Int)
    {        
        val c = toDense
        val t = c(j); c(j) = c(i); c(i) = t
        val x = RleVectorL (c)
        v = x.v        
    } // swap2
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Check whether 'this' vector is nonnegative (has no negative elements).
     */
    def isNonnegative: Boolean =
    {
        for (i <- crange if v(i).value < 0l) return false
        true
    } // isNonnegative
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compare 'this' vector with vector 'b'.
     *  @param b  the other vector
     */
    def tryCompareTo [B >: RleVectorL] (b: B)
        (implicit view_1: (B) => PartiallyOrdered [B]): Option [Int] =
    {
        throw new NoSuchMethodException ("tryCompareTo is not implemented yet") 
    } // tryCompareTo
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Converts the compressed vector to a dense vector.
     *  @return Vector of base type `Long` containing uncompressed data
     */
    def toDense: VectorL =
    { 
        val ar = v()
        VectorL (for (i <- 0 until v.length; j <- 0 until ar(i).count) yield ar(i).value)
    } // toDense
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Override equals to determine whether 'this' vector equals vector 'b.
     *  @param b  the vector to compare with this
     */
    override def equals (b: Any): Boolean =
    {
        b.isInstanceOf [RleVectorL] && (v == b.asInstanceOf [RleVectorL].v)
    } // equals
          
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Must also override hashCode for 'this' vector to be compatible with equals.
     */
    override def hashCode (): Int = v.hashCode
   
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Convert 'this' vector to a string.
     */
    override def toString: String = 
    {   
        var sb = new StringBuilder ("RleVectorL (")
        if (dim == 0) return sb.append (")").mkString
        for (i <- crange) { sb.append ("("); sb.append (fString.format (v(i).value))
                            sb.append (v(i).count); sb.append (",")
                            sb.append (v(i).startPos); sb.append ("), ")
        } // for
        sb.replace (sb.length-2, sb.length-1, ")").mkString
    } // toString   
     
} // RleVectorL class


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `RleVectorLTest` object tests the operations provided by `RleVectorL`.
 *  > run-main scalation.linalgebra.RleVectorLTest
 */
object RleVectorLTest extends App 
{
    val a = VectorL (10)
    val c = RleVectorL (a)
    println ("c = " + c)
    println ("c.dim = " + c.dim)
    println ("c.sum = " + c.sum)
    println ("c.max () = " + c.max ())
    println ("c.toDense = " + c.toDense)

} // RleVectorLTest object


