
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller, Hao Peng, Zhe Jin
 *  @version 1.2
 *  @date    Mon Jul 27 01:27:00 EDT 2015
 *  @see     LICENSE (MIT style license file).
 */

package scalation.analytics.classifier

import scala.collection.mutable.{Map, Set => SET}
import scala.util.control.Breaks.{break, breakable}

import scalation.graphalytics.Pair
import scalation.graphalytics.mutable.{MGraph, MinSpanningTree}
import scalation.linalgebra.{MatrixD, MatriI, MatrixI, VectorD, VectoI, VectorI}
import scalation.linalgebra.gen.{HMatrix3, HMatrix4, HMatrix5}
import scalation.relalgebra.Relation

import BayesClassifier.me_default

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `SelTAN` class implements an Integer-Based Tree Augmented Selective
 *  Naive Bayes Classifier,  which is a combinations of two commonly used classifiers
 *  for discrete input data.  The classifier is trained using a data matrix 'x' and a
 *  classification vector 'y'.  Each data vector in the matrix is classified into one
 *  of 'k' classes numbered 0, ..., k-1.  Prior probabilities are calculated based on
 *  the population of each class in the training-set.  Relative posterior probabilities
 *  are computed by multiplying these by values computed using conditional probabilities.
 *  The classifier supports limited dependency between features/variables. The classifier
 *  also uses backward elimination algorithm in an attempt to find the most important
 *  subset of features/variables.
 *  -----------------------------------------------------------------------------
 *  @param x      the integer-valued data vectors stored as rows of a matrix
 *  @param y      the class vector, where y(l) = class for row l the matrix x, x(l)
 *  @param fn     the names for all features/variables
 *  @param k      the number of classes
 *  @param cn     the names for all classes
 *  @param fset   the `Boolean` array indicating the selected features
 *  @param vc     the value count (number of distinct values) for each feature
 *  @param me     use m-estimates (me == 0 => regular MLE estimates)
 *  @param thres  the correlation threshold between 2 features for possible parent-child relationship
 */
class SelTAN (x: MatriI, y: VectoI, fn: Array [String], k: Int, cn: Array [String], private var fset: Array [Boolean] = null,
              thres: Double = 0.3, me: Int = me_default, private var vc: VectoI = null)
      extends BayesClassifier (x, y, fn, k, cn)
{
    private val DEBUG  = false                          // debug flag
    private val TOL    = 0.01                           // tolerance indicating negligible improvement adding features
    private val cor    = calcCorrelation                // feature correlation matrix
    private var parent = new VectorI (n)                // allocate the parent vector
    private val vcp    = new VectorI (n)                // value count for the parent

    private val popC  = new VectorI (k)                 // frequency counts for classes 0, ..., k-1
    private val probC = new VectorD (k)                 // probabilities for classes 0, ..., k-1
    private val popX  = new HMatrix4 [Int] (k, n)       // conditional frequency counts for variable/feature j
    private val probX = new HMatrix4 [Double] (k, n)    // conditional probabilities for variable/feature j

    if (vc == null) vc = vc_fromData                    // set to default for binary data (2)
    if (fset == null) fset = Array.fill (n)(true)       // set to default, all features included
    computeParent ()                                    // initialize the parent of each feature
    computeVcp ()                                       // initialize the value count of each parent feature

    popX.alloc (fset, vc, vcp)
    probX.alloc (fset, vc, vcp)

    if (DEBUG) {
        println ("feature set fset   = " + fset.deep)
        println ("parents parent        = " + parent)
        println ("value count vc     = " + vc)
        println ("value count vcp    = " + vcp)
        println ("correlation matrix = " + cor)
    } // if

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the parent of each feature based on the correlation matrix.
     *  Feature x_i is only a possible candidate for parent of feature x_j if i < j.
     */
    def computeParent ()
    {
        val countC   = new VectorD (k)                                           // count the number of classes in each class
        val countXYC = new HMatrix5 [Double] (k, n, n, vc.toArray, vc.toArray)   // countXYC count the number where X=x,Y=y,C=c
        val countXC  = new HMatrix3 [Double] (k, n, vc.toArray)                  // countXC count the number where X=x,C=c
        val ch       = Array.ofDim [SET [Int]] (n)
        val elabel   = Map [Pair, Double] ()
//      parent(0)    = -1                                                        // feature 0 does not have a parent

        //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        /** Compute frequency counts for each value in each variable.
         */
        def frequencies ()
        {
            for (i <- 0 until m) {
                countC(y(i)) += 1
                for (j <- 0 until n) {
                    countXC(y(i), j, x(i, j)) += 1
                    for (t <- 1 until n) countXYC(y(i), j, t, x(i, j), x(i, t)) += 1
                } // for
            } // for
        } // frequencies

        //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        /** Compute frequency counts for each value in each variable.
         */
        def probabilities: (VectorD, HMatrix3 [Double], HMatrix5 [Double]) =
           {
            for (i <- 0 until k; j <- 0 until n; t <- 0 until vc(j)) countXC(i, j, t) = (countXC(i, j, t) + (me.toDouble) / m) / (m + me)
            for (i <- 0 until n; j <- 0 until n; t <- 0 until k; p <- 0 until vc(i); q <- 0 until vc(j)) {
                countXYC(t, i, j, p, q) = (countXYC(t, i, j, p, q) + (me.toDouble / m)) / (m + me)
            } // for
            (countC / m, countXC, countXYC)
        } // frequencies

        //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        /** Create MinSpanningTree from conditional mutual information
         */
        def minSpanningTree (ip: MatrixD): MinSpanningTree =
        {
            for (i <- 0 until n) ch(i) = SET((i + 1 until n): _*)
            for (i <- 0 until n; j <- i + 1 until n) elabel += new Pair(i, j) -> ip(i, j)
            val g = new MGraph (ch, Array.ofDim(n), elabel)
            new MinSpanningTree(g, false, false)
        } // minSpanningTree

        frequencies ()
        val probs = probabilities
        val ipxyz = condMutualInformation (probs._1, probs._2, probs._3)
        parent    = VectorI (minSpanningTree(ipxyz).makeITree ())

        //println("temppar=   "+ temppar.deep)
        //for(i<- 0 until n) parent(i) = temppar(i)

    } // computeParent

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Compute the value count of each parent feature based on the parent vector.
     */
    def computeVcp ()
    {
        vcp.set (1)                                 //set default value count to 1
        for (j <- 0 until n if (fset(j) && parent(j) > -1)) vcp(j) = vc(parent(j))
    } // computeVcp

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count the frequencies for 'y' having class 'i' and 'x' for cases 0, 1, ...
     *  @param testStart  starting index of test region (inclusive)
     *  @param testEnd    ending index of test region (exclusive)
     */
    private def frequencies (testStart: Int, testEnd: Int)
    {
        for (l <- 0 until m if l < testStart || l >= testEnd) {
            // l = lth row of data matrix x
            val i = y(l)                                    // get the class
            popC(i) += 1                                    // increment ith class
            for (j <- 0 until n if fset(j)) {
                if (parent(j) > -1) popX(i, j, x(l, j), x(l, parent(j))) += 1
                else popX(i, j, x(l, j), 0) += 1
            } // for
        } // for

        if (DEBUG) {
            println ("popC = " + popC)                      // #(C = i)
            println ("popX = " + popX)                      // #(X_j = x & C = i)
        } // if
    } // frequencies

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the classifier by computing the probabilities for C, and the
     *  conditional probabilities for X_j.
     *  @param testStart  starting index of test region (inclusive)
     *  @param testEnd    ending index of test region (exclusive)
     */
    def train (testStart: Int = 0, testEnd: Int = 0)
    {
        frequencies (testStart, testEnd)                    // compute frequencies skipping test region

        for (i <- 0 until k) {                              // for each class i
            val pci = popC(i).toDouble                      // population of class i
            probC(i) = pci / md                             // probability of class i

            for (j <- 0 until n if fset(j)) {               // for each feature j in fset
                val me_vc = me / vc(j).toDouble
                for (xj <- 0 until vc(j); xp <- 0 until vcp(j)) {
                    // for each value for feature j: xj, parent(j): xp
                    probX(i, j, xj, xp) = (popX(i, j, xj, xp) + me_vc) / (pci + me)
                } // for
            } // for
        } // for

        if (DEBUG) {
            println("probC = " + probC)                     // P(C = i)
            println("probX = " + probX)                     // P(X_j = x | C = i)
        } // if
    } // train

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count the frequencies for 'y' having class 'i' and 'x' for cases 0, 1, ...
     *  @param itrain  indices of the instances considered train data
     */
    private def frequencies (itrain: Array [Int])
    {
        for (l <- itrain) {                                  // l = lth row of data matrix x
            val i = y(l)                                     // get the class
            popC(i) += 1                                     // increment ith class
            for (j <- 0 until n if fset(j)) {
                if (parent(j) > -1) popX(i, j, x(l, j), x(l, parent(j))) += 1
                else popX(i, j, x(l, j), 0) += 1
            } // for
        } // for

        if (DEBUG) {
            println ("popC = " + popC)                       // #(C = i)
            println ("popX = " + popX)                       // #(X_j = x & C = i)
        } // if
    } // frequencies

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the classifier by computing the probabilities for C, and the
     *  conditional probabilities for X_j.
     *  @param itrain  indices of the instances considered train data
     */
    override def train (itrain: Array [Int])
    {
        frequencies (itrain)                                 // compute frequencies skipping test region

        for (i <- 0 until k) {                               // for each class i
            val pci = popC(i).toDouble                       // population of class i
            probC(i) = pci / md                              // probability of class i

            for (j <- 0 until n if fset(j)) {                // for each feature j in fset
                val me_vc = me / vc(j).toDouble
                for (xj <- 0 until vc(j); xp <- 0 until vcp(j)) {
                    // for each value for feature j: xj, parent(j): xp
                    probX(i, j, xj, xp) = (popX(i, j, xj, xp) + me_vc) / (pci + me)
                } // for
            } // for
        } // for

        if (DEBUG) {
            println ("probC = " + probC)                     // P(C = i)
            println ("probX = " + probX)                     // P(X_j = x | C = i)
        } // if
    } // train


    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a discrete data vector 'z', classify it returning the class number
     *  (0, ..., k-1) with the highest relative posterior probability.
     *  Return the best class, its name and its realtive probability.
     *  @param z  the data vector to classify
     */
    def classify (z: VectoI): (Int, String, Double) =
    {
        val prob = new VectorD(k)
        for (i <- 0 until k) {
            prob(i) = probC(i) // P(C = i)                                   // P(C = i)
            for (j <- 0 until n if fset(j)) {
                prob(i) *= (if (parent(j) > -1) probX(i, j, z(j), z(parent(j))) // P(X_j = z_j | C = i), parent
                else probX(i, j, z(j), 0)) // P(X_j = z_j | C = i), no parent
            } //for
        } // for
        if (DEBUG) println ("prob = " + prob)
        val best = prob.argmax ()             // class with the highest relative posterior probability
        (best, cn(best), prob(best))          // return the best class, its name and its probability
    } // classify

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reset or re-initialize the frequency tables and the probability tables
     *  with the updated parent vector.
     */
    def reset ()
    {
        computeParent ()
        computeVcp ()
        popC.set (0)
        probC.set (0)
        popX.clear ()
        probX.clear ()
        popX.alloc (fset, vc, vcp)
        probX.alloc (fset, vc, vcp)
    } // reset

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Build the Tree Augmented Selective Naive Bayes classier model by using backward-elimination
     *  Selective algorithm. Limited dependencies between variables/features are also supported.
     *  @param testStart  starting index of test region (inclusive)
     *  @param testEnd    ending index of test region (exclusive)
     */
    def buildModel (testStart: Int = 0, testEnd: Int = 0): (Array [Boolean], DAG) =
    {
        for (j <- 0 until n) fset(j) = true           // set the feature set to all features included
        //initialize the model using n-fold cross validation and obtaining the accuracy without removing any features
        var accuracy = crossValidateRand ()
        if (DEBUG) println ("Initial accuracy with no feature removed: " + accuracy)

        // keep removing one feature at a time until no more feature should be removed
        breakable { while (true) {
            var accuracyDiff = 0.0
            var minDiff      = 1.0
            var toRemove     = 0
            if (DEBUG) println ("Try to removing each feature and achieve best accuracy...")

            for (j <- 0 until n if fset(j)) {
                if (DEBUG) println("Test by temporarily removing feature " + j)
                fset(j) = false
                accuracyDiff = accuracy - crossValidateRand()
                if (accuracyDiff <= minDiff) { minDiff = accuracyDiff; toRemove = j }
                fset(j) = true
            } // for
            accuracy -= minDiff

            //only remove the feature if the minimum accuracy drop is less than a small TOL value (acceptable accuracy reduction)
            if (fset(toRemove) && minDiff < TOL) {
                if (DEBUG) println ("Feature " + toRemove + " has been removed from the model.")
                fset(toRemove) = false
                if (DEBUG) println ("Re-train model by removing feature " + toRemove)
                crossValidateRand()
                if (DEBUG) println ("The new accuracy is " + accuracy + " after removing feature " + toRemove)
            } else {
                if (DEBUG) println ("No more features to removed: Re-train the model without removing any features")
                crossValidateRand ()
                if (DEBUG) {
                    println ("Final parent  = " + parent)
                    println ("Final fset = " + fset.deep)
                } // if
                break
            } // if
        }} // while
        computeParent ()
        val pp: Traversable [Array [Int]] = for (p <- parent) yield Array (p)
        (fset, new DAG(pp.toArray))
    } // buildModel class

} // SelTAN class


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** `SelTAN` is the companion object for the `SelTAN` class.
 */
object SelTAN
{
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a `SelTAN object, passing 'x' and 'y' together in one table.
     *  @param xy     the data vectors along with their classifications stored as rows of a matrix
     *  @param fn     the names of the features/variables
     *  @param k      the number of classes
     *  @param cn     the names for all classes
     *  @param fset   the `Boolean` array indicating the selected features
     *  @param vc     the value count (number of distinct values) for each feature
     *  @param me     use m-estimates (me == 0 => regular MLE estimates)
     *  @param thres  the correlation threshold between 2 features for possible parent-child relationship
     */
    def apply (xy: MatriI, fn: Array [String], k: Int, cn: Array [String],
               fset: Array [Boolean] = null, thres: Double = 0.3, me: Int = me_default, vc: VectoI = null) =
    {
        new SelTAN (xy(0 until xy.dim1, 0 until xy.dim2 - 1), xy.col(xy.dim2 - 1), fn, k, cn,
                    fset, thres, me, vc)
    } // apply

} // SelTAN object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `SelTANTest` object is used to test the `SelTAN` class.
 *  Classify whether a car is more likely to be stolen (1) or not (1).
 *  @see www.inf.u-szeged.hu/~ormandi/ai2/06-SelTAN-example.pdf
 *  > run-main scalation.analytics.classifier.SelTANTest
 */
object SelTANTest extends App
{
    // x0: Color:   Red (1), Yellow (0)
    // x1: Type:    SUV (1), Sports (0)
    // x2: Origin:  Domestic (1), Imported (0)
    // features:                x0 x1 x2
    val x = new MatrixI((10, 3), 1, 0, 1,               // data matrix
                                 1, 0, 1,
                                 1, 0, 1,
                                 0, 0, 1,
                                 0, 0, 0,
                                 0, 1, 0,
                                 0, 1, 0,
                                 0, 1, 1,
                                 1, 1, 0,
                                 1, 0, 0)

    val y = VectorI (1, 0, 1, 0, 1, 0, 1, 0, 0, 1)      // classification vector: 0(No), 1(Yes))
    val fn = Array("Color", "Type", "Origin")           // feature/variable names
    val cn = Array("No", "Yes")                         // class names

    println("xy = " + (x :^+ y))
    println("---------------------------------------------------------------")

    val asnb = new SelTAN (x, y, fn, 2, cn)             // create the classifier

    // train the classifier ---------------------------------------------------
    // asnb.train ()
    asnb.buildModel (3)

    // test sample ------------------------------------------------------------
    val z1 = VectorI (1, 0, 1)                         // new data vector to classify
    val z2 = VectorI (1, 1, 1)                         // new data vector to classify
    println ("classify (" + z1 + ") = " + asnb.classify (z1) + "\n")
    println ("classify (" + z2 + ") = " + asnb.classify (z2) + "\n")

    asnb.crossValidate ()                              // cross validate the classifier

} // SelTANTest object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `SelTANTest2` object is used to test the `SelTAN` class.
 *  Given whether a person is Fast and/or Strong, classify them as making C = 1
 *  or not making C = 0 the football team.
 *  > run-main scalation.analytics.classifier.SelTANTest2
 */
object SelTANTest2 extends App
{
    // training-set -----------------------------------------------------------
    // x0: Fast
    // x1: Strong
    // y:  Classification (No/0, Yes/1)
    // features:                  x0 x1  y
    val xy = new MatrixI((10, 3), 1, 1, 1,
                                  1, 1, 1,
                                  1, 0, 1,
                                  1, 0, 1,
                                  1, 0, 0,
                                  0, 1, 0,
                                  0, 1, 0,
                                  0, 1, 1,
                                  0, 0, 0,
                                  0, 0, 0)

    val fn = Array ("Fast", "Strong")
    val cn = Array ("No", "Yes")

    println ("xy = " + xy)
    println ("---------------------------------------------------------------")

    val asnb = SelTAN (xy, fn, 2, cn)                  // create the classifier

    // train the classifier ---------------------------------------------------
    // asnb.train ()
    asnb.buildModel ()

    // test sample ------------------------------------------------------------
    val z = VectorI (1, 0)                             // new data vector to classify
    println ("classify (" + z + ") = " + asnb.classify (z) + "\n")

    asnb.crossValidate()                               // cross validate the classifier

} // SelTANTest2 object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `SelTanTest3` object is used to test the `SelTAN` class
 *  > run-main scalation.analytic.classifier.SelTANTest3
 */
object SelTANTest3 extends App
{
    val filename = BASE_DIR + "breast-cancer.arff"
    var data = Relation (filename, -1, null)
    val xy = data.toMatriI2 (null)
    val fn = data.colName.toArray
    val cn = Array ("0", "1")                          // class names
    val k  = 2

    println("---------------------------------------------------------------")
    val anb = SelTAN (xy, fn, k, cn)                   // create the classifier
    anb.buildModel ()
    anb.train ()
    anb.crossValidate ()

} // SelTANTest3 object

