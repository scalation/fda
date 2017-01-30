
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller, Zhe Jin
 *  @version 1.2
 *  @date    Sat Sep  8 13:53:16 EDT 2012
 *  @see     LICENSE (MIT style license file).
 *
 *  @see eric.univ-lyon2.fr/~ricco/tanagra/fichiers/en_Tanagra_Naive_Bayes_Classifier_Explained.pdf
 */

package scalation.analytics.classifier

import scala.math.log

import scalation.linalgebra.{MatriI, VectoI, VectorD, VectorI}
import scalation.linalgebra.gen.HMatrix3
import scalation.relalgebra.Relation
import scalation.util.{banner, time}

import BayesClassifier.me_default

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `NaiveBayes` class implements an Integer-Based Naive Bayes Classifier,
 *  which is a commonly used such classifier for discrete input data.  The
 *  classifier is trained using a data matrix 'x' and a classification vector 'y'.
 *  Each data vector in the matrix is classified into one of 'k' classes numbered
 *  0, ..., k-1.  Prior probabilities are calculated based on the population of
 *  each class in the training-set.  Relative posterior probabilities are computed
 *  by multiplying these by values computed using conditional probabilities.  The
 *  classifier is naive, because it assumes feature independence and therefore
 *  simply multiplies the conditional probabilities.
 *  -----------------------------------------------------------------------------
 *  @param x   the integer-valued data vectors stored as rows of a matrix
 *  @param y   the class vector, where y(l) = class for row l of the matrix x, x(l)
 *  @param fn  the names for all features/variables
 *  @param k   the number of classes
 *  @param cn  the names for all classes
 *  @param vc  the value count (number of distinct values) for each feature
 *  @param me  use m-estimates (me == 0 => regular MLE estimates)
 */
class NaiveBayes (x: MatriI, y: VectoI, fn: Array [String], k: Int, cn: Array [String],
                  private var vc: VectoI = null, me: Int = me_default)
      extends BayesClassifier (x, y, fn, k, cn)
{
    private val DEBUG = true                               // debug flag
    private val cor   = calcCorrelation                    // feature correlation matrix

    private val f_Y  = new VectorI (k)                     // frequency counts for class yi
    private val f_YX = new HMatrix3 [Int] (k, n)           // frequency counts for class yi & feature xj

    private var p_Y: VectorD = _                           // prior probabilities for class yi
    private val p_X_Y = new HMatrix3 [Double] (k, n)       // conditional probabilities for feature xj given class yi

    if (vc == null) {
        shiftToZero; vc = vc_fromData                      // set value counts from data
    } // if
    f_YX.alloc (vc().toArray)
    p_X_Y.alloc (vc().toArray)

    if (DEBUG) {
        println ("distinct value count vc = " + vc)
        println ("correlation matrix      = " + cor)
    } // if

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Build a model.
     *  @param testStart  starting index of test region (inclusive) used in cross-validation
     *  @param testEnd    ending index of test region (exclusive) used in cross-validation
     */
    def buildModel (testStart: Int, testEnd: Int): (Array [Boolean], DAG) =
    {
        (Array.fill (n)(true), new DAG (Array.ofDim [Int] (n, 0)))
    } // buildModel

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Increment frequency counters based on the 'i'th row of the data matrix.
     *  @param i  the index for current data row
     */
    private def increment (i: Int)
    {
        val yi   = y(i)                                       // get the class for ith row
        f_Y(yi) += 1                                          // increment frequency for class yi
        for (j <- 0 until n) f_YX(yi, j, x(i, j)) += 1        // increment frequency for class yi and all X-features
    } // increment

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count the frequencies for 'y' having class 'yi' and value 'x' for cases 0, 1, ...
     *  Only the test region from 'testStart' to 'testEnd' is skipped, the rest is
     *  training data.
     *  @param testStart  starting index of test region (inclusive) used in cross-validation
     *  @param testEnd    ending index of test region (exclusive) used in cross-validation
     */
    private def frequencies (testStart: Int, testEnd: Int)
    {
        if (DEBUG) banner ("frequencies (testStart, testEnd)")
        for (i <- 0 until m if i < testStart || i >= testEnd) increment (i)

        if (DEBUG) {
            println ("f_Y  = " + f_Y)                         // #(Y = yi)
            println ("f_YX = " + f_YX)                        // #(Y = yi & X_j = x)
        } // if
    } // frequencies

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Count the frequencies for 'y' having class 'yi' and value 'x' for cases 0, 1, ...
     *  Only the test region from 'testStart' to 'testEnd' is skipped, the rest is
     *  training data.
     *  @param itrain  indices of the instances considered train data
     */
    private def frequencies (itrain: Array [Int])
    {
        if (DEBUG) banner ("frequencies (itrain)")
        for (i <- itrain) increment (i)

        if (DEBUG) {
            println ("f_Y  = " + f_Y)                         // #(Y = yi)
            println ("f_YX = " + f_YX)                        // #(Y = yi & X_j = x)
        } // if
    } // frequencies

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the classifier by computing the probabilities for Y, and the
     *  conditional probabilities for X_j.
     *  @param testStart  starting index of test region (inclusive) used in cross-validation.
     *  @param testEnd    ending index of test region (exclusive) used in cross-validation.
     */
    def train (testStart: Int, testEnd: Int)
    {
        frequencies (testStart, testEnd)                      // compute frequencies skipping test region
        if (DEBUG) banner ("train (testStart, testEnd)")
        train2 ()
    } // train

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the classifier by computing the probabilities for Y, and the
     *  conditional probabilities for X_j.
     *  @param itrain indices of the instances considered train data
     */
    override def train (itrain: Array [Int])
    {
        frequencies (itrain)                                  // compute frequencies skipping test region
        if (DEBUG) banner ("train (itrain)")
        train2 ()
    } // train

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Train the classifier by computing the probabilities for Y, and the
     *  conditional probabilities for X_j.
     */
    private def train2 ()
    {
        p_Y = f_Y.toDouble / md                               // "so called" prior probability for class yi
        for (i <- 0 until k; j <- 0 until n) {                // for each class yi & feature xj
            val me_vc = me / vc(j).toDouble
            for (xj <- 0 until vc(j)) {                       // for each value for feature j: xj
                p_X_Y(i, j, xj) = (f_YX(i, j, xj) + me_vc) / (f_Y(i) + me)
            } // for
        } // for

        if (DEBUG) {
            println ("p_Y   = " + p_Y)                        // P(Y = yi)
            println ("p_X_Y = " + p_X_Y)                      // P(X_j = x | Y = yi)
        } // if
    } // train2

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a discrete data vector 'z', classify it returning the class number
     *  (0, ..., k-1) with the highest relative posterior probability.
     *  Return the best class, its name and its relative probability.
     *  This version multiplies probabilities.
     *  @param z  the data vector to classify
     */
    def classify (z: VectoI): (Int, String, Double) =
    {
//      if (DEBUG) banner ("classify (z)")
        val prob = new VectorD (p_Y)
        for (i <- 0 until k) {
            for (j <- 0 until n) prob(i) *= p_X_Y(i, j, z(j))   // P(X_j = z_j | Y = i)
        } // for
        if (DEBUG) println ("prob = " + prob)
        val best = prob.argmax ()                    // class with the highest relative posterior probability
        (best, cn(best), prob(best))                 // return the best class and its name
    } // classify

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given a discrete data vector 'z', classify it returning the class number
     *  (0, ..., k-1) with the highest relative posterior probability.
     *  Return the best class, its name and its relative log-probability.
     *  This version adds log-probabilities.
     *  @param z  the data vector to classify
     */
    def lclassify (z: VectoI): (Int, String, Double) =
    {
//      if (DEBUG) banner ("lclassify (z)")
        val prob = vlog (new VectorD (p_Y))
        for (i <- 0 until k) {
            for (j <- 0 until n) prob(i) += log (p_X_Y(i, j, z(j)))   // P(X_j = z_j | Y = i)
        } // for
        if (DEBUG) println ("prob = " + prob)
        val best = prob.argmax ()                    // class with the highest relative posterior probability
        (best, cn(best), prob(best))                 // return the best class and its name
    } // lclassify

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Take the log of the given probability vector.
     *  @param p  the given probability vector 
     */
    private def vlog (p: VectorD): VectorD = p.map (log (_))

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reset or re-initialize all the population and probability vectors and
     *  hypermatrices to 0.
     */
    def reset ()
    {
        f_Y.set (0)
        f_YX.set (0)
        p_Y.set (0)
        p_X_Y.set (0)
    } // reset

} // NaiveBayes class


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** `NaiveBayes` is the companion object for the `NaiveBayes` class.
 */
object NaiveBayes
{
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a `NaiveBayes` object, passing 'x' and 'y' together in one matrix.
     *  @param xy  the data vectors along with their classifications stored as rows of a matrix
     *  @param fn  the names of the features
     *  @param k   the number of classes
     *  @param vc  the value count (number of distinct values) for each feature
     *  @param me  use m-estimates (me == 0 => regular MLE estimates)
     */
    def apply (xy: MatriI, fn: Array [String], k: Int, cn: Array [String],
               vc: VectoI = null, me: Int = me_default) =
    {
        new NaiveBayes (xy(0 until xy.dim1, 0 until xy.dim2 - 1), xy.col(xy.dim2 - 1), fn, k, cn, vc, me)
    } // apply

} // NaiveBayes object

import scalation.linalgebra.MatrixI

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `NaiveBayesTest` object is used to test the `NaiveBayes` class.
 *  Classify whether a car is more likely to be stolen (1) or not (1).
 *  @see www.inf.u-szeged.hu/~ormandi/ai2/06-naiveBayes-example.pdf
 *  > run-main scalation.analytics.classifier.NaiveBayesTest
 */
object NaiveBayesTest extends App
{
    // x0: Color:   Red (1), Yellow (0)
    // x1: Type:    SUV (1), Sports (0)
    // x2: Origin:  Domestic (1), Imported (0)
    // x3: Mpg:     High (1), Low (0)
    // features:                 x0 x1 x2 x3
    val x = new MatrixI ((10, 4), 1, 0, 1, 1,               // data matrix
                                  1, 0, 1, 0,
                                  1, 0, 1, 1,
                                  0, 0, 1, 1,
                                  0, 0, 0, 1,
                                  0, 1, 0, 0,
                                  0, 1, 0, 0,
                                  0, 1, 1, 1,
                                  1, 1, 0, 0,
                                  1, 0, 0, 0)

    val y  = VectorI (1, 0, 1, 0, 1, 0, 1, 0, 0, 1)         // classification vector: 0(No), 1(Yes))
    val fn = Array ("Color", "Type", "Origin", "Mpg")       // feature/variable names
    val cn = Array ("No", "Yes")                            // class names

    println ("x = " + x)
    println ("y = " + y)
    println ("---------------------------------------------------------------")

    val nb = new NaiveBayes (x, y, fn, 2, cn)               // create the classifier            

    // train the classifier ---------------------------------------------------
    nb.train()

    // test sample ------------------------------------------------------------
    val z1 = VectorI (1, 0, 1, 1)                           // existing data vector to classify
    val z2 = VectorI (1, 1, 1, 0)                           // new data vector to classify
    println ("classify (" + z1 + ") = " + nb.classify (z1) + "\n")
    println ("classify (" + z2 + ") = " + nb.classify (z2) + "\n")

//  nb.crossValidateRand ()                                 // cross validate the classifier

} // NaiveBayesTest object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `NaiveBayesTest2` object is used to test the 'NaiveBayes' class.
 *  Given whether a person is Fast and/or Strong, classify them as making Y = 1
 *  or not making Y = 0 the football team.
 *  > run-main scalation.analytics.classifier.NaiveBayesTest2
 */
object NaiveBayesTest2 extends App
{
    // training-set -----------------------------------------------------------
    // x0: Fast
    // x1: Strong
    // y:  Classification (No/0, Yes/1)
    // features:                  x0 x1  y
    val xy = new MatrixI ((10, 3), 1, 1, 1,
                                   1, 1, 1,
                                   1, 0, 1,
                                   1, 0, 1,
                                   1, 0, 0,
                                   0, 1, 0,
                                   0, 1, 0,
                                   0, 1, 1,
                                   0, 0, 0,
                                   0, 0, 0)

    val fn = Array ("Fast", "Strong")                     // feature names
    val cn = Array ("No", "Yes")                          // class names

    println ("xy = " + xy)
    println ("---------------------------------------------------------------")

    val nb = NaiveBayes (xy, fn, 2, cn, null)             // create the classifier

    // train the classifier ---------------------------------------------------
    nb.train ()

    // test sample ------------------------------------------------------------
    val z = VectorI (1, 0)                                // new data vector to classify
    println ("classify (" + z + ") = " + nb.classify (z) + "\n")

    nb.crossValidate ()                                   // cross validate the classifier

} // NaiveBayesTest2 object


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `NaiveBayesTest3` object is used to test the 'NaiveBayes' class.
 *  @see archive.ics.uci.edu/ml/datasets/Lenses
 *  @see docs.roguewave.com/imsl/java/7.3/manual/api/com/imsl/datamining/NaiveBayesClassifierEx2.html
 *  > run-main scalation.analytics.classifier.NaiveBayesTest3
 */
object NaiveBayesTest3 extends App
{
   // training-set -----------------------------------------------------------
    // y:  Classification (1): hard contact lenses, (2) soft contact lenses, (3) no contact lenses
    // x0. Age of the patient: (1) young, (2) pre-presbyopic, (3) presbyopic
    // x1. Spectacle prescription:  (1) myope, (2) hypermetrope
    // x2. Astigmatic:     (1) no, (2) yes
    // x3. Tear production rate:  (1) reduced, (2) normal
    // features:                  x0  x1  x2  x3   y
    val xy = new MatrixI ((24, 5), 1,  1,  1,  1,  3,           // 1
                                   1,  1,  1,  2,  2,           // 2
                                   1,  1,  2,  1,  3,           // 3
                                   1,  1,  2,  2,  1,           // 4
                                   1,  2,  1,  1,  3,           // 5
                                   1,  2,  1,  2,  2,           // 6
                                   1,  2,  2,  1,  3,           // 7
                                   1,  2,  2,  2,  1,           // 8
                                   2,  1,  1,  1,  3,           // 9
                                   2,  1,  1,  2,  2,           // 10
                                   2,  1,  2,  1,  3,           // 11
                                   2,  1,  2,  2,  1,           // 12
                                   2,  2,  1,  1,  3,           // 13
                                   2,  2,  1,  2,  2,           // 14
                                   2,  2,  2,  1,  3,           // 15
                                   2,  2,  2,  2,  3,           // 16
                                   3,  1,  1,  1,  3,           // 17
                                   3,  1,  1,  2,  3,           // 18
                                   3,  1,  2,  1,  3,           // 19
                                   3,  1,  2,  2,  1,           // 20
                                   3,  2,  1,  1,  3,           // 21
                                   3,  2,  1,  2,  2,           // 22
                                   3,  2,  2,  1,  3,           // 23
                                   3,  2,  2,  2,  3)           // 24

    xy -= 1

    val fn = Array ("Age", "Spectacle", "Astigmatic", "Tear")   // feature names
    val cn = Array ("Hard", "Soft", "Neither")                       // class names

    println ("xy = " + xy)
    println ("---------------------------------------------------------------")

    val nb = NaiveBayes (xy, fn, 3, cn, null, 0)                // create the classifier

    // train the classifier ---------------------------------------------------
    nb.train ()

    // test all rows ------------------------------------------------------------
    for (i <- xy.range1) {
        val z  = xy(i).slice (0, 4)                             // x-values
        val y  = xy(i, 4)                                       // y-value
        val yp = nb.classify (z)                                // y predicted
        println (s"yp = classify ($z) = $yp,\t y = $y,\t ${cn(y)}")
    } // for
                                  
} // NaiveBayesTest3 object

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `NaiveBayesTest4` object is used to test the 'NaiveBayes' class.
 *  > run-main scalation.analytics.classifier.NaiveBayesTest4
 */
object NaiveBayesTest4 extends App
{
    val filename = BASE_DIR + "breast-cancer.arff"
    var data = Relation (filename, -1, null)
    val xy = data.toMatriI2 (null)
    val fn = data.colName.toArray
    val cn = Array ("0", "1")                              // class names
    val nb = NaiveBayes (xy, fn, 2, cn, null, 0)           // create the classifier
    nb.train ()
    nb.crossValidate ()

} // NaiveBayesTest4 object
