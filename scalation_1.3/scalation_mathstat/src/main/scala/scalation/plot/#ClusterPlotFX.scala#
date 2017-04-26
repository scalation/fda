//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael E. Cotterell & Alec Molinaro
  *  @version 1.2
  *  @date    THR Apr 20 00:05:16 EST 2015
  *  @see     LICENSE (MIT style license file).
  */


package scalation.plot

import javafx.application.Application
import javafx.event.{ActionEvent, EventHandler}
import javafx.geometry.Insets
import javafx.scene.{Group, Node, Scene}
import javafx.scene.chart._
import javafx.scene.paint.Color
import javafx.stage.Stage
import javafx.scene.control.ToggleButton
import javafx.scene.layout.{HBox, VBox}

import scala.math._
import scalation.linalgebra.{MatriD, VectoD}
import scalation.linalgebra.{VectoD, VectoI}
import scalation.linalgebra.{MatrixD, VectorD}
import scalation.scala2d.{Panel, VizFrame}
import scala.math.pow


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `ClusterPlotFX` class takes an 'x' vector and a 'y' matrix of data values and plots
  *  the '(x, y_i)' data points for each row 'y_i' of the matrix.
  *  @param x       the x vector of data values (horizontal)
  *  @param y       the y matrix of data values where y(i) is the i-th vector (vertical)
  *  @param _title  the title of the plot
  */
class ClusterPlotFX (x: MatriD, cl: Array [Int], k: Int, _title: String = "ClusterPlotFX")
  extends VizFrame (_title, null) {

  // creating canvas
  val c = new CanvasFX(_title, 900, 800)

  // defining the axes
  val xAxis: NumberAxis = new NumberAxis()
  val yAxis: NumberAxis = new NumberAxis()

  // creating the chart
  var scatterChart: ScatterChart[Number, Number] =
  new ScatterChart[Number, Number](xAxis, yAxis)

  // creating labels
  xAxis.setLabel("X Axis")
  yAxis.setLabel("Y Axis")
  scatterChart.setTitle("Data Clustering")

  // setting format and dimensions
  scatterChart.setPrefSize(900, 700);
  val scene: Scene = new Scene(new Group())
  val vbox: VBox = new VBox()
  val hbox: HBox = new HBox()
  hbox.setSpacing(10)

  val series = Array.ofDim[XYChart.Series[Number, Number]] (k)

  for (i <- 0 until k ) {

    series(i) = new XYChart.Series[Number, Number]()
    series(i).setName(s"Cluster $i")

  } // for

  // establishing clusters
  for(i <- 0 until x.dim1) {

    val point = new XYChart.Data[Number, Number] (x(i, 0), x(i, 1))
    series(cl(i)).getData.add(point)

  } // for

  for (i <- 0 until k) scatterChart.getData.add(series(i)) // for

  // putting it all together
  vbox.getChildren.addAll(scatterChart, hbox)
  hbox.setPadding(new Insets(10, 10, 10, 50))
  c.nodes.getChildren.addAll(vbox)
  c.show()

  println(x)

} // ClusterPlotFX


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `ClusterPlotFX` companion object provides a builder method for plotting several
  *  'y' vectors versus an 'x' vector.
  */
object ClusterPlotFX
{

  import scalation.linalgebra.{MatrixD, VectorD}

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /** Create a plot of several 'y' vectors versus an 'x' vector.
    *  @param x  the x vector of data values (horizontal)
    *  @param y  one or more vectors of values where y(i) is the i-th vector (vertical)
    */
  def apply (x: VectoD, y: VectorD*)
  {
    val yy = new MatrixD (y.length, x.dim)
    for (i <- 0 until y.length) yy(i) = y(i)
    new PlotMFX (x, yy)
  } // apply

} // ClusterPlotFX object


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `ClusterPlotFXTest` object is used to test the `ClusterPlotFX` class.
  */
object ClusterPlotFXTest extends App {

  import scalation.linalgebra.{MatrixD, VectorD}
  import scalation.random.{Normal, Bernoulli}

  val n     = 100
  val coin  = Bernoulli ()
  val dist1 = Normal (2.0, 1.0)
  val dist2 = Normal (8.0, 1.0)
  val v     = new MatrixD (n, 2)    // data points
  val cls   = Array.ofDim [Int] (n) // basis for cluster assignment
  val k     = 4                     // number of clusters

  for (i <- v.range1) {
    val c1 = coin.gen.toInt
    val c2 = coin.gen.toInt
    val d1 = dist1.gen
    val d2 = dist2.gen
    println (s"c1 = $c1; c2 = $c2; d1 = $d1; d2 = $d2")

    (c1, c2) match {
      case _ if (c1 == 0 && c2 == 0) => {cls(i) = 0; v(i, 0) = dist1.gen; v(i, 1) = dist1.gen; println("1")}
      case _ if (c1 == 0 && c2 == 1) => {cls(i) = 1; v(i, 0) = dist1.gen; v(i, 1) = dist2.gen;println("2")}
      case _ if (c1 == 1 && c2 == 0) => {cls(i) = 2; v(i, 0) = dist2.gen; v(i, 1) = dist1.gen;println("3")}
      case _ if (c1 == 1 && c2 == 1) => {cls(i) = 3; v(i, 0) = dist2.gen; v(i, 1) = dist2.gen;println("4")}
      //case _ => {}

    } // match
  } // for

  println(s"!!!!!!!!!!!!!!!!!!!!!!!!!!! ${dist1.gen}")

  println(v)
  val plot = new ClusterPlotFX (v, cls, k)
    println ("plot = " + plot)

} // ClusterPlotFXTest object