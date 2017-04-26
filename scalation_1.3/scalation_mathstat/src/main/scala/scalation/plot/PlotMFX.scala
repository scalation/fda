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
import javafx.scene.chart.NumberAxis
import javafx.scene.chart.LineChart
import javafx.scene.chart.XYChart
import javafx.scene.chart.Axis
import javafx.scene.paint.Color
import javafx.stage.Stage
import javafx.scene.control.ToggleButton
import javafx.scene.layout.{HBox, VBox}

import scala.math.{ceil, floor, min, pow, round}
import scalation.linalgebra.{MatriD, VectoD}
import scalation.linalgebra.{VectoD, VectoI}
import scalation.linalgebra.{MatrixD, VectorD}
import scalation.scala2d.{Panel, VizFrame}
import scala.math.pow


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `PlotM` class takes an 'x' vector and a 'y' matrix of data values and plots
  *  the '(x, y_i)' data points for each row 'y_i' of the matrix.
  *  @param x       the x vector of data values (horizontal)
  *  @param y       the y matrix of data values where y(i) is the i-th vector (vertical)
  *  @param _title  the title of the plot
  */
class PlotMFX (x: VectoD, y: MatriD, _title: String = "PlotM y_i vs. x for each i")
  extends VizFrame (_title, null)
{

  // creating canvas
  val c = new CanvasFX(_title, 800, 700)

  // defining the axes
  val xAxis: NumberAxis = new NumberAxis()
  val yAxis: NumberAxis = new NumberAxis()

  // creating the chart
  var lineChart: LineChart[Number, Number] =
  new LineChart[Number, Number](xAxis, yAxis)

  // creating labels
  xAxis.setLabel("X Axis")
  yAxis.setLabel("Y Axis")
  lineChart.setTitle("Data Overview")

  // setting format and dimensions
  lineChart.setPrefSize(800, 600);
  val scene: Scene = new Scene(new Group())
  val vbox: VBox = new VBox()
  val hbox: HBox = new HBox()
  hbox.setSpacing(10)


  // establishing information and controls
  for(i <- 0 until y.dim1) {

    var series: XYChart.Series[Number, Number]=
	    new XYChart.Series[Number, Number]()

    series.setName(s"row $i")

	  for (j <- 0 until x.dim) {
	    val xy = new XYChart.Data[Number, Number] (x(j), y(i,j))
	    series.getData.add(xy)
	  } // for

  	lineChart.getData.add(series)

    val btn: ToggleButton = new ToggleButton()
    btn.setText("Hide " + series.getName())
    btn.setOnAction(new EventHandler[ActionEvent]() {
      override def handle(event: ActionEvent): Unit = {
        val source: ToggleButton =
          event.getSource.asInstanceOf[ToggleButton]
        if (source.isSelected()) {
          setT(i, lineChart)
          btn.setText("Show " + series.getName())
        } else {
          setC(i, lineChart)
          btn.setText("Hide " + series.getName())
        } // if
      } // handle
    }) // btn

    hbox.getChildren.addAll(btn)

  } // for


  // putting it all together
  vbox.getChildren.addAll(lineChart, hbox)
  hbox.setPadding(new Insets(10, 10, 10, 50))
  c.nodes.getChildren.addAll(vbox)
  c.show()

  // setting data transparent
  def setT(n: Int, x: LineChart[Number, Number]): Unit = {

    x.lookup(".default-color" + n + ".chart-series-line").setStyle("-fx-stroke: transparent;")

  } //setT

  // returning data to original color
  def setC(n: Int, x: LineChart[Number, Number]): Unit = {

    x.lookup(".default-color" + n + ".chart-series-line").setStyle("-fx-stroke: DEFAULT_COLOR+(nextClearBit%8);")

  } //setT

} // PlotMFX class


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `PlotM` companion object provides a builder method for plotting several
  *  'y' vectors versus an 'x' vector.
  */
object PlotMFX
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

} // PlotMFX object

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `PlotMTest` object is used to test the `PlotM` class.
  */
object PlotMFXTest extends App
{
  import scalation.linalgebra.{MatrixD, VectorD}

  val x = new VectorD (100)
  val y = new MatrixD (5, 100)

  for (i <- 0 until 100) {
    x(i)    = (i - 100) / 10.0
    y(0, i) = 10.0 * x(i)
    y(1, i) = pow (x(i), 2)
    y(2, i) = .1 * pow (x(i), 3)
    y(3, i) = .01 * pow (x(i), 4)
    y(4, i) = .001 * pow (x(i), 5)
  } // for
val plot = new PlotMFX (x, y)
  println ("plot = " + plot)

} // PlotMTest object

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `PlotMTest2` object is used to test the `PlotM` class.
  */
object PlotMFXTest2 extends App
{
  import scalation.linalgebra.{MatrixD, VectorD}

  val x = new VectorD (10)
  val y = new MatrixD (2, 10)

  for (i <- 0 until 10) {
    x(i)    = i
    y(0, i) = 0
    y(1, i) = math.sin(i)  

  } // for
val plot = new PlotMFX (x, y)
  println ("plot = " + plot)

} // PlotMTest2 object

