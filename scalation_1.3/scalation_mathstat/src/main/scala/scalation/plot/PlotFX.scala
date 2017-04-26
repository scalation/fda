/**
  * Created by alec on 2/6/17.
  */

package scalation.plot

import javafx.application.Application
import javafx.scene.Scene
import javafx.scene.chart.NumberAxis
import javafx.scene.chart.LineChart
import javafx.scene.chart.XYChart
import javafx.scene.chart.Axis
import javafx.scene.paint.Color
import javafx.stage.Stage
import javafx.scene.Node
import javafx.scene.control.ToggleButton
import javafx.scene.Group;
import javafx.scene.layout.HBox
import javafx.scene.layout.VBox
import javafx.event.ActionEvent
import javafx.event.EventHandler
import javafx.geometry.Insets



import scalation.linalgebra.{VectoD, VectoI}
import scalation.scala2d.{Panel, VizFrame}
import scala.math.pow


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `PlotFX` class takes 'x' and 'y' vectors of data values and plots the '(x, y)'
  *  data points. For more vertical vectors use `PlotM`.
  *  @param x       the x vector of data values (horizontal)
  *  @param y       the y vector of data values (primary vertical)
  *  @param z       the z vector of data values (secondary vertical) to compare with y
  *  @param _title  the title of the plot
  */
class PlotFX (x: VectoD, y: VectoD, z: VectoD, _title: String = "Plot y vs. x")
  extends VizFrame (_title, null)
{
 // getContentPane ().add (new Canvas (x, y, z, getW, getH))
 // setVisible (true)

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

  // defining the series
  var series1: XYChart.Series[Number, Number] =
    new XYChart.Series[Number, Number]()

  //if (z != null) {
    var series2: XYChart.Series[Number, Number] =
      new XYChart.Series[Number, Number]()
  //} // if


  series1.setName("Quadratic")
  series2.setName("Cubic")

  for (i <- 0 until x.dim) {
    val xy = new XYChart.Data[Number, Number](x(i), y(i))
    series1.getData.add(xy)

  } // for

  for (i <- 0 until x.dim) {
    val xz = new XYChart.Data[Number, Number](x(i), z(i))
    series2.getData.add(xz)

  } // for

  lineChart.getData.add(series1)
  lineChart.getData.add(series2)

  lineChart.setPrefSize(800, 600);
  val scene: Scene = new Scene(new Group())
  val vbox: VBox = new VBox()
  val hbox: HBox = new HBox()

  val btn1: ToggleButton = new ToggleButton()
  btn1.setText("Hide " + series1.getName())
  btn1.setOnAction(new EventHandler[ActionEvent]() {
    override def handle(event: ActionEvent): Unit = {
      val source: ToggleButton =
        event.getSource.asInstanceOf[ToggleButton]
      if (source.isSelected()) {
        setT(0, lineChart)
        btn1.setText("Show " + series1.getName())
      } else {
        setC(0, lineChart)
        btn1.setText("Hide " + series1.getName())
      }
    }
  })

  val btn2: ToggleButton = new ToggleButton()
  btn2.setText("Hide " + series2.getName())
  btn2.setOnAction(new EventHandler[ActionEvent]() {
    override def handle(event: ActionEvent): Unit = {
      val source: ToggleButton =
        event.getSource.asInstanceOf[ToggleButton]
      if (source.isSelected()) {
        setT(1, lineChart)
        btn2.setText("Show " + series2.getName())
      } else {
        setC(1, lineChart)
        btn2.setText("Hide " + series2.getName())
      }
    }
  })


  hbox.setSpacing(10)
  hbox.getChildren.addAll(btn1)
  hbox.getChildren.addAll(btn2)
  vbox.getChildren.addAll(lineChart, hbox)
  hbox.setPadding(new Insets(10, 10, 10, 50))


  c.nodes.getChildren.addAll(vbox)
  c.show()


  def setT(n: Int, x: LineChart[Number, Number]): Unit = {

    x.lookup(".default-color" + n + ".chart-series-line").setStyle("-fx-stroke: transparent;")

  } //setT

  def setC(n: Int, x: LineChart[Number, Number]): Unit = {

    x.lookup(".default-color" + n + ".chart-series-line").setStyle("-fx-stroke: DEFAULT_COLOR+(nextClearBit%8);")

  } //setT


} // PlotFX class


object PlotFXTest extends App {
  import scalation.linalgebra.VectorD
  /*
      val x = VectorD (0.0, 1.0, 2.0, 3.0,  4.0,  5.0,  6.0, 7.0, 8.0, 9.0, 10.0)
      val y = VectorD (0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 16.0, 9.0, 4.0, 1.0,  0.0)
   */
    val n = 40
  val x = new VectorD (n)
  val y = new VectorD (n)
  val z = new VectorD (n)
  for (i <- 0 until n) { x(i) = i / 10.0; y(i) = pow (x(i) - 5, 2) }
  for (i <- 0 until n) { x(i) = i / 10.0; z(i) = pow (x(i) - 5, 3) }
  val plot = new PlotFX (x, y, z)
  println ("plot = " + plot)

} // PlotFXTest
