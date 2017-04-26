//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael E. Cotterell & Alec Molinaro
  *  @version 1.2
  *  @date    THR Apr 20 00:05:16 EST 2015
  *  @see     LICENSE (MIT style license file).
  */


package scalation.plot

import javafx.application.Application
import javafx.scene.Scene
import javafx.scene.chart.NumberAxis
import javafx.scene.chart.LineChart
import javafx.scene.chart.XYChart
import javafx.scene.chart.Axis
import javafx.scene.shape.{Line, Rectangle}
import javafx.scene.paint.Color

import java.io.File

import scalation.linalgebra.{VectoD, VectoI}
import scalation.scala2d.{Panel, VizFrame}
import scala.math.pow

import javafx.application.Application
import javafx.event.{ActionEvent, EventHandler}
import javafx.geometry.Insets
import javafx.scene.{Group, Node, Scene}
import javafx.scene.chart._
import javafx.scene.paint.Color
import javafx.stage.Stage
import javafx.scene.control.{ToggleButton, ScrollPane}

import javafx.scene.control.ScrollPane.ScrollBarPolicy
import javafx.scene.layout.{Pane, HBox, VBox}

import scala.math._
import scalation.linalgebra.{MatriD, VectoD}
import scalation.linalgebra.{VectoD, VectoI}
import scalation.linalgebra.{MatrixD, VectorD}
import scalation.scala2d.{Panel, VizFrame}
import scala.math.pow


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `HeatMapFX` class takes 'x' and 'y' vectors of data values and plots the '(x, y)'
  *  data points. For more vertical vectors use `PlotM`.
  *  @param x       the x matrix of data values (horizontal)
  *  @param cl      the cl array of data values (primary vertical)
  *  @param k       the k value of data values (secondary vertical) to compare with y
  *  @param _title  the title of the plot
  */
class HeatMapFX (x: MatriD, cl: Array [Int], k: Int, _title: String = "HeatMapFX")
  extends VizFrame (_title, null) {

    val WIDTH = 200
    val HEIGHT = 10
    val WAVE_LOW = 380.0
    val WAVE_HIGH = 780.0

    val min = x.min()
    val max = x.max()

    // creating canvas
    val c = new CanvasFX(_title, 700, 600)

    // creating the scene
    val scene: Scene = new Scene(new Group())
    val vbox: Pane = new Pane()
    val hbox: HBox = new HBox()
    hbox.setSpacing(10)

    // creating a scroll pane
    val sp = new ScrollPane()
    sp.setPrefSize(700, 578)


    var map = Array.ofDim[Rectangle](x.dim1, x.dim2)

    for (i <- 0 until x.dim1; j <- 0 until x.dim2) {

      val r = new Rectangle()
      r.setWidth(WIDTH)
      r.setHeight(HEIGHT)
      r.setX(WIDTH*j+j)
      r.setY(HEIGHT*i+i)
      val ci = wToRGB(linearInterp(x(i, j)))
      r.setFill(ci)

      map(i)(j) = r

      vbox.getChildren.addAll(r)

    } // for

    // putting it all together
    vbox.getChildren.addAll(hbox)
    hbox.setPadding(new Insets(10, 10, 10, 50))

    // sizing the scroll bar
    sp.setHbarPolicy(ScrollBarPolicy.AS_NEEDED)
    sp.setVbarPolicy(ScrollBarPolicy.AS_NEEDED)
    sp.setPannable(true)
    sp.setFitToHeight(true)

    // packing it up and sending it off
    sp.setContent(vbox)
    c.nodes.getChildren.addAll(sp)
    c.show()

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /** Linearally interpolates a plot of several 'v' vectors versus an 'x' vector.
    *  @param v  the v vector of data values (horizontal)
    */
    def linearInterp(v: Double): Double = {
      val wavelength = ( (v-min)/(max-min) )  * (WAVE_HIGH - WAVE_LOW) + WAVE_LOW
      return wavelength
    } // linearInterp


  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /** Translates a plot of a 'w' vector versus an 'x' vector.
    *  @param w  the w vector of data values (horizontal)
    *  @see http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm
    */
    def wToRGB (w: Double): Color = {

      val GAMMA = 0.80
      val INTENSITY_MAX = 255.0
      val INTENSITY_MIN = 0.0

      var r = 0.0
      var b = 0.0
      var g = 0.0
      var factor = 0.0


      def adjust (color: Double, factor: Double): Int = {
        val result = if (color == 0) 0
                     else math.round(INTENSITY_MAX * math.pow(color * factor, GAMMA))

        result.toInt
      }   //adjust


      if (380 <= w && w < 440 ) {
        r = -(w - 440) / (440 - 380)
        g = 0
        b = 1
      } else if (440 <= w && w < 490){
        r = 0
        g = (w - 440) / (490 - 440)
        b = 1
      } else if (490 <= w && w < 510) {
        r = 0
        g = 1
        b = - (w - 510) / (510 - 490)
      } else if (510 <= w && w < 580) {
        r = (w - 510) / (580 - 510)
        g = 1
        b = 0
      } else if (580 <= w && w < 645) {
        r = 1
        g = - (w - 645) / (645 - 580)
        b = 0
      } else if (645 <= w && w <= 780) {
        r = 1
        g = 0
        b = 0
      } else  {
        r = 0
        g = 0
        b = 0
      } // if


      if      (380 <= w && w <  420) factor = 0.3 + 0.7 * (w - 380) / (420 - 380)
      else if (420 <= w && w <  701) factor = 1
      else if (701 <= w && w <= 780) factor = 0.3 + 0.7 * (780 - w) / (780 - 700)
      else                           factor = 0

      r = adjust(r, factor)
      g = adjust(g, factor)
      b = adjust(b, factor)

      Color.rgb(r.toInt, g.toInt, b.toInt)

    } // wToXYZ
} // HeatMapFX class


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `HeatMapFXTest` object is used to test the `HeatMapFXTest` class.
  */
  object HeatMapFXTest extends App {

    val f = new File("nzdata.csv")
    println(f.getAbsolutePath())

    val m = MatrixD (f.getAbsolutePath())

    val hmap = new HeatMapFX(m, null, 0, null)

  } // HeatMapFXTest


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `HeatMapFXTest2` object is used to test the `HeatMapFXTest2` class.
  */
  object HeatMapFXTest2 extends App {

    val f = new File("random.csv")
    println(f.getAbsolutePath())

    val m = MatrixD (f.getAbsolutePath())

    val hmap = new HeatMapFX(m, null, 0, null)

  } // HeatMapFXTest2


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `HeatMapFXTest3` object is used to test the `HeatMapFXTest3` class.
  */
  object HeatMapFXTest3 extends App {

    val f = new File("randomSorted.csv")
    println(f.getAbsolutePath())

    val m = MatrixD (f.getAbsolutePath())

    val hmap = new HeatMapFX(m, null, 0, null)

  } // HeatMapFXTest3
