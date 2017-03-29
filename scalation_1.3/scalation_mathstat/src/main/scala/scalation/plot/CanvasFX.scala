

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael E. Cotterell
  *  @edited Alec Molinaro
  *  @version 1.2
  *  @date    Sat Feb 21 00:05:16 EST 2015
  *  @see     LICENSE (MIT style license file).
  */

package scalation.plot

import javafx.application.Platform
import javafx.embed.swing.JFXPanel
import javafx.scene.{Group, Node, Scene}


import javax.swing.JFrame
import javax.swing.JPanel
import javax.swing.JScrollPane
import java.awt.Dimension
import javax.swing.SwingUtilities
import javax.swing.WindowConstants._

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** A simple canvas for displaying JavaFX nodes.
  * @param title   canvas title
  * @param width   canvas width
  * @param height  canvas height
  */
class CanvasFX (val title: String, val width: Int, val height: Int) extends Runnable {

  private var scene: Scene = null

  /** Underlying `JFrame` for this `Canvas`. */
  val frame   = new JFrame (title)

    val sp = new JScrollPane ()

  /** Underlying `JFXPanel` for this `Canvas`. */
  val fxPanel = new JFXPanel ()

  /** Underlying `Group` that is the root of the scene
    * graph for this `Canvas`.
    */
  val nodes   = new Group ()       //

  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /** Launches the `Canvas` contents inside of a JavaFX Application Thread.
    */
  protected def run (): Unit = {
    Platform.runLater (new Runnable() {
      def run(): Unit =
      {

        val minimumSize = new Dimension (width, height)
        fxPanel.setMinimumSize (minimumSize)


        scene = new Scene (nodes, width, height)
        fxPanel.setScene (scene)

        frame.add (new JScrollPane(fxPanel))
        frame.setSize (width, height)
        frame.setVisible (true)
        frame.setDefaultCloseOperation (DISPOSE_ON_CLOSE)

        nodes.setAutoSizeChildren(true);

      } // run
    }) // runLater
  } // run

  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /** Display the canvas.
    */
  def show (): Unit = SwingUtilities.invokeLater (this)

  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /** Adds a node to the canvas's scene graph.
    * @param node  a `Node` to add
    */
  def add (node: Node): Unit = nodes.getChildren.add (node)

} // CanvasFX

object CanvasTest extends App {

  import javafx.scene.shape.Circle

  val canvas = new CanvasFX ("Test0", 640, 480)
  canvas.add (new Circle (40, 40, 30))
  canvas.show ()

  val canvas1 = new CanvasFX ("Test1", 640, 480)
  canvas1.add (new Circle (80, 80, 30))
  canvas1.show ()

} // App
