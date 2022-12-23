# 2D-sketch-to-3D-mesh-generation-in-real-time
A sketching interface that transforms 2D freeform strokes drawn on the artboard to a 3D model mesh in OpenGL. According to the research Teddy, all the processes are implemented (Igarashi et al., 1999).

**Dependencies required** 
1. [OpenGL](https://www.opengl.org/)
2. [ImGui](https://github.com/ocornut/imgui)
3. [Poly2Tri](https://github.com/greenm01/poly2tri)
4. [glm](https://glm.g-truc.net/)
5. [GLFW](https://www.glfw.org/)

We have Already Included Poly2Tri and ImGui

**How to run this code**
1. Clone the git repository or download the zip folder.
2. Make sure that you possess the OpenGL library installed
3. Build the source code using the command: make
4. Use the command ./2Dto3D to start sketching.
 
**How the code works**
1. Drag your mouse while pressing the left mouse button to sketch; as soon as you release the button, your shape will be rescaled to clear extra strokes and make a clean shape. And then, the CDT  ( constrained Delaunay triangulation) will happen, and you will be able to see the triangular mesh.
2. Press the ctrl key to apply to prune and implement chordial axis.
3. Click the right mouse button to sew triangles and then elevate the edges proportional to the distance between the vertex and the external vertice.
4. Press Alt key to rotate the 3D mesh model.

