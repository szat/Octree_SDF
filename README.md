# OctreeSDF
An octree for ray / mesh intersection queries in pure C++STL. There is still some work to do. The main idea to use an octree where each leaf stores a triangle, represented by its barycenter. Each node, including the leaves, store a box that bounds all the vertices of all the triangles under it (thus the barycenters will certainly be inside). 

Youtube link: https://www.youtube.com/watch?v=4NMcUA6rHvQ

# Installation
To make this work, you only need to drag and drop "octree.h" into your project. Note that the code is expecting your data to be in the "vector<array<T,3>>" format, where each array represents a point or vector. Faces are also represented by "vector<array<int,3>>", storing the indices of the vertices which are the points of the triangles. 

If you are using EIGEN, there is a small templated function in main.cpp to convert "MatrixXd" into "vector<array<double,3>>". 

# libigl 
The installation of libigl, which is used only for visualisation, is kinda annoying. Here is a small tutorial. TODO.

# TODO:
-libigl installation tutorial.
-Handling of very large meshes
-Automatic visualisation
-Multithreading
