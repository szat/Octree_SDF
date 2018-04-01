# OctreeSDF
An octree for ray / mesh intersection queries in pure C++. There is still some work to do. The main idea to use an octree where each leaf stores a triangle, represented by its barycenter. 

# Installation
To make this work, you only need to drag and drop "octree.h" into your project. Note that the code is expecting your data to be in the "vector<array<T,3>>" format, where each array represents a point or vector. Faces are also represented like that. If you are using EIGEN, there is a small templated function in main.cpp to convert MatrixXd into vector<array<double,3>>. 

To use libigl for visualisation, the installation of libigl is kinda annoying. Here is a small tutorial. 

# TODO:
-Handling of very large meshes
-Automatic visualisation
-Multithreading
