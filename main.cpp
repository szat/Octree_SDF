///only move sematic, not redundant
//https://docs.microsoft.com/en-gb/cpp/cpp/explicitly-defaulted-and-deleted-functions
//http://blog.matejzavrsnik.com/using_smart_pointers_in_building_trees_in_which_child_nodes_need_an_access_to_the_parent.html
//http://www.opengl-tutorial.org/miscellaneous/clicking-on-objects/picking-with-custom-ray-obb-function/

//#include "stdafx.h"
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <igl/median.h>
#include <igl/igl_inline.h>
#include <igl/matrix_to_list.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>

#include <igl/doublearea.h>
#include <igl/internal_angles.h>
#include <igl/is_irregular_vertex.h>

#include "tutorial_shared_path.h"
#include <iostream>
//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
#include <float.h>
#include <string>
#include <array>

#include "octree.h"
#include "testing.h"

Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
	//test_RayBox();
	//test_RayTriangle();

	vector<array<int, 2>> vec = { {1,2},{3,4},{5,4},{7,8} };
	accumulate(vec.begin(), vec.end(), 0);
	cin.ignore();
	////Read data
	//MatrixXd V;
	//MatrixXi F;
	//MatrixXd N;
	//igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", V, F);
	//igl::per_face_normals(V, F, N);

	//MatrixXd C(F.rows(), 3);

	//// Plot the mesh
	//igl::opengl::glfw::Viewer viewer;
	//viewer.data().set_mesh(V, F);

	//// Launch the viewer
	//viewer.launch();
}