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

//#include "octree.h"

Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

bool dud1(const array<double,3>& in, array<double,3>& out) {
	for (size_t  i = 0; i < in.size(); ++i) {
		out[i] = 2*in[i];
	}
	return true;
}

double dot(const array<double, 3> &A, const array<double, 3> &B) {
	double out = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
	return out;
}

array<double,3> cross(const array<double,3> &A, const array<double,3> &B) {
	//a = A.y * B.z - A.z * B.y;
	//b = A.z * B.x - A.x * B.z;
	//c = A.x * B.y - A.y * B.x;
	array<double,3> out;
	out[0] = A[1] * B[2] - A[2] * B[1];
	out[1] = A[2] * B[0] - A[0] * B[2];
	out[2] = A[0] * B[1] - A[1] * B[0];
	return out;
}

int main(int argc, char *argv[])
{
	
	array<double, 3> A = { 3,4,5 };
	array<double, 3> B = { 4,3,5 };
	array<double, 3> C = cross(A, B);

	array<array<double, 3>, 2>  box;
	box[0] = { 1,2,3 };
	box[1] = { 7,8,9 };
	cout << "box " << box[0][0] << " " << box[0][1] << " " << box[0][2] << endl;
	
	for (auto x : C) cout << x << " ";

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