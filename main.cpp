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
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <float.h>
#include <string>
#include <array>

#include "octree.h"
#include "testing.h"

Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

template <class T>
vector<array<T, 3>> Eigen2CPP(Matrix<T, Dynamic, 3> & mat) {
	T* T_ptr = mat.data();
	vector<T> col0(T_ptr, T_ptr + mat.rows());
	vector<T> col1(T_ptr + mat.rows(), T_ptr + 2 * mat.rows());
	vector<T> col2(T_ptr + 2 * mat.rows(), T_ptr + 3 * mat.rows());
	vector<array<T, 3>> out;
	for (size_t i = 0; i < mat.rows(); ++i) {
		out.push_back({ col0.at(i), col1.at(i), col2.at(i) });
	}
	return out;
}

void test_Eigen2CPP() {
	Matrix<int, Dynamic, 3> F = (MatrixXi(4,3) << 1, 2, 3, 
		1, 2, 3, 
		5, 6, 7, 
		7, 8, 9).finished();
	vector<array<int, 3>> matF = Eigen2CPP(F);
	for (int row = 0; row < F.rows(); ++row) {
		for (int col = 0; col < 3; ++col) {
			cout << matF.at(row)[col] << " ";
		}
		cout << endl;
	}
	cin.ignore();
}

int main(int argc, char *argv[])
{
	//test_RayBox();
	//test_RayTriangle();
	//test_Eigen2CPP();
	//test_build();


	//Read data
	Matrix<double, Dynamic,3> V;
	Matrix<int, Dynamic, 3> F;
	Matrix<double, Dynamic, 3> N;
	Matrix<double, Dynamic, 3> B;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", V, F);
	igl::per_face_normals(V, F, N);
	igl::barycenter(V, F, B);
	
	cout << "V size = " << V.rows() << endl;
	cout << "F size = " << F.rows() << endl;
	cout << "BC size " << B.rows() << endl;

	vector<array<double,3>> vecV = Eigen2CPP(V);
	vector<array<int, 3>> vecF = Eigen2CPP(F);
	vector<array<double, 3>> vecN = Eigen2CPP(N);
	vector<array<double, 3>> vecB = Eigen2CPP(B);


	cout << "vecV size = " << vecV.size() << endl;
	cout << "vecF size = " << vecF.size() << endl;
	cout << "vecB size " <<vecB.size() << endl;

	SDF tree(vecV, vecF, vecB);
	tree.init();
	tree.build();
	tree.test();


	cout << "almost" << endl;
	cin.ignore();

	//MatrixXd C(F.rows(), 3);

	//// Plot the mesh
	//igl::opengl::glfw::Viewer viewer;
	//viewer.data().set_mesh(V, F);
	//viewer.data().add_points(BC, Eigen::RowVector3d(0, 0, 0));
	//// Launch the viewer
	//viewer.launch();
}