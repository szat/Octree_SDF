///only move sematic, not redundant
//https://docs.microsoft.com/en-gb/cpp/cpp/explicitly-defaulted-and-deleted-functions
//http://blog.matejzavrsnik.com/using_smart_pointers_in_building_trees_in_which_child_nodes_need_an_access_to_the_parent.html
//http://www.opengl-tutorial.org/miscellaneous/clicking-on-objects/picking-with-custom-ray-obb-function/

//#include "stdafx.h"
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <igl/igl_inline.h>
//#include <igl/matrix_to_list.h>
#include <igl/per_vertex_normals.h>
#include <igl/ply.h>
#include <igl/readPLY.h>
#include <igl/readMESH.h>
#include <igl/readOFF.h>

#include "tutorial_shared_path.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
//#include <float.h>
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
	Matrix<double, Dynamic, 3> temp;

	igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", V, F);
	//igl::readOBJ(TUTORIAL_SHARED_PATH "/bumpy-cube.obj", V, F);
	//igl::readOBJ(TUTORIAL_SHARED_PATH "/arm.obj", V, F);
	//igl::readMESH(TUTORIAL_SHARED_PATH "/octopus-high.mesh", V, temp,F);
	//igl::readMESH(TUTORIAL_SHARED_PATH "/big-sigcat.mesh", V, temp, F);
	//igl::readMESH(TUTORIAL_SHARED_PATH "/hand.mesh", V, temp, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/lion.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/camelhead.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);

	igl::per_vertex_normals(V, F, N);
	igl::barycenter(V, F, B);
	
	cout << "Welcome to octreeSDF, an octree implementation of the Surface Diameter Function." << endl;
	cout << "This is still under development, TODO: very big meshes, visualization." << endl;
	cout << "The mesh has V.rows(): " << V.rows() << ", F.rows(): " << F.rows() << ", B.rows():" << B.rows() << endl;

	vector<array<double,3>> vecV = Eigen2CPP(V);
	vector<array<int, 3>> vecF = Eigen2CPP(F);
	vector<array<double, 3>> vecN = Eigen2CPP(N);
	vector<array<double, 3>> vecB = Eigen2CPP(B);

	double total1 = 0;
	for (int i = 0; i < 30; ++i) {
		chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
		SDF tree(vecV, vecF, vecB);
		tree.init();
		tree.build();
		chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
		chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
		total1 += time_span.count();
	}
	total1 /= 30;
	cout << "Average running time (30 trials) for tree.build() is " << total1 << " seconds." << endl;

	//Doing the queries
	for (size_t i = 0; i < vecN.size(); ++i) {
		vecN.at(i)[0] = -vecN.at(i)[0];
		vecN.at(i)[1] = -vecN.at(i)[1];
		vecN.at(i)[2] = -vecN.at(i)[2];
	}
	SDF tree(vecV, vecF, vecB);
	tree.init();
	tree.build();
	vector<double> sdf = tree.query(vecN);

	//Adjust the values for visualization
	double mymax = *max_element(sdf.begin(), sdf.end());
	double mean = 0;
	for (size_t i = 0; i < sdf.size(); ++i) {
		mean += sdf.at(i);
	}
	mean /= sdf.size();
	for (size_t i = 0; i < sdf.size(); ++i) {
		sdf.at(i) = log(sdf.at(i) / mymax);
		//sdf.at(i) = log(2*sdf.at(i)+2*mean / mymax); //for octopus, sig-cat, fertility, hand mesh, lion, cheburashka, camelhead
	}

	double* ptr = &sdf[0];
	Map<VectorXd> Z(ptr, sdf.size());
	MatrixXd C(V.rows(), 3);
	igl::jet(Z, true, C);

	//// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	//// Launch the viewer

	viewer.data().set_colors(C);
	viewer.data().show_lines = !viewer.data().show_lines;
	viewer.launch();
}