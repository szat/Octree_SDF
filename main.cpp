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
	igl::per_vertex_normals(V, F, N);
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

	/*
	double total = 0;
	for (int i = 0; i < 30; ++i) {
		chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
		SDF tree(vecV, vecF, vecB);
		tree.init();
		tree.build();
		chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
		chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
		total += time_span.count();
	}
	total /= 30;
	cout << "Average running time (30 trials) for tree.build() is " << total << " seconds." << endl;*/

	for (size_t i = 0; i < vecN.size(); ++i) {
		vecV.at(i)[0] = -vecV.at(i)[0];
		vecV.at(i)[1] = -vecV.at(i)[1];
		vecV.at(i)[2] = -vecV.at(i)[2];
	}
	SDF tree(vecV, vecF, vecB);
	tree.init();
	tree.build();
	vector<double> sdf;
	int counter = 0;
	for (size_t i = 0; i < vecV.size(); ++i) {
		vector<array<double, 3>> intersect = tree.query(vecV.at(i), vecN.at(i));
		set<double> lengths;
		//cout << "penetration distances ";
		for (size_t j = 0; j < intersect.size(); ++j) {
			array<double, 3> vec = { intersect.at(j)[0] - vecV.at(i)[0], intersect.at(j)[1] - vecV.at(i)[1],intersect.at(j)[2] - vecV.at(i)[2] };
			double norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
			//cout << norm << ", " << endl;
			lengths.insert(sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
		}
		if (lengths.size() > 0) {
			sdf.push_back(*lengths.begin());
		}
		else {
			sdf.push_back(0);
			counter++;
		}
		//std::set<double>::iterator it;
		//for (it = lengths.begin(); it != lengths.end(); ++it) {
		//	cout << *it << " "; // Note the "*" here
		//}
		//cout << endl;
	}

	cout << "length of sdf " << sdf.size() << endl;
	array<double, 3> source = { -1000,0,0 };
	array<double, 3> dir = {1, 0, 0};
	//tree.test();
	vector<array<double,3>> inter = tree.query(source, dir);
	double mymax = *max_element(sdf.begin(), sdf.end());
	for (size_t i = 0; i < sdf.size(); ++i) {
		sdf.at(i) = sdf.at(i) / mymax;
	}
	cout << "intersection size " << inter.size() << endl;
	cout << "Let's merrily code ahead!" << endl;
	cin.ignore();

	double* ptr = &sdf[0];
	Map<VectorXd> Z(ptr, sdf.size());
	MatrixXd C(V.rows(), 3);
	igl::jet(Z, true, C);

	//// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	//// Launch the viewer

	viewer.data().set_colors(C);
	viewer.launch();

}