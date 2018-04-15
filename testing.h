#pragma once
#include <random>
#include <chrono>
#include <ctime>
#include "octree.h"

using namespace std;
//
//void test_Eigen2CPP() {
//	/*
//	This function tests Eigen2CPP which serves as a bridge between Eigen and this C++ header project.
//	Input: an Eigen matrix, that NEEDS to be dynamic rows and with 3 cols. The type is templated.
//	Output: a C++ vector of arrays with 3 columns. The type is templated.
//	*/
//	Matrix<int, Dynamic, 3> F = (MatrixXi(4, 3) << 1, 2, 3,
//		1, 2, 3,
//		5, 6, 7,
//		7, 8, 9).finished();
//	vector<array<int, 3>> matF = Eigen2CPP(F);
//	for (int row = 0; row < F.rows(); ++row) {
//		for (int col = 0; col < 3; ++col) {
//			cout << matF.at(row)[col] << " ";
//		}
//		cout << endl;
//	}
//	cin.ignore();
//}

void test_RayBox() {
	cout << "Int test RayBox ..." << endl;
	double epsilon = 0.00000001;
	array<array < double, 3>, 2> box;
	box[0] = { 0,0,0 };
	box[1] = { 2,2,2 };
	array<double, 3> source;
	array<double, 3> dir;
	double low_t = 0;
	double high_t = DBL_MAX;

	//Test no intersection
	source = { 1,1,1 };
	dir = { 1,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] <<"," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 0,0,0 };
	dir = { 1,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { -1,-1,0 };
	dir = { 1,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 1,0,0 };
	dir = { 1,0,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 0,1,0 };
	dir = { 0,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 0,0,1 };
	dir = { 0,0,1 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 2,2,2 };
	dir = { 1,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { -2,0,0 };
	dir = { 1,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	box[0] = { box[0][0] - epsilon, box[0][1] - epsilon, box[0][2] - epsilon };
	box[1] = { box[1][0] + epsilon, box[1][1] + epsilon, box[1][2] + epsilon };

	source = { 0,1,0 };
	dir = { 0,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 0,0,1 };
	dir = { 0,0,1 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 2,2,2 };
	dir = { 1,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { -2,0,0 };
	dir = { 1,1,0 };
	cout << "Box = [" << box[0][0] << "," << box[0][1] << "," << box[0][2] << "],[" << box[1][0] << "," << box[1][1] << "," << box[1][2] << "]; ";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayBox(source, dir, low_t, high_t, box)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(-1000, 1000);
	double total = 0;
	for (int i = 0; i < 1000; ++i) {
		box[0] = { (double)dis(gen), (double)dis(gen), (double)dis(gen) };
		box[1] = { (double)dis(gen), (double)dis(gen), (double)dis(gen) };
		source = { (double)dis(gen), (double)dis(gen), (double)dis(gen) };
		dir =    { (double)dis(gen), (double)dis(gen), (double)dis(gen) };
		chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
		RayBox(source, dir, low_t, high_t, box);
		chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
		chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
		total += time_span.count();
	}
	total /= 1000;
	cout << "Average running time for RayBox is " << total << " seconds." << endl;
}

void test_RayTriangle() {
	cout << "Int test RayBox ..." << endl;
	double epsilon = 0.00000001;
	array<array < double, 3>, 3> tri;
	tri[0] = { 0,0,0 };
	tri[1] = { 3,0,0 };
	tri[2] = { 0,3,0 };

	array<double, 3> source;
	array<double, 3> dir;
	array<double, 3> intersect;

	//Test
	source = { 0,0,0 };
	dir = { 1,0,0 };
	cout << "Tri = [" << tri[0][0] << "," << tri[0][1] << "," << tri[0][2] << "],[" << tri[1][0] << "," << tri[1][1] << "," << tri[1][2] << "],[" << tri[2][0] << "," << tri[2][1] << "," << tri[2][2] << "]";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayTriangle(source, dir, tri, intersect)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { -1,0,0 };
	dir = { 1,0,0 };
	cout << "Tri = [" << tri[0][0] << "," << tri[0][1] << "," << tri[0][2] << "],[" << tri[1][0] << "," << tri[1][1] << "," << tri[1][2] << "],[" << tri[2][0] << "," << tri[2][1] << "," << tri[2][2] << "]";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayTriangle(source, dir, tri, intersect)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 0,0,0 };
	dir = { 1,1,0 };
	cout << "Tri = [" << tri[0][0] << "," << tri[0][1] << "," << tri[0][2] << "],[" << tri[1][0] << "," << tri[1][1] << "," << tri[1][2] << "],[" << tri[2][0] << "," << tri[2][1] << "," << tri[2][2] << "]";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayTriangle(source, dir, tri, intersect)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { -1,0,0 };
	dir = { -1,-1,0 };
	cout << "Tri = [" << tri[0][0] << "," << tri[0][1] << "," << tri[0][2] << "],[" << tri[1][0] << "," << tri[1][1] << "," << tri[1][2] << "],[" << tri[2][0] << "," << tri[2][1] << "," << tri[2][2] << "]";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayTriangle(source, dir, tri, intersect)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 1,1,1 };
	dir = { 0,0,-1 };
	cout << "Tri = [" << tri[0][0] << "," << tri[0][1] << "," << tri[0][2] << "],[" << tri[1][0] << "," << tri[1][1] << "," << tri[1][2] << "],[" << tri[2][0] << "," << tri[2][1] << "," << tri[2][2] << "]";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayTriangle(source, dir, tri, intersect)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { 1,1,0 };
	dir = { 0,0,1 };
	cout << "Tri = [" << tri[0][0] << "," << tri[0][1] << "," << tri[0][2] << "],[" << tri[1][0] << "," << tri[1][1] << "," << tri[1][2] << "],[" << tri[2][0] << "," << tri[2][1] << "," << tri[2][2] << "]";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayTriangle(source, dir, tri, intersect)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;

	source = { -1,-1,0 };
	dir = { 0,0,1 };
	cout << "Tri = [" << tri[0][0] << "," << tri[0][1] << "," << tri[0][2] << "],[" << tri[1][0] << "," << tri[1][1] << "," << tri[1][2] << "],[" << tri[2][0] << "," << tri[2][1] << "," << tri[2][2] << "]";
	cout << "Ray, source = [" << source[0] << "," << source[1] << "," << source[2] << "], dir = [" << dir[0] << "," << dir[1] << "," << dir[2] << "]; ";
	if (true == RayTriangle(source, dir, tri, intersect)) cout << "RayBox(source, dir, box) == true" << endl;
	else cout << "RayBox(source, dir, box) == false" << endl;


}

void test_build() {
	return;
}
