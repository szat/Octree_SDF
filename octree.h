#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <array>

using namespace std;

double dot(const array<double, 3> &A, const array<double, 3> &B) {
	double out = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
	return out;
}

array<double, 3> cross(const array<double, 3> &A, const array<double, 3> &B) {
	//a = A.y * B.z - A.z * B.y;
	//b = A.z * B.x - A.x * B.z;
	//c = A.x * B.y - A.y * B.x;
	array<double, 3> out;
	out[0] = A[1] * B[2] - A[2] * B[1];
	out[1] = A[2] * B[0] - A[0] * B[2];
	out[2] = A[0] * B[1] - A[1] * B[0];
	return out;
}

bool halfRayTriangleIntersect(const array<double,3> & source, const array<double,3> & dir, const array<double,9> & tri, array<double,3> & intersection) {
	//https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
	//This function returns a coordinate in 3d 
	//A vector of points in 3d
	const double epsilon = 0.0000001;
	array<double, 3> v0 = { tri[0], tri[1], tri[2] };
	array<double, 3> v1 = { tri[3], tri[4], tri[5] };
	array<double, 3> v2 = { tri[6], tri[7], tri[8] }; 
	array<double, 3> h, s, q;
	double a, f, u, v;
	array<double, 3> edge1 = { v1[0] - v0[0], v1[1] - v0[1] , v1[2] - v0[2] };
	array<double, 3> edge2 = { v2[0] - v0[0], v2[1] - v0[1] , v2[2] - v0[2] };
	h = cross(dir, edge2);
	a = dot(edge1, h);
	if (a > -epsilon && a < epsilon)
		return false;
	f = 1 / a;
	s = { source[0] - v0[0], source[1] - v0[1], source[2] - v0[2] };
	u = f * dot(s, h);
	if (u < 0.0 || u > 1.0)
		return false;
	q = cross(s, edge1);
	v = f * dot(dir, q);
	if (v < 0.0 || u + v > 1.0)
		return false;
	// At this stage we can compute t to find out where the intersection point is on the line.
	double t = f * dot(edge2, q);
	if (t > epsilon) { //ray intersection
		intersection[0] = source[0] + t * dir[0];
		intersection[1] = source[1] + t * dir[1];
		intersection[2] = source[2] + t * dir[2];
		return true;
	}
	else
		return false;
}

bool halfRayBoxIntersect(const array<double,3> & source, const array<double,3> & dir, double low_t, double high_t, const array<array<double,3>,2> & box) {
	//https://www.cs.utah.edu/~awilliam/box/box.pdf
	array<double, 3> inv_dir = { 1 / dir[0], 1 / dir[1], 1 / dir[2] };
	array<int, 3> sign;
	sign[0] = inv_dir[0] < 0;
	sign[1] = inv_dir[1] < 0;
	sign[2] = inv_dir[2] < 0;

	double tmin, tmax, tymin, tymax, tzmin, tzmax;
	tmin =  (box[sign[0]][0] -   source[0]) * inv_dir[0];
	tmax =  (box[1-sign[0]][0] - source[0]) * inv_dir[0];
	tymin = (box[sign[1]][1] -   source[1]) * inv_dir[1];
	tymax = (box[1 - sign[1]][1] - source[1]) * inv_dir[1];
	if ((tmin > tymax) || (tymin > tmax))
		return false;
	if (tymin > tmin)
		tmin = tymin;
	if (tymax < tmax)
		tmax = tymax;
	tzmin = (box[sign[2]][2] -   source[2]) * inv_dir[2];
	tzmax = (box[1-sign[2]][2] - source[2]) * inv_dir[2];
	if ((tmin > tzmax) || (tzmin > tmax))
		return false;
	if (tzmin > tmin)
		tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;
	return ((tmin < high_t) && (tmax > low_t));

	return true;
}