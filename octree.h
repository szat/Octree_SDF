#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <array>

using namespace std;

class SDF {
private:
	//Data
	const vector<array<double, 3>> & V;
	const vector<array<int, 3>> & F;
	const vector<array<double, 3>> & bary;

	//Node
	array<array<double, 3>, 2> box;
	vector<int> indices; //of triangles
	array<unique_ptr<SDF>, 8> children;
public:
	SDF(const vector<array<double, 3>> & V, const vector<array<int, 3>> & F, const vector<array<double, 3>> & bary)
		: V(V), F(F), bary(bary) {};
	SDF(const SDF& other) = delete;
	SDF& operator=(const SDF& rhs) = delete;

	void init();
	void build();
};

void SDF::init() {
	for (size_t i = 0; i < this->F.size(); ++i) {
		this->indices.push_back((int)i);
	}
}

void SDF::build() {
	if (this->indices.size() == 0) {
		return;
	}

	//Compute box
	array<double, 3> bound_down;
	array<double, 3> bound_up;
	double epsilon = 0.0000000001;
	array<double, 3> first;
	first = this->V.at(this->F.at(this->indices.at(0))[0]);
	bound_up = { first[0] + epsilon, first[1] + epsilon, first[2] + epsilon };
	bound_down = { first[0] - epsilon, first[1] - epsilon, first[2] - epsilon };

	for (size_t i = 0; i < this->indices.size(); ++i) {
		for (size_t j = 0; j < 3; ++j) {
			array<double, 3> pt = this->V.at(this->F.at(this->indices.at(i))[j]);
			for (size_t d = 0; d < 3; ++d) {
				if (bound_up[d] < pt[d])
					bound_up[d] = pt[d];
				if (bound_down[d] > pt[d])
					bound_down[d] = pt[d];
			}
		}
	}
	bound_up = { bound_up[0] + epsilon, bound_up[1] + epsilon, bound_up[2] + epsilon };
	bound_down = { bound_down[0] - epsilon, bound_down[1] - epsilon, bound_down[2] - epsilon };
	this->box = { bound_down, bound_up };

	//More than one triangle remaining in the box
	if (this->indices.size() > 1) {
		//Compute center
		array<double, 3> mean;
		for (size_t i = 0; i < this->indices.size(); ++i) {
			mean[0] += mean[0] + bary.at(this->indices.at(i))[0];
			mean[1] += mean[1] + bary.at(this->indices.at(i))[1];
			mean[2] += mean[2] + bary.at(this->indices.at(i))[2];
		}
		mean[0] = mean[0] / this->indices.size();
		mean[1] = mean[1] / this->indices.size();
		mean[2] = mean[2] / this->indices.size();

		//Compute new indices
		array<vector<int>, 8> sub_indices;

		for (size_t i = 0; i < this->indices.size(); ++i) {
			array<double, 3> pt = bary.at(indices.at(i));
			if (pt[0] > mean[0]) {
				if (pt[1] > mean[1]) {
					if (pt[2] > mean[2])
						sub_indices[0].push_back(indices.at(i));
					else //pt[2] <= mean[2]
						sub_indices[1].push_back(indices.at(i));
				}
				else {   //pt[1] <= mean[1]
					if (pt[2] > mean[2])
						sub_indices[2].push_back(indices.at(i));
					else //pt[2] <= mean[2]
						sub_indices[3].push_back(indices.at(i));
				}
			}
			else {		 //pt[0] <= mean[0]
				if (pt[1] > mean[1]) {
					if (pt[2] > mean[2])
						sub_indices[4].push_back(indices.at(i));
					else //pt[2] <= mean[2]
						sub_indices[5].push_back(indices.at(i));
				}
				else {   //pt[1] <= mean[1]
					if (pt[2] > mean[2])
						sub_indices[6].push_back(indices.at(i));
					else //pt[2] <= mean[2]
						sub_indices[7].push_back(indices.at(i));
				}
			}
		}

		//Create children and recurse
		for (size_t i = 0; i < 8; ++i) {
			if (sub_indices[i].size() > 0) {
				this->children[i] = unique_ptr<SDF>(new SDF(this->V, this->F, this->bary));
				this->children[i]->indices = sub_indices[i];
				this->children[i]->build();
			}
		}
	}
	return;
}

double dot(const array<double, 3> &A, const array<double, 3> &B) {
	double out = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
	return out;
}

array<double, 3> cross(const array<double, 3> &A, const array<double, 3> &B) {
	array<double, 3> out;
	out[0] = A[1] * B[2] - A[2] * B[1];
	out[1] = A[2] * B[0] - A[0] * B[2];
	out[2] = A[0] * B[1] - A[1] * B[0];
	return out;
}

bool RayTriangle(const array<double, 3> & source, const array<double, 3> & dir, const array<array<double, 3>, 3> & tri, array<double, 3> & intersection) {
	//https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
	const double epsilon = 0.00000001;
	array<double, 3> v0 = tri[0];
	array<double, 3> v1 = tri[1];
	array<double, 3> v2 = tri[2];
	array<double, 3> h, s, q;
	double a, f, u, v;
	array<double, 3> edge1 = { v1[0] - v0[0], v1[1] - v0[1] , v1[2] - v0[2] };
	array<double, 3> edge2 = { v2[0] - v0[0], v2[1] - v0[1] , v2[2] - v0[2] };
	h = cross(dir, edge2);
	a = dot(edge1, h);
	if (a > -epsilon && a < epsilon) {
		//DOES NOT HANDLE CASE WHEN LINE IS PARALLEL TO TRIANGLE PLANE
		return false;
	}
	f = 1 / a;
	s = { source[0] - v0[0], source[1] - v0[1], source[2] - v0[2] };
	u = f * dot(s, h);
	if (u < 0.0 - epsilon || u > 1.0 + epsilon)
		return false;
	q = cross(s, edge1);
	v = f * dot(dir, q);
	if (v < 0.0 - epsilon || u + v > 1.0 + epsilon)
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

bool RayBox(const array<double, 3> & source, const array<double, 3> & dir, double low_t, double high_t, const array<array<double, 3>, 2> & box) {
	//https://www.cs.utah.edu/~awilliam/box/box.pdf
	array<double, 3> inv_dir = { 1 / dir[0], 1 / dir[1], 1 / dir[2] };
	array<int, 3> sign;
	sign[0] = inv_dir[0] < 0;
	sign[1] = inv_dir[1] < 0;
	sign[2] = inv_dir[2] < 0;
	double tmin, tmax, tymin, tymax, tzmin, tzmax;
	tmin = (box[sign[0]][0] - source[0]) * inv_dir[0];
	tmax = (box[1 - sign[0]][0] - source[0]) * inv_dir[0];
	tymin = (box[sign[1]][1] - source[1]) * inv_dir[1];
	tymax = (box[1 - sign[1]][1] - source[1]) * inv_dir[1];
	if ((tmin > tymax) || (tymin > tmax))
		return false;
	if (tymin > tmin)
		tmin = tymin;
	if (tymax < tmax)
		tmax = tymax;
	tzmin = (box[sign[2]][2] - source[2]) * inv_dir[2];
	tzmax = (box[1 - sign[2]][2] - source[2]) * inv_dir[2];
	if ((tmin > tzmax) || (tzmin > tmax))
		return false;
	if (tzmin > tmin)
		tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;
	return ((tmin < high_t) && (tmax > low_t));
}