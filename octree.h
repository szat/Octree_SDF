#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <array>
#include <cmath>

using namespace std;

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

bool TriangleBox(const array<array<double, 3>, 3> & triangle, const array<array<double, 3>, 2> & box) {
	//	this->box[0] = bound_down;
	//  this->box[1] = bound_up;
	for (size_t i = 0; i < 3; ++i) {
		array<double, 3> pt = triangle[i];
		if (box[0][0] > pt[0] || pt[0] > box[1][0]) {
			return false;
		}
		if (box[0][1] > pt[1] || pt[1] > box[1][1]) {
			return false;
		}
		if (box[0][2] > pt[2] || pt[2] > box[1][2]) {
			return false;
		}
	}
	return true;
}

class SDF {
private:
	//Data
	vector<array<double, 3>> & V;
	vector<array<int, 3>> & F;
	vector<array<double, 3>> & bary;

	//Node
	int leaves;
	array<array<double, 3>, 2> box;
	vector<int> indices; //of triangles
	array<unique_ptr<SDF>, 8> children;

	//For testing
	bool is_leaf() const;
	bool test1() const;
	int test2() const;
	vector<int> test3() const;
	vector<int> test4() const;
public:
	SDF(vector<array<double, 3>> & V, vector<array<int, 3>> & F, vector<array<double, 3>> & bary)
		: V(V), F(F), bary(bary) {};
	SDF(const SDF& other) = delete;
	SDF& operator=(const SDF& rhs) = delete;

	void init();
	void build();
	void test() const;
	vector<array<double, 3>> query(array<double, 3> & source, array<double, 3> & dir) const;
	vector<double> query(vector<array<double, 3>> & v_normals) const;
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
	this->leaves = indices.size();
	//cout << "leaves nb " << this->leaves << endl;

	//Compute box
	array<double, 3> bound_down;
	array<double, 3> bound_up;
	double epsilon = 0.0000000001;
	array<double, 3> first;
	first = this->V.at(this->F.at(this->indices.at(0))[0]);
	bound_up = { first[0] + epsilon, first[1] + epsilon, first[2] + epsilon };
	bound_down = { first[0] - epsilon, first[1] - epsilon, first[2] - epsilon };
	//cout << "bound up : " << bound_up[0] << ", " << bound_up[1] << ", " << bound_up[2] << endl;
	//cout << "bound down : " << bound_down[0] << ", " << bound_down[1] << ", " << bound_down[2] << endl;

	for (size_t i = 0; i < this->indices.size(); ++i) {
		for (size_t j = 0; j < 3; ++j) {
			array<double, 3> pt = this->V.at(this->F.at(this->indices.at(i))[j]);
			for (size_t d = 0; d < 3; ++d) {
				if (bound_up[d] < pt[d]) {
					bound_up[d] = pt[d];
				}
				if (bound_down[d] > pt[d]) {
					bound_down[d] = pt[d];
				}
			}
		}
	}
	bound_up = { bound_up[0] + epsilon, bound_up[1] + epsilon, bound_up[2] + epsilon };
	bound_down = { bound_down[0] - epsilon, bound_down[1] - epsilon, bound_down[2] - epsilon };
	//cout << "bound up : " << bound_up[0] << ", " << bound_up[1] << ", " << bound_up[2] << endl;
	//cout << "bound down : " << bound_down[0] << ", " << bound_down[1] << ", " << bound_down[2] << endl;

	this->box[0] = bound_down;
	this->box[1] = bound_up;

	//More than one triangle remaining in the box
	if (this->indices.size() > 1) {
		//Compute center
		array<double, 3> mean;
		int nb = this->indices.size();
		mean[0] = bary.at(this->indices.at(0))[0];
		mean[1] = bary.at(this->indices.at(0))[1];
		mean[2] = bary.at(this->indices.at(0))[2];
		for (size_t i = 1; i < nb; ++i) {
			mean[0] = mean[0] + (bary.at(this->indices.at(i))[0] - mean[0]) / (i + 1);
			mean[1] = mean[1] + (bary.at(this->indices.at(i))[1] - mean[1]) / (i + 1);
			mean[2] = mean[2] + (bary.at(this->indices.at(i))[2] - mean[2]) / (i + 1);
		}
		//cout << "center " << mean[0] << ", " << mean[1] << ", " << mean[2] << endl;

		//Compute new indices
		array<vector<int>, 8> sub_indices;

		for (size_t i = 0; i < this->indices.size(); ++i) {
			array<double, 3> pt = bary.at(indices.at(i));
			if (pt[0] > mean[0]) {
				if (pt[1] > mean[1]) {
					if (pt[2] > mean[2]) {
						sub_indices[0].push_back(indices.at(i));
					}
					else {//pt[2] <= mean[2] 
						sub_indices[1].push_back(indices.at(i));
					}
				}
				else {   //pt[1] <= mean[1]
					if (pt[2] > mean[2]) {
						sub_indices[2].push_back(indices.at(i));
					}
					else {//pt[2] <= mean[2]
						sub_indices[3].push_back(indices.at(i));
					}
				}
			}
			else {		 //pt[0] <= mean[0]
				if (pt[1] > mean[1]) {
					if (pt[2] > mean[2]) {
						sub_indices[4].push_back(indices.at(i));
					}
					else {//pt[2] <= mean[2]
						sub_indices[5].push_back(indices.at(i));
					}
				}
				else {   //pt[1] <= mean[1]
					if (pt[2] > mean[2]) {
						sub_indices[6].push_back(indices.at(i));
					}
					else //pt[2] <= mean[2]
						sub_indices[7].push_back(indices.at(i));
				}
			}
		}
		indices.clear(); //free memory

		//Create children and recurse
		for (size_t i = 0; i < 8; ++i) {
			if (sub_indices[i].size() > 0) {
				this->children[i] = unique_ptr<SDF>(new SDF(this->V, this->F, this->bary));
				this->children[i]->indices = sub_indices[i];
				this->children[i]->build();
			}
			else {
				this->children[i] = nullptr;
			}
		}
	}
	return;
}

vector<array<double, 3>> SDF::query(array<double, 3> & source, array<double, 3> & dir) const {
	//Returns all the intersection points of the mesh with the ray defined by source and dir
	if (this->is_leaf()) {
		array<double, 3> intersection;
		array<array<double, 3>, 3> tri;

		int tri_idx = (this->indices).at(0);
		array<int, 3> pts_idx = (this->F).at(tri_idx);
		tri[0] = (this->V).at(pts_idx[0]);
		tri[1] = (this->V).at(pts_idx[1]);
		tri[2] = (this->V).at(pts_idx[2]);

		bool is_hit = RayTriangle(source, dir, tri, intersection);
		if (is_hit == true) {
			vector<array<double, 3>> out = { intersection };
			return out;
		}
		else {
			vector<array<double, 3>> empty;
			return empty;
		} //can add epsilon wiggle after
	}
	else {
		vector<array<double, 3>> out;
		for (size_t i = 0; i < 8; ++i) {
			if (this->children[i] != nullptr) { 
				if (true == RayBox(source, dir, 0, DBL_MAX, this->children[i]->box)) { //leap of faith
					vector<array<double,3>> temp = this->children[i]->query(source, dir);
					out.insert(end(out), begin(temp), end(temp));
				}
			}
		}
		return out; //may be empty
	}
};

vector<double> SDF::query(vector<array<double, 3>> & inv_normals) const {
	//Needs normals that point outside->in
	//Returns the SDF for all the vertices of the mesh, needs the input of the normals
	//NORM SQUARED
	vector<double> sdf;
	for (size_t i = 0; i < inv_normals.size(); ++i) {
		vector<array<double, 3>> intersect = this->query((this->V).at(i), inv_normals.at(i));
		double my_min = DBL_MAX;
		for (size_t j = 0; j < intersect.size(); ++j) {
			array<double, 3> diff = { intersect.at(j)[0] - (this->V).at(i)[0], intersect.at(j)[1] - (this->V).at(i)[1],intersect.at(j)[2] - (this->V).at(i)[2] };
			double norm_sq = dot(diff, diff);
			if (my_min > norm_sq && norm_sq > 0.0005) {
				my_min = norm_sq;
			}
		}
		if (my_min != DBL_MAX) {
			sdf.push_back(my_min);
		}
		else {
			sdf.push_back(0); //in case there was "no interection"
		}
	}
	//Take out if not visualizing
	//double my_max = *max_element(sdf.begin(), sdf.end());
	//for (size_t i = 0; i < inv_normals.size(); ++i) {
	//	if (sdf.at(i) == 0) {
	//		sdf.at(i) = my_max;
	//	}
	//}
	return sdf;
}

bool SDF::is_leaf() const {
	bool is_leaf = true;
	for (int i = 0; i < 8; ++i) {
		if (this->children[i] != nullptr) {
			is_leaf = false;
		}
	}
	return is_leaf;
}

void SDF::test() const {
	cout << "Beginnig to test the octree..." << endl;
	bool t1 = this->test1();
	if (t1 == true) {
		cout << "->test1(): passed, all leaves have exactly one triangle index." << endl;
	}
	else {
		cout << "->test1(): failed, some leaves do not have exactly one triangle index." << endl;
	}

	int t2 = this->test2();
	if (t2 == this->F.size()) {
		cout << "->test2(): passed, triangle nb == leaf nb." << endl;
	}
	else {
		cout << "->test2(): failed, triangle nb != leaf nb" << endl;
	}

	vector<int> idx_o;
	for (size_t i = 0; i < this->F.size(); ++i) {
		idx_o.push_back(i);
	}
	set<int> idxset_orig(idx_o.begin(), idx_o.end());

	vector<int> t3 = this->test3();
	set<int> idxset3(t3.begin(), t3.end());
	if (idxset3 == idxset_orig) {
		cout << "->test3(): passed, all the triangles are contained in the current box." << endl;
	}
	else {
		cout << "->test3(): failed, a triangle is not contained in the current box." << endl;
	}

	vector<int> t4 = this->test4();
	set<int> idxset(t4.begin(), t4.end());
	if (idxset == idxset_orig) {
		cout << "->test4(): passed, the indices in the leaves equals the original indices." << endl;
	}
	else {
		cout << "->test4(): failed, the indices in the leaves does not equal the original indices." << endl;
	}
}

bool SDF::test1() const {
	//Verify that each leaf has only one triangle
	if (this->is_leaf() == true) { //in a leaf
		if (this->indices.size() == 1) {
			return true;
		}
		else {
			cout << "This leaf has " << this->indices.size() << " triangle index size." << endl;
			return false;
		}
	}
	else { //not a leaf
		bool all_ok = true;
		for (int i = 0; i < 8; ++i) {
			if (this->children[i] != nullptr) {
				all_ok = all_ok & this->children[i]->test1();
			}
		}
		return all_ok;
	}
}

int SDF::test2() const {
	//Verify that the number of leaves is the number of total triangles
	if (this->is_leaf() == true) {
		return 1;
	}
	else {
		int count = 0;
		bool all_ok = true;
		for (int i = 0; i < 8; ++i) {
			if (this->children[i] != nullptr) {
				count += this->children[i]->test2();
			}
		}
		return count;
	}
}

vector<int> SDF::test3() const {
	//Verify that the boxes do contain the triangles at every level
	if (this->is_leaf()) { //return the index set only if the indexed triangle is contained in the box
		int tri_idx = (this->indices).at(0);
		array<array<double, 3>, 3> triangle;
		array<int, 3> pts_idx = (this->F).at(tri_idx);
		triangle[0] = (this->V).at(pts_idx[0]);
		triangle[1] = (this->V).at(pts_idx[1]);
		triangle[2] = (this->V).at(pts_idx[2]);
		if (TriangleBox(triangle, this->box) == true) {
			return this->indices;
		}
		else {
			vector<int> empty;
			return empty;
		}
	}
	else { //not a leaf
		vector<int> idx;
		for (size_t i = 0; i < 8; ++i) {
			if (this->children[i] != nullptr) {
				vector<int> temp = this->children[i]->test3(); //<===== recurse
				idx.insert(end(idx), begin(temp), end(temp)); //append indices of the subtrees
			}
		}
		bool all_ok = true;
		for (size_t i = 0; i < idx.size(); ++i) {
			int tri_idx = idx.at(i);
			array<array<double, 3>, 3> triangle;
			array<int, 3> pts_idx = (this->F).at(tri_idx);
			triangle[0] = (this->V).at(pts_idx[0]);
			triangle[1] = (this->V).at(pts_idx[1]);
			triangle[2] = (this->V).at(pts_idx[2]);
			if (TriangleBox(triangle, this->box) == false) {
				all_ok = false;
			}
		}
		if (all_ok == true) {
			return idx; //return all the indices only if each triangle indexed was in the current box
		}
		else { //at least one triangle indexed was not contained in the current box, thus critical fail
			vector<int> empty;
			return empty;
		}
	}
}

vector<int> SDF::test4() const {
	//Verify that the union of the indices in the leaves is the total index set
	if (this->is_leaf()) {
		return this->indices;
	}
	else {
		vector<int> idx;
		for (int i = 0; i < 8; ++i) {
			if (this->children[i] != nullptr) {
				vector<int> temp = this->children[i]->test4();
				idx.insert(end(idx), begin(temp), end(temp));
			}
		}
		return idx;
	}
}
