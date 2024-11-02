#include "aabb.h" 
#include <algorithm>
// @@@@@@ VOTRE CODE ICI
// Implémenter l'intersection d'un rayon avec un AABB dans l'intervalle décrit.
bool AABB::intersect(Ray ray, double t_min, double t_max)  {
	for (int i = 0; i < 3; ++i) { // Iterate over the x, y, and z axes
        double invD = 1.0 / ray.direction[i]; // Inverse of the ray direction for this axis
        double t0 = (min[i] - ray.origin[i]) * invD;
        double t1 = (max[i] - ray.origin[i]) * invD;

        // Swap t0 and t1 if needed, to ensure t0 is the entry point and t1 is the exit point
        if (invD < 0.0) std::swap(t0, t1);

        // Update the interval
        t_min = std::max(t_min, t0);
        t_max = std::min(t_max, t1);

        // If the intervals do not overlap, there is no intersection
        if (t_max <= t_min) return false;
    }
    return true;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction qui permet de trouver les 8 coins de notre AABB.
std::vector<double3> retrieve_corners(AABB aabb) {
	std::vector<double3> boxCorners;
	double3 min = aabb.min;
	double3 max = aabb.max;
	for(int i = 0; i < 8; i++) {
		double a = (i / 4) % 2 == 0 ? min.x : max.x;
        double b = (i / 2) % 2 == 0 ? min.y : max.y;
        double c = (i % 2) == 0 ? min.z : max.z;

        boxCorners.push_back(double3{a, b, c});
	}

	return boxCorners;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.
AABB construct_aabb(std::vector<double3> points) {
	double3 max{-DBL_MAX,-DBL_MAX,-DBL_MAX};
	double3 min{DBL_MAX,DBL_MAX,DBL_MAX};
	for(const auto& point : points) {
		min.x = std::min(min.x,point.x);
		min.y = std::min(min.y,point.y);
		min.z = std::min(min.z,point.z);
		max.x = std::max(max.x,point.x);
		max.y = std::max(max.y,point.y);
		max.z = std::max(max.z,point.z);
	}
	return AABB{min,max};
};

AABB combine(AABB a, AABB b) {
	double3 min = {std::min(a.min.x, b.min.x), std::min(a.min.y, b.min.y), std::min(a.min.z, b.min.z)};
    double3 max = {std::max(a.max.x, b.max.x), std::max(a.max.y, b.max.y), std::max(a.max.z, b.max.z)};
    return AABB{min, max};
};

bool compare(AABB a, AABB b, int axis){
	return a.min[axis] < b.min[axis];
};