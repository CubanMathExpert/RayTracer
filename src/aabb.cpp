#include "aabb.h" 

// @@@@@@ VOTRE CODE ICI
// Implémenter l'intersection d'un rayon avec un AABB dans l'intervalle décrit.
bool AABB::intersect(Ray ray, double t_min, double t_max)  {
       return true;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction qui permet de trouver les 8 coins de notre AABB.
std::vector<double3> retrieve_corners(AABB aabb) {
	std::vector<double3> boxCorners;

	double3 min = aabb.min;
	double3 max = aabb.max;

	boxCorners[0] = double3{min.x, min.y, min.z};
	boxCorners[1] = double3{min.x, min.y, max.z};
	boxCorners[2] = double3{min.x, max.y, min.z};
	boxCorners[3] = double3{min.x, max.y, max.z};
	boxCorners[4] = double3{max.x, min.y, min.z};
	boxCorners[5] = double3{max.x, min.y, max.z};
	boxCorners[6] = double3{max.x, max.y, min.z};
	boxCorners[7] = double3{max.x, max.y, max.z};

	return boxCorners;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.
AABB construct_aabb(std::vector<double3> points) {
	return AABB{double3{-DBL_MAX,-DBL_MAX,-DBL_MAX},double3{DBL_MAX,DBL_MAX,DBL_MAX}};
};

AABB combine(AABB a, AABB b) {
	return AABB{min(a.min,b.min),max(a.max,b.max)};
};

bool compare(AABB a, AABB b, int axis){
	return a.min[axis] < b.min[axis];
};