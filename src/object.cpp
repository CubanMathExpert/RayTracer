#include "object.h"

// Fonction retournant soit la valeur v0 ou v1 selon le signe.
int rsign(double value, double v0, double v1) {
	return (int(std::signbit(value)) * (v1-v0)) + v0;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection d'une sphère.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Sphere::local_intersect(Ray ray, 
							 double t_min, double t_max, 
							 Intersection *hit) 
{

	double a = dot(ray.direction, ray.direction);
	double b = 2 * dot(ray.direction, ray.origin);
	double c = dot(ray.origin, ray.origin) - pow(radius, 2);
	double t;

	double discriminant = b * b - 4 * a * c;

	// no solutions
	if (discriminant < EPSILON) 
	{
		return false;
	}
	// les deux valeurs possibles de t
	double sqrt = std::sqrt(discriminant);
	double temp1 = (-b + sqrt) / (2 * a);
	double temp2 = (-b - sqrt) / (2 * a);

	if (temp1 > t_min && temp1 < t_max && temp2 > t_min && temp2 < t_max){
		t = std::min(temp1, temp2); // two sol
	}else if (temp1 <t_max && temp1 > t_min){
		t = temp1; // one sol 
	}else if (temp2 < t_max && temp2 > t_min){
		t = temp2; // one sol 
	}else{
		return false;
	}

	// update the hit information
	hit->depth = t;
	hit->position = ray.origin + t * ray.direction;
	hit->normal = normalize(hit->position);
	return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
	//c'est la construction d'un AABB pour une sphère
	double3 min = {-radius, -radius, -radius};
    double3 max = {radius, radius, radius};
    return AABB{min, max};
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un quad (rectangle).
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Quad::local_intersect(Ray ray, 
							double t_min, double t_max, 
							Intersection *hit)
{
	// Quad properties
	double3 normal = double3(0, 0, 1);
	double3 origin = double3(0, 0, 0);

	double denominator = dot(normal, ray.direction);
	// no collision possible
	if (denominator == 0)
	{
		//std::cout << "Ray is parallel to the quad" << std::endl;
		return false; // ray is parallel to the quad
	}
	else
	{
		double t = dot(normal, origin - ray.origin) / denominator;
		if (t < t_min || t > t_max)
		{
			return false;
		}

		double3 intersection = ray.origin + t * ray.direction;
		if (intersection.x < -half_size || intersection.x > half_size || intersection.y < -half_size || intersection.y > half_size)
		{
			return false;
		}

		hit->depth = t;
		hit->position = intersection;
		hit->normal = normal;
		hit->uv = double2((intersection.x + half_size) / (2 * half_size), (intersection.y + half_size) / (2 * half_size));
	}

	
	return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
	double epsilon = 1e-6;
	//c'est la construction d'un AABB pour un quad
	double3 min = {-half_size, -half_size, -epsilon}; 
	double3 max = {half_size, half_size, epsilon}; 
	return Object::compute_aabb();
	//return Object::compute_aabb();
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un cylindre.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Cylinder::local_intersect(Ray ray, 
							   double t_min, double t_max, 
							   Intersection *hit)
{
	// Cylinder properties
	double a = pow(ray.direction.x, 2) + pow(ray.direction.z, 2);
	double b = 2 * (ray.direction.x * ray.origin.x + ray.direction.z * ray.origin.z);
	double c = pow(ray.origin.x, 2) + pow(ray.origin.z, 2) - pow(radius, 2);
	double t;
	bool twoHits = false;
	// Calculate discriminant
    double discriminant = b * b - 4 * a * c;
    if (discriminant < EPSILON) {
        return false; // No real roots, ray misses the cylinder
    }

    // Calculate both roots of the quadratic equation
    double sqrtDiscriminant = sqrt(discriminant);
    double t1 = (-b - sqrtDiscriminant) / (2 * a);
    double t2 = (-b + sqrtDiscriminant) / (2 * a);

    // Ensure t1 is the smaller value
    if (t1 > t2) std::swap(t1, t2);

    // Check if t1 is within bounds and if the intersection is within the cylinder's height
    double y1 = ray.origin.y + t1 * ray.direction.y;
    if (t1 > t_min && t1 < t_max && y1 >= -half_height && y1 <= half_height) {
        hit->depth = t1;
        hit->position = ray.origin + t1 * ray.direction;
        hit->normal = normalize(double3(hit->position.x, 0, hit->position.z));
        return true;
    }

    // Check t2 if t1 was out of bounds or not within the height
    double y2 = ray.origin.y + t2 * ray.direction.y;
    if (t2 > t_min && t2 < t_max && y2 >= -half_height && y2 <= half_height) {
        hit->depth = t2;
        hit->position = ray.origin + t2 * ray.direction;
        hit->normal = normalize(double3(hit->position.x, 0, hit->position.z));
        return true;
    }

    return false; // No valid intersection within bounds and height
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il fautscr que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
	//c'est la construction d'un AABB pour un cylindre
	double3 min = {-radius, -radius, -half_height};
    double3 max = {radius, radius, half_height};
	return AABB{min, max};
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un mesh.
//
// Référez-vous au PDF pour la paramétrisation pour les coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
//
bool Mesh::local_intersect(Ray ray,  
						   double t_min, double t_max, 
						   Intersection* hit)
{

	bool hit_anything = false;
    Intersection temp_hit;

    for (auto const& tri : triangles) {
        if (intersect_triangle(ray, t_min, t_max, tri, &temp_hit)) {
            if (temp_hit.depth < t_max) {
                t_max = temp_hit.depth;
                *hit = temp_hit;
                hit_anything = true;
            }
        }
    }

    return hit_anything;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un triangle.
// S'il y a intersection, remplissez hit avec l'information sur la normale et les coordonnées texture.
bool Mesh::intersect_triangle(Ray  ray, 
							  double t_min, double t_max,
							  Triangle const tri,
							  Intersection *hit)
{
	// Extrait chaque position de sommet des données du maillage.
	double3 const &p0 = positions[tri[0].pi]; // ou Sommet A (Pour faciliter les explications)
	double3 const &p1 = positions[tri[1].pi]; // ou Sommet B
	double3 const &p2 = positions[tri[2].pi]; // ou Sommet C

	// Triangle en question. Respectez la convention suivante pour vos variables.
	//
	//     A
	//    / \
	//   /   \
	//  B --> C
	//
	// Respectez la règle de la main droite pour la normale.

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Pour plus de d'informations sur la géométrie, référez-vous à la classe dans object.hpp.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.
	// on va faire un plan par dessus le triangle 
	double3 normal = normalize(cross(p1 - p0, p2 - p0)); 
	double3 origin = {0, 0, 0};
	double3 intersection;

	double denominator = dot(normal, ray.direction);
	// do we intersect plane containing triangle ?
	// make a fake plane around the triangle that is of size 
	// 2 * length (p2 - p1) and 2 * length (p0 - p1)
	if (fabs(denominator) < EPSILON)
	{
		return false; // ray is parallel to the quad
	}

	double t = dot(normal, p0 - ray.origin) / denominator;
	if (t < t_min || t > t_max)
	{
		return false;
	}

	intersection = ray.origin + t * ray.direction;

	double3 edge0 = p1 - p0;
	double3 edge1 = p2 - p1;
	double3 edge2 = p0 - p2;

	double3 C0 = intersection - p0;
	double3 C1 = intersection - p1;
	double3 C2 = intersection - p2;

	if (dot(normal, cross(edge0, C0)) >= 0 &&
	   dot(normal, cross(edge1, C1)) >= 0 &&
	   dot(normal, cross(edge2, C2)) >= 0)
	{
		hit->depth = t;
		hit->position = intersection;
		hit->normal = normal;
		return true;
	}

	return false;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
	return construct_aabb(positions);
}