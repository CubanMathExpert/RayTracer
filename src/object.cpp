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

	double phi = atan2(hit->position.x, hit->position.z);
	if (phi < 0) phi += 2 * M_PI;  

	double theta = acos(hit->position.y / radius);

	double u = phi / (2 * M_PI); 
	double v =  (theta / M_PI);

	hit->uv = double2(u, v);
	return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
	 
	AABB cylinder = construct_aabb({{-radius, -radius, -radius}, {radius,radius,radius}});
	std::vector<double3> local_corners = retrieve_corners(cylinder);
    double3 min, max;// global
    bool first_point = true;
    for (const auto& corner : local_corners) {
        double4 global_corner = mul(transform, {corner, 1});

        if (first_point) {
            min = {global_corner.x, global_corner.y, global_corner.z};
            max = min;
            first_point = false;
        } else {
            min = {std::min(min.x, global_corner.x),
                   std::min(min.y, global_corner.y),
                   std::min(min.z, global_corner.z)};
            max = {std::max(max.x, global_corner.x),
                   std::max(max.y, global_corner.y),
                   std::max(max.z, global_corner.z)};
        }
    }

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
	double3 normal = double3(0, 0, 1);
	double denominator = dot(normal, ray.direction);

	if (fabs(denominator) < EPSILON) {
		return false; 
	}

	double t = dot(normal, double3(0, 0, 0) - ray.origin) / denominator;
	if (t < t_min || t > t_max) {
		return false;
	}

	double3 intersection = ray.origin + t * ray.direction;
	if (fabs(intersection.x) > half_size || fabs(intersection.y) > half_size) {
		return false;
	}
	 // Adjust normal direction based on ray approach direction (making it double-sided)
    if (denominator > 0) {
        normal = -normal;  // Flip normal to face the ray
    }
	double u = (intersection.x + half_size) / (2 * half_size);
	double v = (intersection.y + half_size) / (2 * half_size);
	//update
	hit->depth = t;
	hit->position = intersection;
	hit->normal = normal;
	hit->uv = double2(u,1-v);

	return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
	double epsilon = 1e-6;

	// Define the four corners 
	double3 local_corners[4] = {
		{-half_size, -half_size, 0},
		{half_size, -half_size, 0},
		{half_size, half_size, 0},
		{-half_size, half_size, 0}
	};

	// Initialize first
	double4 first_corner_global = mul(transform, {local_corners[0], 1});
	double3 min = {first_corner_global.x, first_corner_global.y, first_corner_global.z};
	double3 max = min;

	// Transform each corner to global space and update min and max bounds
	for (int i = 1; i < 4; ++i) {
		double4 global_corner = mul(transform, {local_corners[i], 1});
		min = {std::min(min.x, global_corner.x),
		       std::min(min.y, global_corner.y),
		       std::min(min.z, global_corner.z)};
		max = {std::max(max.x, global_corner.x),
		       std::max(max.y, global_corner.y),
		       std::max(max.z, global_corner.z)};
	}

	// error bounds to avoid the ray to miss 
	min.z -= epsilon;
	max.z += epsilon;

	return AABB{min, max};
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
	double a = pow(ray.direction.x, 2) + pow(ray.direction.z, 2);
    double b = 2 * (ray.direction.x * ray.origin.x + ray.direction.z * ray.origin.z);
    double c = pow(ray.origin.x, 2) + pow(ray.origin.z, 2) - pow(radius, 2);
    double t;
    bool twoHits = false;

    double discriminant = b * b - 4 * a * c;
    // no solutions
    if (discriminant < EPSILON) {
        return false;
    }
    // one solution
    else if (discriminant == 0) {
        double temp = -b / (2 * a); // valeur de t pour la seule collision
        if (temp < t_max && temp > t_min) {
            t = temp;
        } else {
            return false;
        }
    }
    // two solutions
    else {
        // les deux valeurs possibles de t
        double temp1 = (-b + sqrt(discriminant)) / (2 * a);
        double temp2 = (-b - sqrt(discriminant)) / (2 * a);

        // les deux t sont validesersection with the inside of the cylinder
        if (temp1 > t_min && temp1 < t_max && temp2 > t_min && temp2 < t_max) {
            if (temp1 < temp2) {
                t = temp1;
            } else {
                t = temp2;
            }
        } else if (temp1 > t_min && temp1 < t_max) {
            t = temp1;
        } else if (temp2 > t_min && temp2 < t_max) {
            t = temp2;
        } else {
            return false;
        }
    }

    // check if the intersection point is within the cylinder's height
    double y = ray.origin.y + t * ray.direction.y;
    if (y >= -half_height && y < half_height) {
        hit->depth = t;
        hit->position = ray.origin + t * ray.direction;
        hit->normal = normalize(double3(hit->position.x, 0, hit->position.z));
        // Calculate UV coordinates (angle and height around y)
		double u = (atan2(hit->position.z, hit->position.x) + M_PI) / (2 * M_PI); 
		double v = 1-(hit->position.y + half_height) / (2 * half_height);         
		hit->uv = double2(u, v);
        return true;
    }

    // check for intersection with the inside of the cylinder
    double temp1 = (-b + sqrt(discriminant)) / (2 * a);
    double temp2 = (-b - sqrt(discriminant)) / (2 * a);
    if (temp1 > t_min && temp1 < t_max) {
        t = temp1;
    } else if (temp2 > t_min && temp2 < t_max) {
        t = temp2;
    } else {
        return false;
    }

    y = ray.origin.y + t * ray.direction.y;
    if (y >= -half_height && y < half_height) {
        hit->depth = t;
        hit->position = ray.origin + t * ray.direction;
        hit->normal = -normalize(double3(hit->position.x, 0, hit->position.z)); // Invert normal for inside
        // Calculate UV coordinates (and height around y)
		double u = (atan2(hit->position.z, hit->position.x) + M_PI) / (2 * M_PI); 
		double v = 1-(hit->position.y + half_height) / (2 * half_height);        
		hit->uv = double2(u, v);
        return true;
    }

    return false;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {

	AABB cylinder = construct_aabb({{-radius, -half_height, -radius}, {radius, half_height,radius}});
	std::vector<double3> local_corners = retrieve_corners(cylinder);
	//make sure it is initialised first
	double3 global_min = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    double3 global_max = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    
    for (const auto& corner : local_corners) {
        double4 transformed_corner = mul(transform, {corner, 1}); // Transform to world space

        global_min = {std::min(global_min.x, transformed_corner.x),
					  std::min(global_min.y, transformed_corner.y),
					  std::min(global_min.z, transformed_corner.z)};
		global_max = {std::max(global_max.x, transformed_corner.x),
					  std::max(global_max.y, transformed_corner.y),
					  std::max(global_max.z, transformed_corner.z)};//get min max
    }
    return AABB{global_min, global_max };
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

	// extrait les coordonnes UV des sommets
	double2 const &uv0 = tex_coords[tri[0].ti];
	double2 const &uv1 = tex_coords[tri[1].ti];
	double2 const &uv2 = tex_coords[tri[2].ti];

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
		//coord barycentrique pour UV
		double air = length(cross(p1 - p0, p2 - p0)); // air du triangle
		double a0 = length(cross(intersection - p1, intersection - p2)) / air;
		double a1 = length(cross(intersection - p2, intersection - p0)) / air;
		double a2 = length(cross(intersection - p0, intersection - p1)) / air;

		hit->uv = a0 * uv0 + a1 * uv1 + a2 * uv2;
      

        
		return true;
	}

	return false;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
	//put all the points in global first
	std::vector<linalg::aliases::double3> global_positions;
	for(const auto& point : positions) {

	 double4 global_point = mul( transform, {point.x, point.y, point.z, 1.0});
	 double3 global_point3 = {global_point.x, global_point.y, global_point.z};
	 global_positions.push_back(global_point3);
	}
	return construct_aabb(global_positions);
}