#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
	return false;
}


// @@@@@@ VOTRE CODE ICI
// - Parcourir tous les objets
// 		- Détecter l'intersection avec l'AABB
//			- Si intersection, détecter l'intersection avec la géométrie.
//				- Si intersection, mettre à jour les paramètres.
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool Naive::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
	bool isHit = false;
    Intersection closestIntersection;

    for (auto obj : objects) {
        // Get the AABB for the current object
        AABB objAABB = obj->compute_aabb();

        // First, check if the ray intersects the AABB
        if (!objAABB.intersect(ray, t_min, t_max)) {
            continue;  // Skip to the next object if no intersection with the AABB
        }

        // If there's an AABB intersection, check for an intersection with the actual object
        if (obj->intersect(ray, t_min, t_max, &closestIntersection)) {
            if (closestIntersection.depth < t_max) {
                t_max = closestIntersection.depth;
                *hit = closestIntersection;
                isHit = true;
            }
        }
    }
    return isHit;
}
