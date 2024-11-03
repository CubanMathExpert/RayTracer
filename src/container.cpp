#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    bool isHit = false;
    std::vector<BVHNode*> toVisit;
    toVisit.push_back(root);

    while (!toVisit.empty()) {
        BVHNode* node = toVisit.back();
        toVisit.pop_back();

        if (node->aabb.intersect(ray, t_min, t_max)) {
            if (node->left == nullptr && node->right == nullptr) {
                Intersection closestIntersection;
                if (objects[node->idx]->intersect(ray, t_min, t_max, &closestIntersection)) {
                    if (closestIntersection.depth < t_max) {
                        t_max = closestIntersection.depth;
                        *hit = closestIntersection;
                        isHit = true;
                    }
                }
            } else {
               if(node->left != nullptr) toVisit.push_back(node->left);
               if(node->right != nullptr) toVisit.push_back(node->right);
            }
        }
    }
    return isHit;
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
        AABB objAABB = obj->compute_aabb();

        // Intersection the AABB
        if (!objAABB.intersect(ray, t_min, t_max)) {
            continue;  // Skip if no intersection 
        }

        // Check actual object
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
