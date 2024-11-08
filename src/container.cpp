#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
#include <stack> 

bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    bool isHit = false;
    std::stack<BVHNode*> toVisit;
    toVisit.push(root);

    while (!toVisit.empty()) {
        BVHNode* node = toVisit.top(); // Get the top node
        toVisit.pop(); // Remove the top node

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
                // Push children to the stack, maintain depth-first traversal
                if (node->right != nullptr) toVisit.push(node->right);
                if (node->left != nullptr) toVisit.push(node->left);
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
