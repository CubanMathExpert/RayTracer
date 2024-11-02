#include "raytracer.h"
#include <iostream>

// Overload the << operator for double3 type
std::ostream& operator<<(std::ostream& os, const double3& vec) {
	os << "{" << vec.x << ", " << vec.y << ", " << vec.z << "}";
	return os;
}

void Raytracer::render(const Scene& scene, Frame* output)
{       
    // Crée le z_buffer.
    double *z_buffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        z_buffer[i] = scene.camera.z_far; //Anciennement DBL_MAX. À remplacer avec la valeur de scene.camera.z_far
    }

	// @@@@@@ VOTRE CODE ICI
	// Calculez les paramètres de la caméra pour les rayons.
	double3 cameraDirection = normalize(scene.camera.center - scene.camera.position); // camera direction vector
	double3 vpRight = normalize(cross(cameraDirection, scene.camera.up)); // right vector for viewport
	double3 vpUp = -normalize(cross(vpRight, cameraDirection)); // up vector for viewport
	double vpDistance = scene.camera.z_near; // viewport distance from camera

	// view volume
	// the view is orthogonal which makes the view volume a box, easier to compute (t,b,l,r)
	double topPlane = vpDistance * tan(deg2rad(scene.camera.fovy / 2));
	double rightPlane = topPlane * scene.camera.aspect; // why multiply with aspect?
	double leftPlane = -rightPlane;
	double bottomPlane = -topPlane;
	
	// distance between pixels
	double pixel_delta_u = (2 * fabs(rightPlane)) / scene.resolution[0];
	double pixel_delta_v = (2 * topPlane) / scene.resolution[1];

	// starting point of the viewport (upper left)
	double3 zeroPlane = (scene.camera.position + scene.camera.z_near * cameraDirection) + (leftPlane * vpRight) + (topPlane * vpUp);
	double3 pixel00 = zeroPlane + (0.5 * pixel_delta_u * vpRight) - (0.5 * pixel_delta_v * vpUp);

    // Itère sur tous les pixels de l'image.
    for(int y = 0; y < scene.resolution[1]; y++) {
		if (y % 40){
			std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
		}

        for(int x = 0; x < scene.resolution[0]; x++) {

			double avg_z_depth = 0; // changed this to double
			double3 avg_ray_color{0,0,0};
			//std::cout << "x: " << x << " y: " << y << std::endl;
			
			for(int iray = 0; iray < scene.samples_per_pixel; iray++) {
				// Génère le rayon approprié pour ce pixel.
				Ray ray;
				// Initialise la profondeur de récursivité du rayon.
				int ray_depth = 0;
				// Initialize la couleur du rayon
				double3 ray_color{0,0,0};
				double z_depth = scene.camera.z_far;

				// @@@@@@ VOTRE CODE ICI
				// Mettez en place le rayon primaire en utilisant les paramètres de la caméra.
				// Lancez le rayon de manière uniformément aléatoire à l'intérieur du pixel dans la zone délimité par jitter_radius. 
				//Faites la moyenne des différentes couleurs obtenues suite à la récursion.
				//double3 rayOrigin = camOrth.minPosition + uVec * x_shift * x + vVec * y_shift * y;
				//auto pixelCenter = scene.camera.position + (x * pixel_delta_u * vpRight) - (y * pixel_delta_v * vpUp);
				double3 jitter{random_in_unit_disk() * scene.jitter_radius, 0};
				double3 pxy = pixel00 + (x + jitter.x) * pixel_delta_u * vpRight - (y + jitter.y) * pixel_delta_v * vpUp;
				auto rayOrigin = scene.camera.position;
				auto rayDirection = normalize(pxy - scene.camera.position);
				ray = Ray(rayOrigin, rayDirection);

				// trace ray and set values to ray_color and ray_depth
				trace(scene, ray, ray_depth, &ray_color, &z_depth);

				avg_ray_color += ray_color;
				avg_z_depth += z_depth;

			}

			avg_z_depth = avg_z_depth / scene.samples_per_pixel;
			avg_ray_color = avg_ray_color / scene.samples_per_pixel;

			//std::cout << avg_ray_color << std::endl;
			//std::cout << avg_z_depth << std::endl;

			// Test de profondeur
			if(avg_z_depth >= scene.camera.z_near && avg_z_depth <= scene.camera.z_far && 
				avg_z_depth < z_buffer[x + y*scene.resolution[0]]) {
				z_buffer[x + y*scene.resolution[0]] = avg_z_depth;

				// Met à jour la couleur de l'image (et sa profondeur)
				output->set_color_pixel(x, y, avg_ray_color);
				output->set_depth_pixel(x, y, (avg_z_depth - scene.camera.z_near) / 
										(scene.camera.z_far-scene.camera.z_near));
			}
        }
    }

    delete[] z_buffer;
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		- Détermine si le rayon intersecte la géométrie.
//      	- Calculer la contribution associée à la réflexion.
//			- Calculer la contribution associée à la réfraction.
//			- Mettre à jour la couleur avec le shading + 
//			  Ajouter réflexion selon material.reflection +
//			  Ajouter réfraction selon material.refraction 
//            pour la couleur de sortie.
//          - Mettre à jour la nouvelle profondeure.
void Raytracer::trace(const Scene& scene,
					  Ray ray, int ray_depth,
					  double3* out_color, double* out_z_depth)
{
	Intersection hit;
	// Fait appel à l'un des containers spécifiées.
	if(scene.container->intersect(ray,EPSILON,*out_z_depth,&hit)) {		
		Material& material = ResourceManager::Instance()->materials[hit.key_material];

		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réflection d'un rayon de manière récursive.
		
		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réfraction d'un rayon de manière récursive.
		// 
		// Assumez que l'extérieur/l'air a un indice de réfraction de 1.
		//
		// Toutes les géométries sont des surfaces et non pas de volumes.



		*out_color = material.color_albedo;
		*out_z_depth = hit.depth;
	} 
	else 
	{
		*out_color = double3{0,0,0};
	}
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		* Calculer la contribution des lumières dans la scène.
//			- Itérer sur toutes les lumières.
//				- Inclure la contribution spéculaire selon le modèle de Blinn en incluant la composante métallique.
//	          	- Inclure la contribution diffuse. (Faites attention au produit scalare. >= 0)
//   	  	- Inclure la contribution ambiante
//      * Calculer si le point est dans l'ombre
//			- Itérer sur tous les objets et détecter si le rayon entre l'intersection et la lumière est occludé.
//				- Ne pas considérer les points plus loins que la lumière.
//			- Par la suite, intégrer la pénombre dans votre calcul
//		* Déterminer la couleur du point d'intersection.
//        	- Si texture est présente, prende la couleur à la coordonnées uv
//			- Si aucune texture, prendre la couleur associé au matériel.

double3 Raytracer::shade(const Scene& scene, Intersection hit)
{
	// Material& material = ResourceManager::Instance()->materials[hit.key_material]; lorsque vous serez rendu à la partie texture.
	return double3{0,0,0};
}
