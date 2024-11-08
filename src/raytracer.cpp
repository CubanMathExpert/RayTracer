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

	double apertureRadius = tan(deg2rad(scene.camera.defocus_angle)) * scene.camera.focus_distance;

    // Itère sur tous les pixels de l'image.
    for(int y = 0; y < scene.resolution[1]; y++) {
		if (y % 40){
			std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
		}

        for(int x = 0; x < scene.resolution[0]; x++) {

			int avg_z_depth = 0; 
			double3 avg_ray_color{0,0,0};
			
			
			for(int iray = 0; iray < scene.samples_per_pixel; iray++) {
				// Génère le rayon approprié pour ce pixel.
				Ray ray;
				// Initialise la profondeur de récursivité du rayon.
				int ray_depth = 16;
				// Initialize la couleur du rayon
				double3 ray_color{0,0,0};
				double z_depth = scene.camera.z_far;

				double3 jitter{random_in_unit_disk() * scene.jitter_radius, 0};
				double3 pxy = pixel00 + (x + jitter.x) * pixel_delta_u * vpRight - (y + jitter.y) * pixel_delta_v * vpUp;
				auto rayOrigin = scene.camera.position;
				auto rayDirection = normalize(pxy - scene.camera.position);

				// depth of field
				if (scene.camera.defocus_angle != 0) {
					double3 focalPoint = rayOrigin + scene.camera.focus_distance * rayDirection;
					rayOrigin += jitter * apertureRadius;
					rayDirection = normalize(focalPoint - rayOrigin);
				}

				ray = Ray(rayOrigin, rayDirection);

				trace(scene, ray, ray_depth, &ray_color, &z_depth);

				avg_ray_color += ray_color;
				avg_z_depth += z_depth;

			}

			avg_z_depth = avg_z_depth / scene.samples_per_pixel;
			avg_ray_color = avg_ray_color / scene.samples_per_pixel;


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
		auto* resourceManager = ResourceManager::Instance();
		Material& material = resourceManager->materials[hit.key_material];
		double3 local_color = shade(scene, hit);
		
		// Réflexion
		double3 reflection_color = {0, 0, 0};
        if (material.k_reflection > 0.0) {
            double3 reflection_direction = ray.direction - 2 * dot(ray.direction, hit.normal) * hit.normal;
            Ray reflection_ray(hit.position + EPSILON * hit.normal, reflection_direction);

            double3 reflected_color = {0, 0, 0};
            double reflection_depth = *out_z_depth;
            trace(scene, reflection_ray, ray_depth - 1, &reflected_color, &reflection_depth);
			
            reflection_color = reflected_color * material.k_reflection;
        }
		// Refraction
		// reference : https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
		//t = eta * i + (eta * cos(theta_i) - sqrt(1 - sin^2(theta_t))) * n
		double3 refraction_color = {0, 0, 0};
		if (material.k_refraction > 0.0) {
    		bool entering = dot(ray.direction, hit.normal) < 0;
			//seulement de l'air à une autre refraction je considere aussi le cas a l'interieur de l'objet
    		double eta = entering ? (1.0 / material.refractive_index) : material.refractive_index;
    		double3 refraction_normal = entering ? hit.normal : -hit.normal;
    		double cos_theta = -dot(refraction_normal, ray.direction);
    		double sin2_theta = eta * eta * (1.0 - cos_theta * cos_theta);

    		if (sin2_theta <= 1.0) {  // No total internal reflection
        		double3 t = eta * ray.direction + (eta * cos_theta - sqrt(1.0 - sin2_theta)) * refraction_normal;
        		Ray refraction_ray(hit.position - EPSILON * refraction_normal, t);

        		// Trace the refraction ray
        		double3 refracted_color = {0, 0, 0};
        		double refraction_depth = *out_z_depth;
        		trace(scene, refraction_ray, ray_depth - 1, &refracted_color, &refraction_depth);

        		refraction_color = refracted_color * material.k_refraction;
    		}
		}
		*out_color = local_color + reflection_color + refraction_color;
		*out_z_depth = hit.depth;
	} 
	else {
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
	Material& material = ResourceManager::Instance()->materials[hit.key_material];
    double3 texture_color;
	double3 color;
    // Check if the material has a texture
    if (material.texture_albedo.width() > 0 && material.texture_albedo.height() > 0) {
        // get info 
        int width = material.texture_albedo.width();
        int height = material.texture_albedo.height();
        double u = hit.uv.x;
        double v = hit.uv.y;

        // Convert UV coordinates to texture coordinates
        int x = static_cast<int>(u * width);
        int y = static_cast<int>(v * height);
        // clamp the values
        x = std::clamp(x, 0, width - 1);
        y = std::clamp(y, 0, height - 1);
		//convert to double3 the rgb_t 
		rgb_t color_pixel = material.texture_albedo.get_pixel(x, y);
        texture_color = double3(color_pixel.red / 255.0, color_pixel.green / 255.0, color_pixel.blue / 255.0);

		color = texture_color * material.k_ambient * scene.ambient_light;
    }else{
		color = material.color_albedo * material.k_ambient * scene.ambient_light;
	}

	// hit information
	double3 intersection_point = hit.position;
	double3 normal = hit.normal;
	double3 view_dir = normalize(scene.camera.position - intersection_point);
	
	// ALL OF THE LIGHTS
	for (const auto& light : scene.lights)
	{
		// Direction et distance vers la lumière
		double3 light_direction = normalize(light.position - intersection_point);
		double light_distance = length(light.position - intersection_point);
		double light_distance_inverse = 1.0 / light_distance;

		int unoccludedCount = 0;

		int echantillon = 16;
		for (int i = 0; i < echantillon; i++)
		{
			double3 jitter = double3{random_in_unit_disk(), 0} * light.radius;
			double3 jittered_light_position = light.position + jitter;
			double3 jittered_light_direction = normalize(jittered_light_position - intersection_point);
			double jittered_light_distance = length(jittered_light_position - intersection_point);

			Ray shadow_ray(intersection_point, jittered_light_direction);
			Intersection shadow_hit;

			//all of the shadows 
			if (!scene.container->intersect(shadow_ray, EPSILON, jittered_light_distance, &shadow_hit))
			{
				unoccludedCount++;
			}
		}
		
		// Calcul de l'intensité de la lumière seulement si le point n'est pas dans l'ombre
		if (unoccludedCount > 0)
		{
			double shadow_factor = static_cast<double>(unoccludedCount) / echantillon;

			// Diffuse
			double projection = std::max(dot(normal, light_direction), 0.0);
			double3 diffuse = material.k_diffuse * projection * material.color_albedo;

			// Spéculaire avec le modèle de Blinn
			double3 H = normalize(view_dir + light_direction);
			double spec_intensity = pow(std::max(dot(normal, H), 0.0), material.shininess);
			double3 specular = material.k_specular * spec_intensity *
								(material.metallic * material.color_albedo + (1 - material.metallic));

			// Ajout des contributions diffuse et spéculaire avec l'atténuation
			color += (diffuse + specular) * pow(light_distance_inverse, 2) * light.emission * shadow_factor;
		}
	}
	return color;
}
