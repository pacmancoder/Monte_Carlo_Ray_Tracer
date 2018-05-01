#include <mcrt/Scene.h>

#include <string>
#include <random>

#include <dcdr/logging/Logger.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/rotate_vector.hpp>


using namespace Mcrt;
using namespace Dcdr::Logging;

// --- Scene class functions --- //

Scene::Scene(uint16_t frameWidth, uint16_t frameHeight) :
    frameWidth_(frameWidth),
    frameHeight_(frameHeight),
    gen_(std::make_unique<std::mt19937>(std::random_device()())),
    dis_(std::make_unique<std::uniform_real_distribution<float>>(0, 1)),
	cameras_(),
	objects_(),
	lightSources_(),
	materials_(),
	photonMapBuilt_(false) {}

bool Scene::intersect(IntersectionData* id, Ray r)
{
	IntersectionData id_smallest_t;
	id_smallest_t.t = 100000; // Ugly solution

	Object3D* intersecting_object = nullptr;

	for (auto& object : objects_)
	{
		IntersectionData id_local;
		if (object.second->intersect(&id_local,r) && id_local.t < id_smallest_t.t)
		{
			id_smallest_t = id_local;
			intersecting_object = object.second.get();
		}
	}

	if (intersecting_object)
	{
		*id = id_smallest_t;
		return true;
	}

	return false;
}

bool Scene::intersectLamp(LightSourceIntersectionData* light_id, Ray r)
{
	LightSourceIntersectionData lamp_id_smallest_t = {};
	lamp_id_smallest_t.t = 100000; // Ugly solution

	LightSource* intersecting_lamp = nullptr;
	for (const auto& lightSource : lightSources_)
	{
		LightSourceIntersectionData id_local  = {};
		if (lightSource.second->intersect(&id_local,r) && id_local.t < lamp_id_smallest_t.t)
		{
			lamp_id_smallest_t = id_local;
			intersecting_lamp = lightSource.second.get();
		}
	}
	if (intersecting_lamp)
	{
		IntersectionData id_smallest_t;
		id_smallest_t.t = 100000; // Ugly solution

		Object3D* intersecting_object = nullptr;
		for (const auto& object : objects_)
		{
			IntersectionData id_local;
			if (object.second->intersect(&id_local,r) && id_local.t < id_smallest_t.t)
			{
				id_smallest_t = id_local;
				intersecting_object = object.second.get();
			}
		}
		if (intersecting_object && id_smallest_t.t < lamp_id_smallest_t.t)
		{
			return false;
		}
		else
		{
			*light_id = lamp_id_smallest_t;
			return true;	
		}
	}
	return false;
}

SpectralDistribution Scene::traceDiffuseRay(
	Ray r,
	RenderMode render_mode,
	IntersectionData id,
	int iteration)
{

	r.has_intersected = true;
	// Start by adding the local illumination part (shadow rays)
	SpectralDistribution total_diffuse = traceLocalDiffuseRay(r, render_mode, id);
    // Add the indirect illumination part (Monte Carlo sampling)
    total_diffuse = total_diffuse + traceIndirectDiffuseRay(r, render_mode, id, iteration);
    return total_diffuse;

    return SpectralDistribution();
}

SpectralDistribution Scene::traceLocalDiffuseRay(
	Ray r,
	RenderMode /*render_mode*/,
	IntersectionData id)
{
	SpectralDistribution L_local;
	// Cast shadow rays
	// We divide up_ the area light source in to n_samples area parts.
	// Used to define the solid angle
	static const int n_samples = 1;
	for (const auto& lightSource : lightSources_)
	{
		for (int j = 0; j < n_samples; ++j)
		{
			Ray shadow_ray = r;
			glm::vec3 differance = lightSource.second->getPointOnSurface((*dis_)(*gen_),(*dis_)(*gen_)) - shadow_ray.origin;
			shadow_ray.direction = glm::normalize(differance);

			SpectralDistribution brdf; // = id.material.color_diffuse / (2 * M_PI); // Dependent on inclination and azimuth
			float cos_theta = glm::dot(shadow_ray.direction, id.normal);

			LightSourceIntersectionData shadow_ray_id;

			if (id.material.diffuse_roughness)
			{
				brdf = evaluateOrenNayarBRDF(
					-r.direction,
					shadow_ray.direction,
					id.normal,
					id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance),
					id.material.diffuse_roughness);
			}
			else
			{
				brdf = evaluateLambertianBRDF(
						-r.direction,
						shadow_ray.direction,
						id.normal,
						id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance));
			}

			if(intersectLamp(&shadow_ray_id, shadow_ray))
			{
				float cos_light_angle = glm::dot(shadow_ray_id.normal, -shadow_ray.direction);
				float light_solid_angle = shadow_ray_id.area / n_samples * glm::clamp(cos_light_angle, 0.0f, 1.0f) / glm::pow(glm::length(differance), 2) / (M_PI * 2);

				L_local +=
					brdf *
					shadow_ray_id.radiosity *
					cos_theta *
					light_solid_angle;
			}
		}
	}
	return L_local;
}

SpectralDistribution Scene::traceIndirectDiffuseRay(
	Ray r,
	RenderMode render_mode,
	IntersectionData id,
	int iteration)
{
	SpectralDistribution L_indirect;
	static const int n_samples = 1;
	for (int i = 0; i < n_samples; ++i)
	{
		// helper is just a random vector and can not possibly be
		// a zero vector since id.normal is normalized
		glm::vec3 helper = id.normal + glm::vec3(1, 1, 1);
		glm::vec3 tangent = glm::normalize(glm::cross(id.normal, helper));

		// rand1 is a random number from the cosine estimator
		float rand1 = (*dis_)(*gen_);
		float rand2 = (*dis_)(*gen_);

		// Uniform distribution over a hemisphere
		auto inclination = static_cast<float>(std::acos(std::sqrt(rand1)));//glm::acos(1 - rand1);//glm::acos(1 -  2 * (*dis_)(*gen_));
		auto azimuth = static_cast<float>(2 * M_PI * rand2);

		// Change the actual vector
        glm::vec3 random_direction = id.normal;
        random_direction = glm::normalize(glm::rotate(
                random_direction,
                inclination,
                tangent));

        random_direction = glm::normalize(glm::rotate(
                random_direction,
                azimuth,
                id.normal));

		SpectralDistribution brdf;
		if (id.material.diffuse_roughness)
		{
			brdf = evaluateOrenNayarBRDF(
				-r.direction,
				random_direction,
				id.normal,
				id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance),
				id.material.diffuse_roughness);
		}
		else
		{
			brdf = evaluateLambertianBRDF(
				-r.direction,
				random_direction,
				id.normal,
				id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance));
		}

		r.direction = random_direction;
		r.radiance *= brdf * M_PI; // Importance, M_PI is because of the importance sampling
		L_indirect += traceRay(r, render_mode, iteration + 1) * M_PI * brdf;
	}
	return L_indirect / n_samples;
}

SpectralDistribution Scene::traceSpecularRay(
	Ray r,
	RenderMode render_mode,
	IntersectionData id,
	int iteration)
{
	r.has_intersected = true;
	SpectralDistribution specular = SpectralDistribution();

	r.direction = glm::reflect(r.direction, id.normal);
	SpectralDistribution brdf = evaluatePerfectBRDF(id.material.color_specular * id.material.reflectance * id.material.specular_reflectance);
	r.radiance *= brdf;
	// Recursively trace the reflected ray
	specular += traceRay(r, render_mode, iteration + 1) * brdf;
	return specular;
}

SpectralDistribution Scene::traceRefractedRay(
	Ray r,
	RenderMode render_mode,
	IntersectionData id,
	int iteration,
	glm::vec3 offset,
	bool inside)
{
	Ray recursive_ray = r;
	recursive_ray.has_intersected = true;

	glm::vec3 normal = inside ? -id.normal : id.normal;
	glm::vec3 perfect_refraction = glm::refract(
		r.direction,
		normal,
		r.material.refraction_index / id.material.refraction_index);
	glm::vec3 perfect_reflection = glm::reflect(r.direction, id.normal);
	if (perfect_refraction != glm::vec3(0))
	{ // Refraction and reflection
		// Schlicks approximation to Fresnels equations.
		float n1 = r.material.refraction_index;
		float n2 = id.material.refraction_index;
		auto R_0 = static_cast<float>(std::pow((n1 - n2) / (n1 + n2), 2));
		auto R = static_cast<float>(R_0 + (1 - R_0) * pow(1 - glm::dot(normal, -r.direction), 5));

		Ray recursive_ray_reflected = recursive_ray;
		Ray recursive_ray_refracted = recursive_ray;

		if (inside)
			offset = -offset;
		
		// Reflected ray
		// Change the material the ray is travelling in
		recursive_ray_reflected.material = Material::air();
		recursive_ray_reflected.origin = r.origin + id.t * r.direction + offset;
		// Refracted ray
		// Change the material the ray is travelling in
		recursive_ray_refracted.material = id.material;
		recursive_ray_refracted.origin = r.origin + id.t * r.direction - offset;
		
		SpectralDistribution to_return;
		recursive_ray_reflected.direction = perfect_reflection;
		recursive_ray_refracted.direction = perfect_refraction;

		SpectralDistribution brdf_specular = evaluatePerfectBRDF(id.material.color_specular * id.material.reflectance * id.material.specular_reflectance * R);
		SpectralDistribution brdf_refractive = evaluatePerfectBRDF(id.material.color_diffuse * id.material.reflectance * id.material.specular_reflectance * (1 - R));


		recursive_ray_reflected.radiance *= brdf_specular;
		recursive_ray_refracted.radiance *= brdf_refractive;

		// Recursively trace the refracted rays
		SpectralDistribution reflected_part = traceRay(recursive_ray_reflected, render_mode, iteration + 1) * brdf_specular;
		SpectralDistribution refracted_part = traceRay(recursive_ray_refracted, render_mode, iteration + 1) * brdf_refractive;
		return reflected_part + refracted_part;
	}
	else
	{ // Brewster angle reached, complete specular reflection
		if (inside)
			recursive_ray.origin = r.origin + id.t * r.direction - offset;
		else
			recursive_ray.origin = r.origin + id.t * r.direction + offset;

		SpectralDistribution brdf_specular = evaluatePerfectBRDF(id.material.color_specular * id.material.reflectance * id.material.specular_reflectance);
		recursive_ray.direction = perfect_reflection;
		recursive_ray.radiance *= brdf_specular;
		// Recursively trace the reflected ray
		return traceRay(recursive_ray, render_mode, iteration + 1) * brdf_specular;
	}
}

SpectralDistribution Scene::traceRay(Ray r, RenderMode render_mode, int iteration)
{
	const int MAX_ITERATIONS = 5;
	IntersectionData id;
	LightSourceIntersectionData lamp_id;
	if (intersectLamp(&lamp_id, r)) // Ray hit light source
    {
        switch (render_mode)
        {
            case RenderMode::WHITTED_SPECULAR:
                return lamp_id.radiosity / (M_PI * 2);
            default:
                return SpectralDistribution();
        }
    }
	else if (intersect(&id, r))
	{
	    // Ray hit another object
		// Russian roulette
		float random = (*dis_)(*gen_);
		//float non_termination_probability = glm::max((1 - float(iteration) / 10), 0.5f);
		auto non_termination_probability = static_cast<float>(iteration == 0 ? 1.0 : 0.8);
		if (random > non_termination_probability || iteration > MAX_ITERATIONS)
			return SpectralDistribution();

		// To make sure it does not intersect with itself again
		glm::vec3 offset = id.normal * 0.00001f;
		bool inside = false;
		if (glm::dot(id.normal, r.direction) > 0) // The ray is inside an object
        {
            inside = true;
        }
		
		float transmissivity = id.material.transmissivity;
		float specularity = id.material.specular_reflectance;

		SpectralDistribution total;
		if (1 - transmissivity)
		{ // Completely or partly reflected
			Ray recursive_ray = r;
			// New position same in both cases, can be computed once, outside
			// of trace functions
			recursive_ray.origin = r.origin + id.t * r.direction +
				(inside ? -offset : offset);
			SpectralDistribution specular_part =
				specularity ?
					traceSpecularRay(
						recursive_ray,
						render_mode,
						id,
						iteration) :
					SpectralDistribution();
			SpectralDistribution diffuse_part;
			switch (render_mode)
			{
				case RenderMode::PHOTON_MAPPING:
				{
					if (r.has_intersected)
					{
						Photon p;
						p.position = recursive_ray.origin;
						p.direction_in = -r.direction;

						auto photon_area = static_cast<float>(Photon::radius * Photon::radius * M_PI);
						// The projected area should be photon area times cos theta,
						// This is avoided both here and later to avoid numerical problem
						// when dividing with small numbers.
						float projected_area = photon_area;// * glm::dot(p.direction_in, id.normal);
						auto solid_angle = static_cast<float>(M_PI);
			
						p.delta_flux = recursive_ray.radiance / non_termination_probability * projected_area * solid_angle;

						KDTreeNode to_insert;
						to_insert.p = p;

						photon_map_.insert(to_insert);
					}
					break;
				}
				case RenderMode::CAUSTICS:
				{
					KDTreeNode ref_node;
					ref_node.p.position = r.origin + r.direction * id.t + offset;

					std::vector<KDTreeNode> closest_photons;
					photon_map_.find_within_range(ref_node,Photon::radius,std::back_insert_iterator<std::vector<KDTreeNode> >(closest_photons));
					
					SpectralDistribution photon_radiance;
					for (size_t i = 0; i < closest_photons.size(); ++i)
					{
						SpectralDistribution brdf;
						if (id.material.diffuse_roughness)
						{
							brdf = evaluateOrenNayarBRDF(
								-r.direction,
								closest_photons[i].p.direction_in,
								id.normal,
								id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance),
								id.material.diffuse_roughness);
						}
						else
						{
							brdf = evaluateLambertianBRDF(
								-r.direction,
								closest_photons[i].p.direction_in,
								id.normal,
								id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance));
						}

						float distance = glm::length(closest_photons[i].p.position - ref_node.p.position);
						// The area of the photon if its inclination angle
						// is 90 degrees and the surface is flat.
						//float cos_theta = glm::max(glm::dot(closest_photons[i].p.direction_in, id.normal), 0.0f);
						auto photon_area = static_cast<float>(Photon::radius * Photon::radius * M_PI);
						float projected_area = photon_area;// * cos_theta;
						photon_radiance +=
							// flux / area / steradian = radiance
							closest_photons[i].p.delta_flux *
							(glm::length(distance) < Photon::radius ? 1 : 0)
							/ (projected_area * 2 * M_PI)
							//
							* brdf // The brdf is part of the integral of the rendering equation
							* (2 * M_PI); // Integration over the whole hemisphere get us back to radiance
					}

					diffuse_part = photon_radiance;
					break;
				}
				case RenderMode::WHITTED_SPECULAR:
				{
					diffuse_part = SpectralDistribution();
					break;
				}
				case RenderMode::MONTE_CARLO:
				{
					diffuse_part =
						(1 - specularity) ?
							traceDiffuseRay(
								recursive_ray,
								render_mode,
								id,
								iteration) :
							SpectralDistribution();

					break;
				}
			}

			total += (specular_part + diffuse_part) * (1 - transmissivity);
		}
		if (transmissivity)
		{ // Completely or partly transmissive
			SpectralDistribution transmitted_part =
				traceRefractedRay(r, render_mode, id, iteration, offset, inside);
			total += transmitted_part * transmissivity;
		}
		return total / non_termination_probability;
	}
	return SpectralDistribution();
}

void Scene::buildPhotonMap(int n_photons)
{

	if (photonMapBuilt_ || lightSources_.empty())
	{
		return;
	}

	SpectralDistribution total_flux = SpectralDistribution();
	float total_flux_norm = 0;
	for (const auto& lightSource : lightSources_)
	{
		total_flux_norm += lightSource.second->radiosity_.norm() * lightSource.second->getArea();
		total_flux += lightSource.second->radiosity_ * lightSource.second->getArea();
	}
	for (int k = 0; k < 100; ++k)
	{
		//#pragma omp parallel for
		for (int i = 0; i < n_photons / 100; ++i)
		{
			// Pick a light source. Bigger flux => Bigger chance to be picked.
			float accumulating_chance = 0;
			float random = (*dis_)(*gen_);

			auto picked_light_source = lightSources_.cbegin();

			for (auto lightSource = lightSources_.cbegin(); lightSource != lightSources_.cend(); ++lightSource)
			{
				float interval =
						lightSource->second->radiosity_.norm() *
						lightSource->second->getArea() /
						total_flux_norm;
				if (random > accumulating_chance && random < accumulating_chance + interval)
				{
					// This lamp got picked
					picked_light_source = lightSource;
					break;
				}
				else
				{
					accumulating_chance += interval;
				}
			}

			Ray r = picked_light_source->second->shootLightRay();

			r.has_intersected = false;
			// Compute delta_flux based on the flux of the light source
			SpectralDistribution delta_flux = total_flux / n_photons;
			auto photon_area = static_cast<float>(Photon::radius * Photon::radius * M_PI);
			auto solid_angle = static_cast<float>(M_PI * 2);
			r.radiance = delta_flux / (photon_area * solid_angle);
			traceRay(r, RenderMode::PHOTON_MAPPING);
		}
	}
	photon_map_.optimize();

	photonMapBuilt_ = true;
}

size_t Scene::getNumberOfPhotons()
{
	return photon_map_.size();
}

void Scene::AddCamera(Scene::SceneEntityId id, const Scene::CameraPtr& camera)
{
    cameras_[id] = camera;
}

void Scene::AddObject(Scene::SceneEntityId id, const Scene::ObjectPtr& object)
{
    objects_[id] = object;
}

void Scene::AddLightSource(Scene::SceneEntityId id, const Scene::LightSourcePtr& lightSource)
{
    lightSources_[id] = lightSource;
}

void Scene::AddMaterial(Scene::SceneEntityId id, const Scene::MaterialPtr& material)
{
    materials_[id] = material;
}

const Camera& Scene::GetDefaultCamera() const
{
    return *cameras_.cbegin()->second;
}

uint16_t Scene::get_frame_width() const
{
	return frameWidth_;
}

uint16_t Scene::get_frame_height() const
{
	return frameHeight_;
}
