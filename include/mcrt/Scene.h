#pragma once

#include <vector>
#include <map>
#include <random>
#include <cstdint>
#include <memory>

#include <glm/glm.hpp>

#include <mcrt/utils.h>
#include <mcrt/Object3D.h>
#include <mcrt/Camera.h>
#include <external/kdtree++/kdtree.hpp>

namespace Mcrt
{
	class Scene
	{
	public:
		Scene(uint16_t frameWidth, uint16_t frameHeight);

		enum class RenderMode
		{
			PHOTON_MAPPING,
			CAUSTICS,
			WHITTED_SPECULAR,
			MONTE_CARLO,
		};

		using SceneEntityId = uint32_t;
		using ObjectPtr = std::shared_ptr<Object3D>;
		using LightSourcePtr = std::shared_ptr<LightSource>;
		using MaterialPtr = std::shared_ptr<Material>;
		using CameraPtr = std::shared_ptr<Camera>;

	public:
		void AddObject(SceneEntityId id, const ObjectPtr &object);
		void AddLightSource(SceneEntityId id, const LightSourcePtr &lightSource);
		void AddMaterial(SceneEntityId id, const MaterialPtr &material);
		void AddCamera(SceneEntityId id, const CameraPtr &camera);

		const Camera &GetDefaultCamera() const;

	public: // Tracing
		SpectralDistribution traceRay(Ray r, RenderMode render_mode, int iteration = 0);


	public: // Photon map
		void buildPhotonMap(int n_photons);
		size_t getNumberOfPhotons();

	private:
		// Normal path tracing for diffuse ray
		SpectralDistribution traceDiffuseRay(
				Ray r,
				RenderMode render_mode,
				IntersectionData id,
				int iteration);
		SpectralDistribution traceLocalDiffuseRay(
				Ray r,
				RenderMode render_mode,
				IntersectionData id);
		SpectralDistribution traceIndirectDiffuseRay(
				Ray r,
				RenderMode render_mode,
				IntersectionData id,
				int iteration);

		// Specular and refractive tracing
		SpectralDistribution traceSpecularRay(
				Ray r,
				RenderMode render_mode,
				IntersectionData id,
				int iteration);
		SpectralDistribution traceRefractedRay(
				Ray r,
				RenderMode render_mode,
				IntersectionData id,
				int iteration,
				glm::vec3 offset,
				bool inside);

		bool intersect(IntersectionData *id, Ray r);
		bool intersectLamp(LightSourceIntersectionData *light_id, Ray r);

	public:
		uint16_t get_frame_width() const;
		uint16_t get_frame_height() const;

	private:
		uint16_t frameWidth_;
		uint16_t frameHeight_;

		std::unique_ptr<std::mt19937> gen_;
		std::unique_ptr<std::uniform_real_distribution<float>> dis_;

		std::map<SceneEntityId, CameraPtr> cameras_;
		std::map<SceneEntityId, ObjectPtr> objects_;
		std::map<SceneEntityId, LightSourcePtr> lightSources_;
		std::map<SceneEntityId, MaterialPtr> materials_;

		KDTree::KDTree<3, KDTreeNode> photon_map_;
		bool photonMapBuilt_;
	};
}