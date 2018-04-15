#pragma once

#include <vector>
#include <map>
#include <random>
#include <cstdint>
#include <memory>

#include <glm/glm.hpp>

#include "utils.h"
#include "Object3D.h"
#include "../external_libraries/common_include/kdtree++/kdtree.hpp"

class Scene
{
public:
	Scene();
	~Scene();
	
	enum class RenderMode
    {
	    PHOTON_MAPPING,
        CAUSTICS,
        WHITTED_SPECULAR,
        MONTE_CARLO,
	};
	
	SpectralDistribution traceRay(Ray r, RenderMode render_mode, int iteration = 0);

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

    bool intersect(IntersectionData* id, Ray r);
    bool intersectLamp(LightSourceIntersectionData* light_id, Ray r);

private:
	std::unique_ptr<std::mt19937> gen_;
	std::unique_ptr<std::uniform_real_distribution<float>> dis_;

	std::vector<Object3D*> objects_;
	std::vector<LightSource*> lamps_;
	std::map<std::string, Material*> materials_;

	KDTree::KDTree<3, KDTreeNode> photon_map_;
};