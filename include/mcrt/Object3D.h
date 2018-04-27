#pragma once

#include "OctTreeAABB.h"
#include "utils.h"

#include <vector>
#include <array>
#include <memory>

#include <glm/glm.hpp>

namespace Mcrt
{
	class Object3D
	{
	public:
		explicit Object3D(Material *material);
		virtual ~Object3D() = default;

		virtual bool intersect(IntersectionData *id, Ray r) const = 0;
		Material material() const;

	private:
		const Material *material_;
	};

	class Mesh : public Object3D
	{
	public:
		Mesh(
				Material *material,
				glm::mat4 transform,
				std::vector<glm::vec3> &&positions,
				std::vector<glm::vec3> &&normals,
				std::vector<glm::vec2> &&uvs);

		bool intersect(IntersectionData *id, Ray r) const override;

		glm::vec3 getMinPosition() const;
		glm::vec3 getMaxPosition() const;
		glm::mat4 getTransform() const;
		size_t getNumberOfTriangles() const;

	private:
		std::vector<glm::vec3> positions_;
		std::vector<glm::vec2> uvs_;
		std::vector<glm::vec3> normals_;
		std::vector<size_t> indices_;

		glm::mat4 transform_;
		std::unique_ptr<OctTreeAABB> ot_aabb_;

		friend class OctNodeAABB;
	};

	class Sphere : public Object3D
	{
	public:
		Sphere(Material *material, glm::vec3 position, float radius);

		bool intersect(IntersectionData *id, Ray r) const override;
		glm::vec3 getPointOnSurface(float u, float v) const;

	private:
		const glm::vec3 position_;
		const float radius_;
	};

// P0, P1, and P2 defines a paralellogram which is the plane
	class Plane : public Object3D
	{
	public:
		Plane(Material *material, glm::vec3 p0, glm::vec3 p1, glm::vec3 p2);

		bool intersect(IntersectionData *id, Ray r) const override;

		glm::vec3 getPointOnSurface(float u, float v) const;
		float getArea() const;
		glm::vec3 getNormal() const;
		glm::vec3 getFirstTangent() const;

	private:
		const glm::vec3 p0_, p1_, p2_, normal_;
		const float area_;
	};

	class LightSource
	{
	private:
		const Plane emitter_;
	public:
		LightSource(
				glm::vec3 p0,
				glm::vec3 p1,
				glm::vec3 p2,
				float flux, // Gets multiplied with color for total flux [Watts]
				SpectralDistribution color);

		bool intersect(LightSourceIntersectionData *light_id, Ray r);
		glm::vec3 getPointOnSurface(float u, float v) const;
		float getArea() const;
		glm::vec3 getNormal() const;

		Ray shootLightRay();

		const SpectralDistribution radiosity_; // [Watts/m^2]
	};
}