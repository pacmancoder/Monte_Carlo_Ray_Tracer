#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "utils.h"

class Camera
{
public:
	glm::vec3 eye_;
	glm::vec3 center_;
	glm::vec3 up_;
	float fov_;
	glm::mat4 VP_inv;

	const uint16_t width_;
	const uint16_t height_;

	Camera(
		glm::vec3 eye,
		glm::vec3 center,
		glm::vec3 up,
		float fov,
		uint16_t baseWidth,
		uint16_t baseHeight);

	Ray castRay(
		uint16_t frameWidth,
		uint16_t frameHeight,
		uint16_t pixel_x,
		uint16_t pixel_y,
		float parameter_x,
		float parameter_y) const;

public:
	Camera& operator=(const Camera&) = delete;
};