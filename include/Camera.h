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

	const int width_;
	const int height_;

	Camera(
		glm::vec3 eye,
		glm::vec3 center,
		glm::vec3 up,
		float fov,
		int width,
		int height);

	Ray castRay(
		int pixel_x, // [0, WIDTH_ - 1]
		int pixel_y, // [0, HEIGHT_ - 1]
		float parameter_x, // [-0.5, 0.5]
		float parameter_y); // [-0.5, 0.5]
};