#include "../include/Camera.h"

Camera::Camera (
	glm::vec3 eye,
	glm::vec3 center,
	glm::vec3 up,
	float fov,
	int width,
	int height) :

	eye_(eye),
	center_(center),
	up_(up),
	fov_(fov),
	width_(width),
	height_(height)
{
	// View and perspective matrices are used in the unProject() function
	glm::mat4 V = glm::lookAt(eye, center, up);
	float aspect = float(width_) / height_;
	glm::mat4 P = glm::perspective(fov, aspect, 0.1f, 100.0f);
	VP_inv = glm::inverse(V * P);
}

Ray Camera::castRay(
	int pixel_x,
	int pixel_y,
	float parameter_x,
	float parameter_y)
{
	Ray r;
	if (pixel_x < 0 || pixel_x > width_ - 1 ||
		pixel_y < 0 || pixel_y > height_ - 1 ||
		parameter_x < -0.5 || parameter_x > 0.5 ||
		parameter_y < -0.5 || parameter_y > 0.5
		)
	{
		r.origin = glm::vec3(0);
		r.direction = glm::vec3(0);
	}
	else
	{
		/*
		// View and perspective matrices are used in the unProject() function
		glm::mat4 V = glm::lookAt(eye_, center_, up_);
		float aspect = float(width_) / height_;
		glm::mat4 P = glm::perspective(fov_, aspect, 0.1f, 100.0f);
		*/

		// The unProject() function returns a vector in world-space which
		// defines a direction out of the frustum depending on which pixel
		// we shoot the ray from. "from" will be on the near-viewplane
		// and "to" will be on the far-viewplane.
		glm::vec4 from4 = VP_inv * glm::vec4(((pixel_x + parameter_x) / width_ - 0.5) * 2, ((pixel_y + parameter_y) / height_ - 0.5) * 2, 1, 1 ); /*
		glm::unProject(
			glm::vec3(pixel_x + parameter_x, pixel_y + parameter_y, 0.0f),
			V,
			P,
			glm::vec4(0, 0, width_, height_));
			*/
		glm::vec4 to4 = VP_inv * glm::vec4(((pixel_x + parameter_x) / width_ - 0.5) * 2, ((pixel_y + parameter_y) / height_ - 0.5) * 2, -1, 1 );/*glm::unProject(
			glm::vec3(pixel_x + parameter_x, pixel_y + parameter_y, 1.0f),
			V,
			P,
			glm::vec4(0, 0, width_, height_));*/

		glm::vec3 from = glm::vec3(from4) * from4.w;
		glm::vec3 to = glm::vec3(to4) * to4.w;

		glm::vec3 direction = glm::normalize(to - from);

		r.origin = eye_;
		r.direction = direction;
		 // Set air as the starting material for the ray to travel in.
		r.material = Material::air();
		r.radiance = SpectralDistribution(); // Importance
		r.radiance[0] = 1;
		r.radiance[1] = 1;
		r.radiance[2] = 1;
	}
	return r;
}