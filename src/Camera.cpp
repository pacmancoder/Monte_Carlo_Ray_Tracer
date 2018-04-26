#include <mcrt/Camera.h>

Camera::Camera (
	glm::vec3 eye,
	glm::vec3 center,
	glm::vec3 up,
	float fov,
	uint16_t baseWidth,
	uint16_t baseHeight) :

	eye_(eye),
	center_(center),
	up_(up),
	fov_(fov),
	width_(baseWidth),
	height_(baseHeight)
{
	// View and perspective matrices are used in the unProject() function
	glm::mat4 V = glm::lookAt(eye, center, up);
	glm::mat4 P = glm::perspective(fov, float(baseWidth) / baseHeight, 0.1f, 100.0f);
	VP_inv = glm::inverse(V * P);
}

Ray Camera::castRay(
	uint16_t frameWidth,
	uint16_t frameHeight,
	uint16_t pixel_x,
	uint16_t pixel_y,
	float parameter_x,
	float parameter_y) const
{
	Ray r;

	if (pixel_x > width_ - 1 ||
		pixel_y > height_ - 1 ||
		parameter_x < -0.5 || parameter_x > 0.5 ||
		parameter_y < -0.5 || parameter_y > 0.5)
	{
		r.origin = glm::vec3(0);
		r.direction = glm::vec3(0);
	}
	else
	{
		glm::vec4 from4 = VP_inv * glm::vec4(((pixel_x + parameter_x) / frameWidth - 0.5) * 2, ((pixel_y + parameter_y) / frameHeight - 0.5) * 2, 1, 1 );
		glm::vec4 to4 = VP_inv * glm::vec4(((pixel_x + parameter_x) / frameWidth - 0.5) * 2, ((pixel_y + parameter_y) / frameHeight - 0.5) * 2, -1, 1 );
		glm::vec3 from = glm::vec3(from4) * from4.w;
		glm::vec3 to = glm::vec3(to4) * to4.w;
		glm::vec3 direction = glm::normalize(to - from);

		r.origin = eye_;
		r.direction = direction;
		r.material = Material::air();
		r.radiance = SpectralDistribution();
		r.radiance[0] = 1;
		r.radiance[1] = 1;
		r.radiance[2] = 1;
	}

	return r;
}