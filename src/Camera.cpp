#include <mcrt/Camera.h>

#include <cmath>

#include <mcrt/utils.h>

using namespace Mcrt;

Camera::Camera(
		glm::vec3 pos,
		glm::vec3 up,
		glm::vec3 direction,
		float fovY,
		float focalDistance,
		float apertureRadius)
{
	pos_ = pos;

	w_ = direction;
	u_ = glm::normalize(glm::cross(up , w_));
	v_ = glm::normalize(glm::cross(w_ , u_));

	m_ = float(1.0 / std::tan(fovY / 2.0));

	focalDistance_ = focalDistance;
	apertureRadius_ = apertureRadius;
}

glm::vec3 Camera::get_pos() const
{
	return pos_;
}

glm::vec3 Camera::get_direction() const
{
	return w_;
}

Ray Camera::cast_ray(
		float x,
		float y,
		uint16_t width,
		uint16_t height,
		float u,
		float v) const
{
	float aspect = float(width) / float(height);
	float px = ((x + u - float(0.5)) / (float(width) - 1)) * 2 - 1;
	float py = ((y + v - float(0.5)) / (float(height) - 1)) * 2 - 1;
	glm::vec3 rayDir = glm::normalize(u_ * (-px * aspect) + v_ * (-py) + w_ * m_);
	glm::vec3 rayPos = pos_;

	/*
    if (apertureRadius_ > 0)
    {
        Types::Vec3 focal_point = pos_ + rayDir * focalDistance_;
        Types::Real angle = rng.generate_real(0.0, 1.0) * 2 * Constants::Math::Pi;
        Types::Real radius = rng.generate_real(0.0, 1.0) * apertureRadius_;
        rayPos += u_ * (std::cos(angle) * radius) + v_ * (std::sin(angle) * radius);
        rayDir = (focal_point - rayPos).normalize();
    }
    */

	return Ray {rayPos, rayDir, Material::air(), SpectralDistribution(1, 1, 1), false};
}
