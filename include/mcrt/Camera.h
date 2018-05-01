#pragma once

#include <mcrt/utils.h>

namespace Mcrt
{
    class Camera
    {
    public:
        Camera(
                glm::vec3 pos,
                glm::vec3 up,
                glm::vec3 direction,
                float fovY,
                float focalDistance,
                float apertureRadius);

        Ray cast_ray(
                float x,
                float y,
                uint16_t width,
                uint16_t height,
                uint16_t u,
                uint16_t v
        ) const;

        glm::vec3 get_pos() const;
        glm::vec3 get_direction() const;

    private:
        glm::vec3 pos_;
        glm::vec3 u_, v_, w_;
        float m_;
        float focalDistance_;
        float apertureRadius_;
    };
}