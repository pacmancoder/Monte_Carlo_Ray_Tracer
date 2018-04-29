#include <mcrt/SurfaceRenderer.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <ctime>

#include <glm/glm.hpp>
#include <omp.h>

#include <mcrt/Camera.h>
#include <mcrt/Scene.h>

using namespace Mcrt;

std::vector<SpectralDistribution> SurfaceRenderer::RenderChunk(
        Scene& scene, const FrameInfo& frameInfo, const ChunkInfo& chunkInfo) const
{
	static const int SUB_SAMPLING_CAUSTICS = 10;
	static const int SUB_SAMPLING_MONTE_CARLO = 500;
	static const int SUB_SAMPLING_DIRECT_SPECULAR = 100;

	auto c = scene.GetDefaultCamera();

	glm::vec3 camera_plane_normal = glm::normalize(c.center_ - c.eye_);

    std::vector<SpectralDistribution> irradiance_values(
            static_cast<size_t>(chunkInfo.width * chunkInfo.height));

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(-0.5f, 0.5f);

	// Loop through all pixels to calculate their irradiance_values by ray-tracing
	for (uint16_t x = 0; x < chunkInfo.width; ++x)
	{
		// Parallellize the for loop with openMP.
		//#pragma omp parallel for
		for (uint16_t  y = 0; y < chunkInfo.height; ++y)
		{
			int index = (x + y * chunkInfo.width);
			SpectralDistribution sd;
			if (SUB_SAMPLING_DIRECT_SPECULAR > 0)
			{
				for (int i = 0; i < SUB_SAMPLING_DIRECT_SPECULAR; ++i)
				{
                    Ray r = c.castRay(
                            frameInfo.width, frameInfo.height,
                            x,  (frameInfo.height - chunkInfo.y - y - uint16_t(1)),
                            dis(gen), dis(gen));
					sd += scene.traceRay(r, Scene::RenderMode::WHITTED_SPECULAR) * glm::dot(r.direction, camera_plane_normal);
				}
				irradiance_values[index] += sd / SUB_SAMPLING_DIRECT_SPECULAR * (2 * M_PI);				
			}
			sd = SpectralDistribution();
			if (SUB_SAMPLING_CAUSTICS > 0)
				{
				for (int i = 0; i < SUB_SAMPLING_CAUSTICS; ++i)
				{
					Ray r = c.castRay(
					    frameInfo.width, frameInfo.height,
                        x,  (frameInfo.height - chunkInfo.y - y - uint16_t(1)),
						dis(gen), dis(gen));
					sd += scene.traceRay(r, Scene::RenderMode::CAUSTICS) * glm::dot(r.direction, camera_plane_normal);
				}
				irradiance_values[index] += sd / SUB_SAMPLING_CAUSTICS * (2 * M_PI);
			}
			sd = SpectralDistribution();
			if (SUB_SAMPLING_MONTE_CARLO > 0)
				{
				for (int i = 0; i < SUB_SAMPLING_MONTE_CARLO; ++i)
				{
					Ray r = c.castRay(
                        frameInfo.width, frameInfo.height,
                        x,  (frameInfo.height - chunkInfo.y - y - uint16_t(1)),
                        dis(gen), dis(gen));

					sd += scene.traceRay(r, Scene::RenderMode::MONTE_CARLO) * glm::dot(r.direction, camera_plane_normal);
				}
				irradiance_values[index] += sd / SUB_SAMPLING_MONTE_CARLO * (2 * M_PI);
			}
		}
	}

    /* Gamma correction -- place to the server ?
	auto gamma = static_cast<float>(1 / 2.2);
	for (int x = 0; x < c.width_; ++x)
	{
		for (int y = 0; y < c.height_; ++y)
		{
			int index = (x + y * c.width_);
			pixel_values[index * 3 + 0] = static_cast<unsigned char>(
			        glm::clamp(glm::pow(irradiance_values[index][0], gamma), 0.0f, 1.0f) * 255); // Red
            pixel_values[index * 3 + 1] = static_cast<unsigned char>(
                    glm::clamp(glm::pow(irradiance_values[index][1], gamma), 0.0f, 1.0f) * 255); // Green
            pixel_values[index * 3 + 2] = static_cast<unsigned char>(
                    glm::clamp(glm::pow(irradiance_values[index][2], gamma), 0.0f, 1.0f) * 255); // Blue
		}
	}
    */
    return irradiance_values;
}
