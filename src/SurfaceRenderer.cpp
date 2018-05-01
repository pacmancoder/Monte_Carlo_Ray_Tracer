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

#include <dcdr/logging/Logger.h>

using namespace Mcrt;
using namespace Dcdr::Logging;

std::vector<SpectralDistribution> SurfaceRenderer::RenderChunk(
        Scene& scene, const FrameInfo& frameInfo, const ChunkInfo& chunkInfo) const
{
	static const int SUB_SAMPLING_CAUSTICS = 0;
	static const int SUB_SAMPLING_MONTE_CARLO = 10;
	static const int SUB_SAMPLING_DIRECT_SPECULAR = 5;

	auto c = scene.GetDefaultCamera();

	glm::vec3 camera_plane_normal = c.get_direction();

    std::vector<SpectralDistribution> irradiance_values(
            static_cast<size_t>(chunkInfo.width * chunkInfo.height));

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0.0f, 1.0f);

	// Loop through all pixels to calculate their irradiance_values by ray-tracing
	for (uint16_t y = 0; y < chunkInfo.height; ++y)
	{
		// Parallellize the for loop with openMP.
		//#pragma omp parallel for
		for (uint16_t  x = 0; x < chunkInfo.width; ++x)
		{
			int index = (x + y * chunkInfo.width);

			SpectralDistribution sd;
			if (SUB_SAMPLING_DIRECT_SPECULAR > 0)
			{
				for (int i = 0; i < SUB_SAMPLING_DIRECT_SPECULAR; ++i)
				{
                    Ray r = c.cast_ray(
                            chunkInfo.x + x,  chunkInfo.y + y,
                            frameInfo.width, frameInfo.height,
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
                    Ray r = c.cast_ray(
                            chunkInfo.x + x,  chunkInfo.y + y,
                            frameInfo.width, frameInfo.height,
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
					Ray r = c.cast_ray(
                        chunkInfo.x + x,  chunkInfo.y + y,
                        frameInfo.width, frameInfo.height,
                        dis(gen), dis(gen));

					sd += scene.traceRay(r, Scene::RenderMode::MONTE_CARLO) * glm::dot(r.direction, camera_plane_normal);
				}
				irradiance_values[index] += sd / SUB_SAMPLING_MONTE_CARLO * (2 * M_PI);
			}
		}
	}
    return irradiance_values;
}
