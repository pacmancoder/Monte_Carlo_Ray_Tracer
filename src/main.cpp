#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <ctime>

#include <glm/glm.hpp>
#include <omp.h>

#include "../include/Camera.h"
#include "../include/Scene.h"

int main()
{
	time_t time_start, time_now, rendertime_start;
	time(&time_start);

	static const int WIDTH = 1024 / 1;
	static const int HEIGHT = 768 / 1;
	static const int SUB_SAMPLING_CAUSTICS = 10;
	static const int SUB_SAMPLING_MONTE_CARLO = 500;
	static const int SUB_SAMPLING_DIRECT_SPECULAR = 100;
	static const int NUMBER_OF_PHOTONS_EMISSION = 2000000;

	// The camera is used to cast appropriate initial rays
	Camera c(
		glm::vec3(0, 0, 3.2), // Eye (position of camera)
		glm::vec3(0, 0, 0), // Center (position to look at)
		glm::vec3(0, 1, 0), // Up vector
		M_PI / 3, // Field of view in radians
		WIDTH, // pixel width
		HEIGHT); // pixel height
	glm::vec3 camera_plane_normal = glm::normalize(c.center_ - c.eye_);

	// 3D objects are contained in the Scene object
	Scene s;

	// irradiance_values will hold image data
    std::vector<SpectralDistribution> irradiance_values(static_cast<size_t>(c.width_ * c.height_));

	// irradiance_values need to be converted to rgb pixel data for displaying
    std::vector<uint8_t> pixel_values(c.width_ * c.height_ * 3, 0);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(-0.5, 0.5);

	std::cout << "Building photon map." << std::endl;
	s.buildPhotonMap(NUMBER_OF_PHOTONS_EMISSION);

	float rendering_percent_finished = 0;
	std::cout << "Rendering started!" << std::endl;
	std::cout << rendering_percent_finished << " \% finished." << std::endl;

	time(&rendertime_start);

	double prerender_time = difftime(rendertime_start, time_start);

	// Loop through all pixels to calculate their irradiance_values by ray-tracing
	for (int x = 0; x < c.width_; ++x)
	{
		// Parallellize the for loop with openMP.
		#pragma omp parallel for
		for (int y = 0; y < c.height_; ++y)
		{
			int index = (x + y * c.width_);
			SpectralDistribution sd;
			if (SUB_SAMPLING_DIRECT_SPECULAR)
			{
				for (int i = 0; i < SUB_SAMPLING_DIRECT_SPECULAR; ++i)
				{
					Ray r = c.castRay(
						x, // Pixel x
						(c.height_ - y - 1), // Pixel y
						dis(gen), // Parameter x (>= -0.5 and < 0.5), for subsampling
						dis(gen)); // Parameter y (>= -0.5 and < 0.5), for subsampling
					sd += s.traceRay(r, Scene::RenderMode::WHITTED_SPECULAR) * glm::dot(r.direction, camera_plane_normal);
				}
				irradiance_values[index] += sd / SUB_SAMPLING_DIRECT_SPECULAR * (2 * M_PI);				
			}
			sd = SpectralDistribution();
			if (SUB_SAMPLING_CAUSTICS)
				{
				for (int i = 0; i < SUB_SAMPLING_CAUSTICS; ++i)
				{
					Ray r = c.castRay(
						x, // Pixel x
						(c.height_ - y - 1), // Pixel y
						dis(gen), // Parameter x (>= -0.5 and < 0.5), for subsampling
						dis(gen)); // Parameter y (>= -0.5 and < 0.5), for subsampling
					sd += s.traceRay(r, Scene::RenderMode::CAUSTICS) * glm::dot(r.direction, camera_plane_normal);
				}
				irradiance_values[index] += sd / SUB_SAMPLING_CAUSTICS * (2 * M_PI);
			}
			sd = SpectralDistribution();
			if (SUB_SAMPLING_MONTE_CARLO)
				{
				for (int i = 0; i < SUB_SAMPLING_MONTE_CARLO; ++i)
				{
					Ray r = c.castRay(
						x, // Pixel x
						(c.height_ - y - 1), // Pixel y
						dis(gen), // Parameter x (>= -0.5 and < 0.5), for subsampling
						dis(gen)); // Parameter y (>= -0.5 and < 0.5), for subsampling
					sd += s.traceRay(r, Scene::RenderMode::MONTE_CARLO) * glm::dot(r.direction, camera_plane_normal);
				}
				irradiance_values[index] += sd / SUB_SAMPLING_MONTE_CARLO * (2 * M_PI);
			}
		}

		// To show how much time we have left.
		rendering_percent_finished = (x+1) * 100 / float(c.width_);
	  	time(&time_now);
		double rendering_time_elapsed = difftime(time_now, rendertime_start);
		double rendering_time_left = (rendering_time_elapsed / rendering_percent_finished) *
			(100 - rendering_percent_finished);

		double total_time_elapsed = prerender_time + rendering_time_elapsed;
		double total_time_estimate = total_time_elapsed + rendering_time_left;
		double total_time_left = total_time_estimate - total_time_elapsed;

		int hours = total_time_left / (60 * 60);
		int minutes = (int(total_time_left) % (60 * 60)) / 60;
		int seconds = int(total_time_left) % 60;

		std::cout << rendering_percent_finished << " \% of rendering finished." << std::endl;
		std::cout << "Estimated time left is "
			<< hours << "h:"
			<< minutes << "m:"
			<< seconds << "s." << std::endl;
	}

	// To show how much time it actually took to render.
	time(&time_now);
	double time_elapsed = difftime(time_now, time_start);
	int hours_elapsed = time_elapsed / (60 * 60);
	int minutes_elapsed = (int(time_elapsed) % (60 * 60)) / 60;
	int seconds_elapsed = int(time_elapsed) % 60;

	int hours_prerender = prerender_time / (60 * 60);
	int minutes_prerender = (int(prerender_time) % (60 * 60)) / 60;
	int seconds_prerender = int(prerender_time) % 60;

	std::string rendering_time_string =
	 	  std::to_string(hours_elapsed) + "h:"
		+ std::to_string(minutes_elapsed) + "m:"
		+ std::to_string(seconds_elapsed) + "s";

	std::string prerender_time_string =
		  std::to_string(hours_prerender) + "h:"
		+ std::to_string(minutes_prerender) + "m:"
		+ std::to_string(seconds_prerender) + "s";

	std::cout << "Rendering time : " << rendering_time_string << std::endl;

	// Convert to byte data
	// Gamma correction
	auto gamma = static_cast<float>(1 / 2.2);
	for (int x = 0; x < c.width_; ++x)
	{
		for (int y = 0; y < c.height_; ++y)
		{
			int index = (x + y * c.width_);
			pixel_values[index * 3 + 0] = static_cast<unsigned char>(
			        glm::clamp(glm::pow(irradiance_values[index][0],gamma), 0.0f, 1.0f) * 255); // Red
            pixel_values[index * 3 + 1] = static_cast<unsigned char>(
                    glm::clamp(glm::pow(irradiance_values[index][1],gamma), 0.0f, 1.0f) * 255); // Green
            pixel_values[index * 3 + 2] = static_cast<unsigned char>(
                    glm::clamp(glm::pow(irradiance_values[index][2],gamma), 0.0f, 1.0f) * 255); // Blue
		}
	}

    return 0;
}