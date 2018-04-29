#pragma once

#include <mcrt/Scene.h>
#include <mcrt/Camera.h>
#include <mcrt/utils.h>

namespace Mcrt
{
    struct FrameInfo
    {
        uint16_t width;
        uint16_t height;
    };

    struct ChunkInfo
    {
        uint16_t x;
        uint16_t y;
        uint16_t width;
        uint16_t height;
    };

    class SurfaceRenderer
    {
    public:
        std::vector<SpectralDistribution> RenderChunk(
                Scene& scene, const FrameInfo& frameInfo, const ChunkInfo& chunkInfo) const;
    };
}