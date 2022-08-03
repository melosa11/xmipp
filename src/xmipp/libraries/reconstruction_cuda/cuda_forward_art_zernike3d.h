#ifndef CUDA_FORWARD_ART_ZERNIKE3D_H
#define CUDA_FORWARD_ART_ZERNIKE3D_H

// Xmipp includes
#include <core/xmipp_image.h>
#include <core/multidim_array.h>
#include <core/xmipp_image.h>
#include <core/matrix1d.h>
#include <core/multidim_array.h>
// Standard includes
#include <vector>
#include <memory>
#include <atomic>

struct float3;
struct double3;

template<typename T>
struct MultidimArrayCuda {
    size_t xdim;
    size_t ydim;
    size_t yxdim;
    int xinit;
    int yinit;
    int zinit;
    T* data;
};

template<typename PrecisionType = float>
class CUDAForwardArtZernike3D
{
    static_assert(std::is_floating_point<PrecisionType>::value, "Floating point type is required.");

    using PrecisionType3 = std::conditional<std::is_same<PrecisionType, float>::value, float3, double3>;

public:
    /// Constant parameters for the computation
    struct ConstantParameters {
        Image<PrecisionType> &Vrefined;
        MultidimArray<int> &VRecMask, &sphMask;
        Matrix1D<int> &vL1, &vN, &vL2, &vM;
        std::vector<PrecisionType> &sigma;
        int RmaxDef;
        int loopStep;
        PrecisionType rot, tilt, psi;
        std::vector<std::unique_ptr<std::atomic<PrecisionType*>>> &p_busy_elem;
        std::vector<std::unique_ptr<std::atomic<PrecisionType*>>> &w_busy_elem;
    };

public:

    template<bool usesZernike>
    void runForwardKernel(const std::vector<PrecisionType> &clnm,
                          std::vector<Image<PrecisionType>> &P,
                          std::vector<Image<PrecisionType>> &W);

    template<bool usesZernike>
    void runBackwardKernel(const std::vector<PrecisionType> &clnm,
                           const Image<PrecisionType> &Idiff);

    CUDAForwardArtZernike3D(const ConstantParameters parameters);
    ~CUDAForwardArtZernike3D();

private:

    MultidimArrayCuda<PrecisionType> V;

    MultidimArrayCuda<int> VRecMask, sphMask;

    int RmaxDef;

    int loopStep;

    int lastX, lastY, lastZ;

    PrecisionType rot, tilt, psi;

    Matrix1D<int> vL1, vN, vL2, vM;

    std::vector<PrecisionType> sigma;

    std::vector<std::unique_ptr<std::atomic<PrecisionType*>>> p_busy_elem;
    std::vector<std::unique_ptr<std::atomic<PrecisionType*>>> w_busy_elem;

private:

    /// Move data from MultidimArray to struct usable by CUDA kernel
    template<typename T>
    MultidimArrayCuda<T> initializeMultidimArray(MultidimArray<T> &multidimArray) const;

    /// Function inspired by std::find with support for CUDA allowed data types
    size_t findCuda(PrecisionType *begin, size_t size, PrecisionType value);

    void splattingAtPos(PrecisionType pos_x, PrecisionType pos_y, PrecisionType weight,
                        MultidimArrayCuda<PrecisionType> &mP, MultidimArrayCuda<PrecisionType> &mW,
                        std::unique_ptr<std::atomic<PrecisionType *>> *p_busy_elem_cuda,
                        std::unique_ptr<std::atomic<PrecisionType *>> *w_busy_elem_cuda);
};

// Include template implementation
#include "cuda_forward_art_zernike3d.tpp"

#endif// CUDA_FORWARD_ART_ZERNIKE3D_H
