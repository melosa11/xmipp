#ifndef CUDA_FORWARD_ART_ZERNIKE3D_TPP
#define CUDA_FORWARD_ART_ZERNIKE3D_TPP

// Xmipp includes
#include "cuda_forward_art_zernike3d.h"

template<typename PrecisionType>
CUDAForwardArtZernike3D<PrecisionType>::CUDAForwardArtZernike3D(
        const CUDAForwardArtZernike3D<PrecisionType>::ConstantParameters parameters) {
    (void) parameters;
}

template<typename PrecisionType>
CUDAForwardArtZernike3D<PrecisionType>::~CUDAForwardArtZernike3D() {

}

template<typename PrecisionType>
template<bool usesZernike>
void CUDAForwardArtZernike3D<PrecisionType>::runForwardKernel(
        const std::vector<PrecisionType> &clnm,
        std::vector<Image<PrecisionType>> &P,
        std::vector<Image<PrecisionType>> &W) {
    if (usesZernike) {
        return;
    }
}

template<typename PrecisionType>
template<bool usesZernike>
void CUDAForwardArtZernike3D<PrecisionType>::runBackwardKernel(const std::vector<PrecisionType> &clnm,
                                                               const Image<PrecisionType> &Idiff) {
    if (usesZernike) {
        return;
    }
}

#endif// CUDA_FORWARD_ART_ZERNIKE3D_TPP