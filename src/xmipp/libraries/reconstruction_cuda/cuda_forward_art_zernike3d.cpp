// Xmipp includes
#include "cuda_forward_art_zernike3d.h"
// Standard includes
#include <iostream>

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

// We need tell the compiler the types for which we want
// to compile, because we want the definition and implementation
// to be in separate files.

template class CUDAForwardArtZernike3D<float>;
template class CUDAForwardArtZernike3D<double>;
