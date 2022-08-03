#ifndef CUDA_FORWARD_ART_ZERNIKE3D_TPP
#define CUDA_FORWARD_ART_ZERNIKE3D_TPP

// Xmipp includes
#include "cuda_forward_art_zernike3d.h"
#include <core/geometry.h>
#include "data/numerical_tools.h"

// Macros
#define IS_OUTSIDE2D(ImD,i,j) \
    ((j) < STARTINGX((ImD)) || (j) > FINISHINGX((ImD)) || \
     (i) < STARTINGY((ImD)) || (i) > FINISHINGY((ImD)))

template<typename PrecisionType>
CUDAForwardArtZernike3D<PrecisionType>::CUDAForwardArtZernike3D(
        const CUDAForwardArtZernike3D<PrecisionType>::ConstantParameters parameters) {

    V = initializeMultidimArray(parameters.Vrefined());
    VRecMask = initializeMultidimArray(parameters.VRecMask);
    sphMask = initializeMultidimArray(parameters.sphMask);
    sigma = parameters.sigma;
    RmaxDef = parameters.RmaxDef;

    lastZ = FINISHINGZ(parameters.Vrefined());
    lastY = FINISHINGY(parameters.Vrefined());
    lastX = FINISHINGX(parameters.Vrefined());

    vL1 = parameters.vL1;
    vL2 = parameters.vL2;
    vN = parameters.vN;
    vM = parameters.vM;
}

template<typename PrecisionType>
CUDAForwardArtZernike3D<PrecisionType>::~CUDAForwardArtZernike3D() {

}

template<typename PrecisionType>
template<bool usesZernike>
void CUDAForwardArtZernike3D<PrecisionType>::runForwardKernel(
        const std::vector<PrecisionType> &clnm,
        std::vector<Image<PrecisionType>> &P,
        std::vector<Image<PrecisionType>> &W,
        PrecisionType rot,
        PrecisionType tilt,
        PrecisionType psi) {
    auto &mV = V;
    const size_t idxY0 = usesZernike ? (clnm.size() / 3) : 0;
    const size_t idxZ0 = usesZernike ? (2 * idxY0) : 0;
    const PrecisionType RmaxF = usesZernike ? RmaxDef : 0;
    const PrecisionType iRmaxF = usesZernike ? (1.0f / RmaxF) : 0;

    // Rotation Matrix
    const Matrix2D<PrecisionType> R = createRotationMatrix(rot, tilt, psi);

    // Setup data for CUDA kernel
    auto &cudaVRecMask = VRecMask;
    auto cudaMV = mV;
    std::vector<MultidimArrayCuda<PrecisionType>> tempP;
    std::vector<MultidimArrayCuda<PrecisionType>> tempW;
    for (int m = 0; m < P.size(); m++)
    {
        tempP.push_back(initializeMultidimArray(P[m]()));
    }
    for (int m = 0; m < W.size(); m++)
    {
        tempW.push_back(initializeMultidimArray(W[m]()));
    }
    auto cudaP = tempP.data();
    auto cudaW = tempW.data();
    auto cudaVL1 = vL1.vdata;
    auto cudaVN = vN.vdata;
    auto cudaVL2 = vL2.vdata;
    auto cudaVM = vM.vdata;
    auto cudaClnm = clnm.data();
    auto cudaR = R.mdata;
    auto sigma_size = sigma.size();
    auto cudaSigma = sigma.data();
    auto p_busy_elem_cuda = p_busy_elem.data();
    auto w_busy_elem_cuda = w_busy_elem.data();

    const auto lastZ = this->lastZ;
    const auto lastY = this->lastY;
    const auto lastX = this->lastX;
    const int step = loopStep;
    for (int k = STARTINGZ(cudaMV); k <= lastZ; k += step)
    {
        for (int i = STARTINGY(cudaMV); i <= lastY; i += step)
        {
            for (int j = STARTINGX(cudaMV); j <= lastX; j += step)
            {
                // Future CUDA code
                PrecisionType gx = 0.0, gy = 0.0, gz = 0.0;
                if (A3D_ELEM(cudaVRecMask, k, i, j) != 0)
                {
                    int img_idx = 0;
                    if (sigma_size > 1)
                    {
                        PrecisionType sigma_mask = A3D_ELEM(cudaVRecMask, k, i, j);
                        img_idx = findCuda(cudaSigma, sigma_size, sigma_mask);
                    }
                    auto &mP = cudaP[img_idx];
                    auto &mW = cudaW[img_idx];
                    if (usesZernike)
                    {
                        auto k2 = k * k;
                        auto kr = k * iRmaxF;
                        auto k2i2 = k2 + i * i;
                        auto ir = i * iRmaxF;
                        auto r2 = k2i2 + j * j;
                        auto jr = j * iRmaxF;
#define SQRT sqrtf
                        auto rr = SQRT(r2) * iRmaxF;
#undef SQRT
                        for (size_t idx = 0; idx < idxY0; idx++)
                        {
                            auto l1 = cudaVL1[idx];
                            auto n = cudaVN[idx];
                            auto l2 = cudaVL2[idx];
                            auto m = cudaVM[idx];
                            if (rr > 0 || l2 == 0)
                            {
                                PrecisionType zsph = ZernikeSphericalHarmonics(l1, n, l2, m, jr, ir, kr, rr);
                                gx += cudaClnm[idx] * (zsph);
                                gy += cudaClnm[idx + idxY0] * (zsph);
                                gz += cudaClnm[idx + idxZ0] * (zsph);
                            }
                        }
                    }

                    auto r_x = j + gx;
                    auto r_y = i + gy;
                    auto r_z = k + gz;

                    auto pos_x = cudaR[0] * r_x + cudaR[1] * r_y + cudaR[2] * r_z;
                    auto pos_y = cudaR[3] * r_x + cudaR[4] * r_y + cudaR[5] * r_z;
                    PrecisionType voxel_mV = A3D_ELEM(cudaMV, k, i, j);
                    splattingAtPos(pos_x, pos_y, voxel_mV, mP, mW, p_busy_elem_cuda, w_busy_elem_cuda);
                }
                // End of future CUDA code
            }
        }
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

template<typename PrecisionType>
template<typename T>
MultidimArrayCuda<T> CUDAForwardArtZernike3D<PrecisionType>::initializeMultidimArray(MultidimArray<T> &multidimArray) const {
    struct MultidimArrayCuda<T> cudaArray = {
            .xdim = multidimArray.xdim,
            .ydim = multidimArray.ydim,
            .yxdim = multidimArray.yxdim,
            .xinit = multidimArray.xinit,
            .yinit = multidimArray.yinit,
            .zinit = multidimArray.zinit,
            .data = multidimArray.data
    };
    return cudaArray;
}

template<typename PrecisionType>
void CUDAForwardArtZernike3D<PrecisionType>::splattingAtPos(PrecisionType pos_x, PrecisionType pos_y, PrecisionType weight,
                                                MultidimArrayCuda<PrecisionType> &mP, MultidimArrayCuda<PrecisionType> &mW,
                                                std::unique_ptr<std::atomic<PrecisionType *>> *p_busy_elem_cuda,
                                                std::unique_ptr<std::atomic<PrecisionType *>> *w_busy_elem_cuda)
{
    int i = round(pos_y);
    int j = round(pos_x);
    if(!IS_OUTSIDE2D(mP, i, j))
    {
        int idy = (i)-STARTINGY(mP);
        int idx = (j)-STARTINGX(mP);
        int idn = (idy) * (mP).xdim + (idx);
        // Not sure if std::unique_ptr and std::atomic can be used in CUDA code
        while ((*p_busy_elem_cuda[idn]) == &A2D_ELEM(mP, i, j));
        (*p_busy_elem_cuda[idn]).exchange(&A2D_ELEM(mP, i, j));
        (*w_busy_elem_cuda[idn]).exchange(&A2D_ELEM(mW, i, j));
        A2D_ELEM(mP, i, j) += weight;
        A2D_ELEM(mW, i, j) += 1.0;
        (*p_busy_elem_cuda[idn]).exchange(nullptr);
        (*w_busy_elem_cuda[idn]).exchange(nullptr);
    }
}

template<typename PrecisionType>
size_t CUDAForwardArtZernike3D<PrecisionType>::findCuda(PrecisionType *begin, size_t size, PrecisionType value)
{
    if (size <= 0)
    {
        return 0;
    }
    for (size_t i = 0; i < size; i++)
    {
        if (begin[i] == value)
        {
            return i;
        }
    }
    return size - 1;
}

template<typename PrecisionType>
Matrix2D<PrecisionType> CUDAForwardArtZernike3D<PrecisionType>::createRotationMatrix(PrecisionType rot,
                                                                                     PrecisionType tilt,
                                                                                     PrecisionType psi) {
    constexpr size_t matrixSize = 3;
    auto tmp = Matrix2D<PrecisionType>();
    tmp.initIdentity(matrixSize);
    Euler_angles2matrix(rot, tilt, psi, tmp, false);
    return tmp;
}

#endif// CUDA_FORWARD_ART_ZERNIKE3D_TPP
