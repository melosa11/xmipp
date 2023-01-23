// Xmipp includes
#include "cuda_forward_art_zernike3d.h"
#include <core/geometry.h>
#include "cuda_forward_art_zernike3d.cu"
#include "data/numerical_tools.h"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <tuple>
#include <utility>
#include "data/numerical_tools.h"

namespace cuda_forward_art_zernike3D {

// Cuda memory helper function
namespace {

	void processCudaError()
	{
		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess) {
			throw std::runtime_error(cudaGetErrorString(err));
		}
	}

	// Copies data from CPU to the GPU
	template<typename T>
	void transportData(T **dest, const T *source, size_t n)
	{
		if (cudaMalloc(dest, sizeof(T) * n) != cudaSuccess) {
			processCudaError();
		}

		if (cudaMemcpy(*dest, source, sizeof(T) * n, cudaMemcpyHostToDevice) != cudaSuccess) {
			cudaFree(*dest);
			processCudaError();
		}
	}

	// Copies data from GPU to the CPU
	template<typename T>
	void transportDataFromGPU(T *dest, const T *source, size_t n)
	{
		if (cudaMemcpy(dest, source, sizeof(T) * n, cudaMemcpyDeviceToHost) != cudaSuccess) {
			processCudaError();
		}
	}

	template<typename T>
	T *transportMultidimArrayToGpu(const MultidimArray<T> &inputArray)
	{
		T *outputArrayData;
		transportData(&outputArrayData, inputArray.data, inputArray.xdim * inputArray.ydim * inputArray.zdim);
		return outputArrayData;
	}

	template<typename T>
	MultidimArrayCuda<T> *transportVectorOfMultidimArrayToGpu(const std::vector<MultidimArrayCuda<T>> &inputVector)
	{
		MultidimArrayCuda<T> *outputVectorData;
		transportData(&outputVectorData, inputVector.data(), inputVector.size());
		return outputVectorData;
	}

	template<typename T>
	T *transportMatrix1DToGpu(const Matrix1D<T> &inputVector)
	{
		T *outputVector;
		transportData(&outputVector, inputVector.vdata, inputVector.vdim);
		return outputVector;
	}

	template<typename T>
	T *transportStdVectorToGpu(const std::vector<T> &inputVector)
	{
		T *outputVector;
		transportData(&outputVector, inputVector.data(), inputVector.size());
		return outputVector;
	}

	template<typename T>
	T *transportMatrix2DToGpu(const Matrix2D<T> &inputMatrix)
	{
		T *outputMatrixData;
		transportData(&outputMatrixData, inputMatrix.mdata, inputMatrix.mdim);
		return outputMatrixData;
	}

	template<typename T>
	MultidimArrayCuda<T> initializeMultidimArrayCuda(const MultidimArray<T> &multidimArray)
	{
		struct MultidimArrayCuda<T> cudaArray = {
			.xdim = multidimArray.xdim, .ydim = multidimArray.ydim, .yxdim = multidimArray.yxdim,
			.zdim = multidimArray.zdim, .xinit = multidimArray.xinit, .yinit = multidimArray.yinit,
			.zinit = multidimArray.zinit, .data = transportMultidimArrayToGpu(multidimArray)
		};

		return cudaArray;
	}

	template<typename T>
	void updateMultidimArrayWithGPUData(MultidimArray<T> &multidimArray, const MultidimArrayCuda<T> &multidimArrayCuda)
	{
		transportDataFromGPU(
			multidimArray.data, multidimArrayCuda.data, multidimArray.xdim * multidimArray.ydim * multidimArray.zdim);
	}

	template<typename T>
	void updateVectorOfMultidimArrayWithGPUData(std::vector<Image<T>> &image,
												const std::vector<MultidimArrayCuda<T>> vectorMultidimArray)
	{
		assert(image.size() == vectorMultidimArray.size());
		for (int m = 0; m < image.size(); m++) {
			updateMultidimArrayWithGPUData(image[m](), vectorMultidimArray[m]);
		}
	}

	template<typename T>
	std::pair<MultidimArrayCuda<T> *, std::vector<MultidimArrayCuda<T>>> convertToMultidimArrayCuda(
		std::vector<Image<T>> &image)
	{
		std::vector<MultidimArrayCuda<T>> output;
		for (int m = 0; m < image.size(); m++) {
			output.push_back(initializeMultidimArrayCuda(image[m]()));
		}
		return std::make_pair(transportVectorOfMultidimArrayToGpu(output), output);
	}

	template<typename T>
	void freeCommonArgumentsKernel(struct Program<T>::CommonKernelParameters &commonParameters)
	{
		cudaFree(commonParameters.cudaClnm);
		cudaFree(commonParameters.cudaR);
	}

	template<typename T>
	void freeVectorOfMultidimArray(std::vector<MultidimArrayCuda<T>> vector)
	{
		for (int m = 0; m < vector.size(); m++) {
			cudaFree(vector[m].data);
		}
	}

	template<typename T>
	Matrix2D<T> createRotationMatrix(const struct Program<T>::AngleParameters angles)
	{
		auto rot = angles.rot;
		auto tilt = angles.tilt;
		auto psi = angles.psi;
		constexpr size_t matrixSize = 3;
		auto tmp = Matrix2D<T>();
		tmp.initIdentity(matrixSize);
		Euler_angles2matrix(rot, tilt, psi, tmp, false);
		return tmp;
	}

	template<typename T>
	struct Program<T>::CommonKernelParameters getCommonArgumentsKernel(
		const struct Program<T>::DynamicParameters &parameters,
		const bool usesZernike,
		const int RmaxDef) {
		auto clnm = parameters.clnm;
		auto angles = parameters.angles;

		const size_t idxY0 = usesZernike ? (clnm.size() / 3) : 0;
		const size_t idxZ0 = usesZernike ? (2 * idxY0) : 0;
		const T RmaxF = usesZernike ? RmaxDef : 0;
		const T iRmaxF = usesZernike ? (1.0f / RmaxF) : 0;

		const Matrix2D<T> R = createRotationMatrix<T>(angles);

		struct Program<T>::CommonKernelParameters output = {
			.idxY0 = idxY0, .idxZ0 = idxZ0, .iRmaxF = iRmaxF, .cudaClnm = transportStdVectorToGpu(clnm),
			.cudaR = transportMatrix2DToGpu(R), .R = R,
		};

		return output;

	}

	struct Size3D
	blockSizeArchitecture()

	{
		return {4, 4, 4};
	}

	enum BackwardKernelMode initializeBackwardMode(bool computeMask, bool usesZernike)
	{
		if (computeMask) {
			return BackwardKernelMode::computeMask;
		}
		if (!usesZernike) {
			return BackwardKernelMode::sharedBlockMask;
		}
		return BackwardKernelMode::trivial;
	}

	Size3D blockSizeBackwardMode(enum BackwardKernelMode mode)
	{
		switch (mode) {
			case BackwardKernelMode::trivial:
				return {4, 4, 4};
			case BackwardKernelMode::computeMask:
				return {8, 8, 2};
			case BackwardKernelMode::sharedBlockMask:
				return {8, 8, 4};
		}
	}

	template<typename T>
	MultidimArray<T> alingMask(const MultidimArray<T> &mask, const Size3D threadBlockDim, int step = 1)
	{
		auto resizedMask = MultidimArray<int>(mask);

		auto newxDim = mask.xdim + mask.xdim % (threadBlockDim.x * step);
		auto newyDim = mask.ydim + mask.ydim % (threadBlockDim.y * step);
		auto newzDim = mask.zdim + mask.zdim % (threadBlockDim.z * step);

		// the padding will fill with zeroes
		resizedMask.resize(newzDim, newyDim, newxDim);
		return resizedMask;
	}

	VolumeMask initializeVolumeMask(const MultidimArray<int> &mask, enum BackwardKernelMode mode)
	{
		const auto blockSize = blockSizeBackwardMode(mode);
		const auto cudaMask =
			mode == BackwardKernelMode::computeMask
				? std::nullopt
				: std::optional<MultidimArrayCuda<int>>{initializeMultidimArrayCuda(alingMask(mask, blockSize))};
		return {
			.mask = cudaMask,
			.blockSize = blockSize,
			.gridSize = {mask.xdim / blockSize.x, mask.ydim / blockSize.y, mask.zdim / blockSize.z},
		};
	}

	template<typename PrecisionType>
	int *initializeBlockMask(const MultidimArrayCuda<PrecisionType> &cudaMV,
							 const VolumeMask &mask,
							 const MultidimArray<int> &resizedMask)
	{
		std::vector<int> v(mask.gridSize.z * mask.gridSize.y * mask.gridSize.x, 0);
		for (int k = 0; k < mask.gridSize.z; ++k) {
			for (int j = 0; j < mask.gridSize.y; ++j) {
				for (int i = 0; i < mask.gridSize.x; ++i) {

					int count = 0;
					for (int l = 0; l < mask.blockSize.z; ++l) {
						for (int m = 0; m < mask.blockSize.y; ++m) {
							for (int n = 0; n < mask.blockSize.x; ++n) {
								int z = STARTINGZ(cudaMV) + l + k * mask.blockSize.z;
								int y = STARTINGY(cudaMV) + m + j * mask.blockSize.y;
								int x = STARTINGX(cudaMV) + n + i * mask.blockSize.x;
								if (A3D_ELEM(resizedMask, z, y, x) != 0) {
									++count;
								}
							}
						}
					}
					const int mask_bit = count > 0 ? 1 : 0;
					v.at(i + j * mask.gridSize.x + k * mask.gridSize.x * mask.gridSize.y) = mask_bit;
				}
			}
		}
		return transportStdVectorToGpu(v);
	}

}  // namespace

template<typename PrecisionType>
Program<PrecisionType>::Program(const Program<PrecisionType>::ConstantParameters parameters)
	: cudaMV(initializeMultidimArrayCuda(parameters.Vrefined())),
	  VRecMaskF(
		  initializeMultidimArrayCuda(alingMask(parameters.VRecMaskF, blockSizeArchitecture(), parameters.loopStep))),
	  backwardMode(initializeBackwardMode(parameters.computeBackwardMask, parameters.usesZernike)),
	  backwardMask(initializeVolumeMask(parameters.VRecMaskB, backwardMode)),
	  sigma(parameters.sigma),
	  RmaxDef(parameters.RmaxDef),
	  lastX(FINISHINGX(parameters.Vrefined())),
	  lastY(FINISHINGY(parameters.Vrefined())),
	  lastZ(FINISHINGZ(parameters.Vrefined())),
	  loopStep(parameters.loopStep),
	  cudaVL1(transportMatrix1DToGpu(parameters.vL1)),
	  cudaVL2(transportMatrix1DToGpu(parameters.vL2)),
	  cudaVN(transportMatrix1DToGpu(parameters.vN)),
	  cudaVM(transportMatrix1DToGpu(parameters.vM)),
	  blockXStep(blockSizeArchitecture().x),
	  blockYStep(blockSizeArchitecture().y),
	  blockZStep(blockSizeArchitecture().z),
	  gridXStep(VRecMaskF.xdim / loopStep / blockXStep),
	  gridYStep(VRecMaskF.ydim / loopStep / blockYStep),
	  gridZStep(VRecMaskF.zdim / loopStep / blockZStep),
	  cudaBlockBackwardMask(initializeBlockMask(cudaMV,
												backwardMask,
												alingMask(parameters.VRecMaskB, blockSizeBackwardMode(backwardMode))))
{}

template<typename PrecisionType>
Program<PrecisionType>::~Program()
{
	cudaFree(VRecMaskF.data);
	if (backwardMask.mask.has_value()) {
		cudaFree(backwardMask.mask.value().data);
	}
	cudaFree(cudaMV.data);

	cudaFree(const_cast<int *>(cudaVL1));
	cudaFree(const_cast<int *>(cudaVL2));
	cudaFree(const_cast<int *>(cudaVN));
	cudaFree(const_cast<int *>(cudaVM));
}

template<typename PrecisionType>
template<bool usesZernike>
void Program<PrecisionType>::runForwardKernel(struct DynamicParameters &parameters)

{
	// Unique parameters
	MultidimArrayCuda<PrecisionType> *cudaP, *cudaW;
	std::vector<MultidimArrayCuda<PrecisionType>> pVector, wVector;
	std::tie(cudaP, pVector) = convertToMultidimArrayCuda(parameters.P);
	std::tie(cudaW, wVector) = convertToMultidimArrayCuda(parameters.W);
	auto sigma_size = sigma.size();
	auto cudaSigma = transportStdVectorToGpu(sigma);
	const int step = loopStep;

	// Common parameters
	auto commonParameters = getCommonArgumentsKernel<PrecisionType>(parameters, usesZernike, RmaxDef);

	forwardKernel<PrecisionType, usesZernike>
		<<<dim3(gridXStep, gridYStep, gridZStep), dim3(blockXStep, blockYStep, blockZStep)>>>(cudaMV,
																							  VRecMaskF,
																							  cudaP,
																							  cudaW,
																							  lastZ,
																							  lastY,
																							  lastX,
																							  step,
																							  sigma_size,
																							  cudaSigma,
																							  commonParameters.iRmaxF,
																							  commonParameters.idxY0,
																							  commonParameters.idxZ0,
																							  cudaVL1,
																							  cudaVN,
																							  cudaVL2,
																							  cudaVM,
																							  commonParameters.cudaClnm,
																							  commonParameters.cudaR);

	cudaDeviceSynchronize();

	updateVectorOfMultidimArrayWithGPUData(parameters.P, pVector);
	updateVectorOfMultidimArrayWithGPUData(parameters.W, wVector);

	freeVectorOfMultidimArray(pVector);
	freeVectorOfMultidimArray(wVector);
	cudaFree(cudaP);
	cudaFree(cudaW);
	cudaFree(cudaSigma);
	freeCommonArgumentsKernel<PrecisionType>(commonParameters);
}

template<typename PrecisionType>
template<bool usesZernike>
void Program<PrecisionType>::runBackwardKernel(struct DynamicParameters &parameters)
{
	// Unique parameters
	auto &mId = parameters.Idiff();
	auto cudaMId = initializeMultidimArrayCuda(mId);

	// Common parameters
	auto commonParameters = getCommonArgumentsKernel<PrecisionType>(parameters, usesZernike, RmaxDef);
	const auto commonArguments = CommonBackwardKernelArguments{
		.cudaMId = cudaMId,
		.r0 = commonParameters.R.mdata[0],
		.r1 = commonParameters.R.mdata[1],
		.r2 = commonParameters.R.mdata[2],
		.r3 = commonParameters.R.mdata[3],
		.r4 = commonParameters.R.mdata[4],
		.r5 = commonParameters.R.mdata[5],
	};
	const auto zernikeParameters = ZernikeParameters{
		.iRmaxF = commonParameters.iRmaxF,
		.idxY0 = commonParameters.idxY0,
		.idxZ0 = commonParameters.idxY0,
		.cudaVL1 = cudaVL1,
		.cudaVN = cudaVN,
		.cudaVL2 = cudaVL2,
		.cudaVM = cudaVM,
		.cudaClnm = commonParameters.cudaClnm,
	};

	const auto gridSize = dim3(backwardMask.gridSize.x, backwardMask.gridSize.y, backwardMask.gridSize.z);
	const auto blockSize = dim3(backwardMask.blockSize.x, backwardMask.blockSize.y, backwardMask.blockSize.z);

	if (backwardMode == BackwardKernelMode::sharedBlockMask) {
		const int diagonal = static_cast<int>(
			std::ceil(std::sqrt(blockSize.x * blockSize.x + blockSize.y * blockSize.y + blockSize.z * blockSize.z)));
		const int sharedMIdDim = 2 * ((diagonal / 2) + 1) + 1;
		const int sharedMemorySize = sharedMIdDim * sharedMIdDim * sizeof(PrecisionType);

		sharedBackwardKernel<PrecisionType>
			<<<gridSize, blockSize, sharedMemorySize>>>(cudaMV, cudaBlockBackwardMask, sharedMIdDim, commonArguments);
	} else if (backwardMode == BackwardKernelMode::computeMask) {
		computeBackwardKernel<PrecisionType, usesZernike>
			<<<gridSize, blockSize>>>(cudaMV, RmaxDef, commonArguments, zernikeParameters);
	} else {
		trivialBackwardKernel<PrecisionType, usesZernike>
			<<<gridSize, blockSize>>>(cudaMV, backwardMask.mask.value(), commonArguments, zernikeParameters);
	}
	cudaDeviceSynchronize();

	cudaFree(cudaMId.data);
	freeCommonArgumentsKernel<PrecisionType>(commonParameters);
}

template<typename PrecisionType>
void Program<PrecisionType>::recoverVolumeFromGPU(Image<PrecisionType> &Vrefined)
{
	updateMultidimArrayWithGPUData(Vrefined(), cudaMV);
}

// explicit template instantiation
template class Program<float>;
template void Program<float>::runForwardKernel<true>(struct DynamicParameters &);
template void Program<float>::runForwardKernel<false>(struct DynamicParameters &);
template void Program<float>::runBackwardKernel<true>(struct DynamicParameters &);
template void Program<float>::runBackwardKernel<false>(struct DynamicParameters &);
}  // namespace cuda_forward_art_zernike3D
