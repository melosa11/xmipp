// Xmipp includes
#include "cuda_forward_art_zernike3d.h"
#include <core/geometry.h>
#include "cuda_forward_art_zernike3d.cu"
#include "data/numerical_tools.h"

#include <algorithm>
#include <cassert>
#include <cmath>
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
			.xdim = static_cast<unsigned>(multidimArray.xdim), .ydim = static_cast<unsigned>(multidimArray.ydim),
			.yxdim = static_cast<unsigned>(multidimArray.yxdim), .xinit = multidimArray.xinit,
			.yinit = multidimArray.yinit, .zinit = multidimArray.zinit,
			.data = transportMultidimArrayToGpu(multidimArray)
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
			.idxY0 = idxY0, .idxZ0 = idxZ0, .iRmaxF = iRmaxF, .cudaClnm = transportStdVectorToGpu(clnm), .R = R,
		};

		return output;

	}

	enum BackwardKernelMode
	initializeBackwardMode(bool computeMask, bool usesZernike)

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

	template<typename T>
	cudaTextureObject_t initTextureMultidimArray(MultidimArrayCuda<T> &array, size_t zdim)

	{
		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeLinear;
		resDesc.res.linear.devPtr = array.data;
		resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
		resDesc.res.linear.desc.x = 32;
		resDesc.res.linear.sizeInBytes = zdim * array.yxdim * sizeof(T);
		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.readMode = cudaReadModeElementType;
		cudaTextureObject_t tex = 0;
		cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);
		return tex;
	}

	template<typename T>
	bool checkStep(MultidimArray<T> &mask, int step, size_t position)

	{
		if (position % mask.xdim % step != 0) {
			return false;
		}
		if (position / mask.xdim % mask.ydim % step != 0) {
			return false;
		}
		if (position / mask.yxdim % step != 0) {
			return false;
		}
		return mask[position] != 0;
	}

	template<typename T>
	std::tuple<unsigned *, size_t> filterMaskTransportCoordinates(MultidimArray<T> &mask, int step)

	{
		std::vector<unsigned> coordinates;
		for (unsigned i = 0; i < static_cast<unsigned>(mask.yxdim * mask.zdim); i++) {
			if (checkStep(mask, step, static_cast<size_t>(i))) {
				coordinates.push_back(i);
			}
		}
		unsigned *coordinatesCuda = transportStdVectorToGpu(coordinates);
		return std::make_tuple(coordinatesCuda, coordinates.size());
	}

	template<typename T>
	std::tuple<unsigned *, size_t, int *> filterMaskTransportCoordinates(MultidimArray<T> &mask,
																		 int step,
																		 bool transportValues)

	{
		std::vector<unsigned> coordinates;
		std::vector<T> values;
		for (unsigned i = 0; i < static_cast<unsigned>(mask.yxdim * mask.zdim); i++) {
			if (checkStep(mask, step, static_cast<size_t>(i))) {
				coordinates.push_back(i);
				if (transportValues) {
					values.push_back(mask[i]);
				}
			}
		}
		unsigned *coordinatesCuda = transportStdVectorToGpu(coordinates);
		int *valuesCuda = transportStdVectorToGpu(values);
		return std::make_tuple(coordinatesCuda, coordinates.size(), valuesCuda);
	}

}  // namespace

template<typename PrecisionType>
Program<PrecisionType>::Program(const Program<PrecisionType>::ConstantParameters parameters)
	: cudaMV(initializeMultidimArrayCuda(parameters.Vrefined())),
	  backwardMode(initializeBackwardMode(parameters.computeBackwardMask, parameters.usesZernike)),
	  backwardMask(initializeVolumeMask(parameters.VRecMaskB, backwardMode)),
	  sigma(parameters.sigma),
	  cudaSigma(transportStdVectorToGpu(parameters.sigma)),
	  RmaxDef(parameters.RmaxDef),
	  lastX(FINISHINGX(parameters.Vrefined())),
	  lastY(FINISHINGY(parameters.Vrefined())),
	  lastZ(FINISHINGZ(parameters.Vrefined())),
	  loopStep(parameters.loopStep),
	  cudaVL1(transportMatrix1DToGpu(parameters.vL1)),
	  cudaVL2(transportMatrix1DToGpu(parameters.vL2)),
	  cudaVN(transportMatrix1DToGpu(parameters.vN)),
	  cudaVM(transportMatrix1DToGpu(parameters.vM)),
	  cudaBlockBackwardMask(initializeBlockMask(cudaMV,
												backwardMask,
												alingMask(parameters.VRecMaskB, blockSizeBackwardMode(backwardMode)))),
	  xdimB(static_cast<unsigned>(parameters.VRecMaskB.xdim)),
	  ydimB(static_cast<unsigned>(parameters.VRecMaskB.ydim)),
	  xdimF(parameters.VRecMaskF.xdim),
	  ydimF(parameters.VRecMaskF.ydim)
{
	std::tie(cudaCoordinatesF, sizeF, VRecMaskF) =
		filterMaskTransportCoordinates(parameters.VRecMaskF, parameters.loopStep, true);
	const auto optimalizedSize = ceil(sizeF / BLOCK_SIZE) * BLOCK_SIZE;
	blockXStep = std::__gcd(BLOCK_SIZE, static_cast<int>(optimalizedSize));
	gridXStep = optimalizedSize / blockXStep;
}

template<typename PrecisionType>
Program<PrecisionType>::~Program()
{
	cudaFree(VRecMaskF);
	if (backwardMask.mask.has_value()) {
		cudaFree(backwardMask.mask.value().data);
	}
	cudaFree(cudaMV.data);
	cudaFree(cudaCoordinatesF);
	cudaFree(const_cast<PrecisionType *>(cudaSigma));

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
	unsigned sigma_size = static_cast<unsigned>(sigma.size());

	// Common parameters
	auto commonParameters = getCommonArgumentsKernel<PrecisionType>(parameters, usesZernike, RmaxDef);

	forwardKernel<PrecisionType, usesZernike><<<gridXStep, blockXStep>>>(cudaMV,
																		 VRecMaskF,
																		 cudaCoordinatesF,
																		 xdimF,
																		 ydimF,
																		 static_cast<unsigned>(sizeF),
																		 cudaP,
																		 cudaW,
																		 sigma_size,
																		 cudaSigma,
																		 commonParameters.iRmaxF,
																		 static_cast<unsigned>(commonParameters.idxY0),
																		 static_cast<unsigned>(commonParameters.idxZ0),
																		 cudaVL1,
																		 cudaVN,
																		 cudaVL2,
																		 cudaVM,
																		 commonParameters.cudaClnm,
																		 commonParameters.R.mdata[0],
																		 commonParameters.R.mdata[1],
																		 commonParameters.R.mdata[2],
																		 commonParameters.R.mdata[3],
																		 commonParameters.R.mdata[4],
																		 commonParameters.R.mdata[5]);

	cudaDeviceSynchronize();

	updateVectorOfMultidimArrayWithGPUData(parameters.P, pVector);
	updateVectorOfMultidimArrayWithGPUData(parameters.W, wVector);

	freeVectorOfMultidimArray(pVector);
	freeVectorOfMultidimArray(wVector);
	cudaFree(cudaP);
	cudaFree(cudaW);
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
