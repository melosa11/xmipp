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
	void transportData(T **dest, const T *source, size_t n, cudaStream_t &stream)
	{
		if (cudaMalloc(dest, sizeof(T) * n) != cudaSuccess) {
			processCudaError();
		}

		if (cudaMemcpyAsync(*dest, source, sizeof(T) * n, cudaMemcpyHostToDevice) != cudaSuccess) {
			cudaFree(*dest);
			processCudaError();
		}
	}

	// Copies data from GPU to the CPU
	template<typename T>
	void transportDataFromGPU(T *dest, const T *source, size_t n, cudaStream_t &stream)
	{
		if (cudaMemcpyAsync(dest, source, sizeof(T) * n, cudaMemcpyDeviceToHost) != cudaSuccess) {
			processCudaError();
		}
	}

	template<typename T>
	T *transportMultidimArrayToGpu(const MultidimArray<T> &inputArray, cudaStream_t &stream)
	{
		T *outputArrayData;
		transportData(&outputArrayData, inputArray.data, inputArray.xdim * inputArray.ydim * inputArray.zdim, stream);
		return outputArrayData;
	}

	template<typename T>
	MultidimArrayCuda<T> *transportVectorOfMultidimArrayToGpu(const std::vector<MultidimArrayCuda<T>> &inputVector,
															  cudaStream_t &stream)
	{
		MultidimArrayCuda<T> *outputVectorData;
		transportData(&outputVectorData, inputVector.data(), inputVector.size(), stream);
		return outputVectorData;
	}

	template<typename T>
	T *transportMatrix1DToGpu(const Matrix1D<T> &inputVector, cudaStream_t &stream)
	{
		T *outputVector;
		transportData(&outputVector, inputVector.vdata, inputVector.vdim, stream);
		return outputVector;
	}

	template<typename T>
	T *transportStdVectorToGpu(const std::vector<T> &inputVector, cudaStream_t &stream)
	{
		T *outputVector;
		transportData(&outputVector, inputVector.data(), inputVector.size(), stream);
		return outputVector;
	}

	template<typename T>
	T *transportMatrix2DToGpu(const Matrix2D<T> &inputMatrix, cudaStream_t &stream)
	{
		T *outputMatrixData;
		transportData(&outputMatrixData, inputMatrix.mdata, inputMatrix.mdim, stream);
		return outputMatrixData;
	}

	template<typename T>
	MultidimArrayCuda<T> initializeMultidimArrayCuda(const MultidimArray<T> &multidimArray, cudaStream_t &stream)
	{
		struct MultidimArrayCuda<T> cudaArray = {
			.xdim = multidimArray.xdim, .ydim = multidimArray.ydim, .yxdim = multidimArray.yxdim,
			.xinit = multidimArray.xinit, .yinit = multidimArray.yinit, .zinit = multidimArray.zinit,
			.data = transportMultidimArrayToGpu(multidimArray, stream)
		};

		return cudaArray;
	}

	template<typename T>
	MultidimArrayCuda<T> initializePinnedMultidimArrayCuda(const MultidimArray<T> &multidimArray, cudaStream_t &stream)
	{
		T *pinnedArray;
		if (cudaMallocHost((void **)&pinnedArray,
						   multidimArray.xdim * multidimArray.ydim * multidimArray.zdim * sizeof(T))
			!= cudaSuccess) {
			processCudaError();
		}

		struct MultidimArrayCuda<T> cudaArray = {
			.xdim = multidimArray.xdim, .ydim = multidimArray.ydim, .yxdim = multidimArray.yxdim,
			.xinit = multidimArray.xinit, .yinit = multidimArray.yinit, .zinit = multidimArray.zinit,
		};

		transportData(
			&cudaArray.data, pinnedArray, multidimArray.xdim * multidimArray.ydim * multidimArray.zdim, stream);

		return cudaArray;
	}

	template<typename T>
	void updateMultidimArrayWithGPUData(MultidimArray<T> &multidimArray,
										const MultidimArrayCuda<T> &multidimArrayCuda,
										cudaStream_t &stream)
	{
		transportDataFromGPU(multidimArray.data,
							 multidimArrayCuda.data,
							 multidimArray.xdim * multidimArray.ydim * multidimArray.zdim,
							 stream);
	}

	template<typename T>
	void updateVectorOfMultidimArrayWithGPUData(std::vector<Image<T>> &image,
												const std::vector<MultidimArrayCuda<T>> vectorMultidimArray,
												cudaStream_t &stream)
	{
		assert(image.size() == vectorMultidimArray.size());
		for (int m = 0; m < image.size(); m++) {
			updateMultidimArrayWithGPUData(image[m](), vectorMultidimArray[m], stream);
		}
	}

	template<typename T>
	std::pair<MultidimArrayCuda<T> *, std::vector<MultidimArrayCuda<T>>> convertToMultidimArrayCuda(
		std::vector<Image<T>> &image,
		cudaStream_t &stream)
	{
		std::vector<MultidimArrayCuda<T>> output;
		for (int m = 0; m < image.size(); m++) {
			output.push_back(initializeMultidimArrayCuda(image[m](), stream));
		}
		return std::make_pair(transportVectorOfMultidimArrayToGpu(output, stream), output);
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
		const int RmaxDef,
		cudaStream_t &stream) {
		auto clnm = parameters.clnm;
		auto angles = parameters.angles;

		const size_t idxY0 = usesZernike ? (clnm.size() / 3) : 0;
		const size_t idxZ0 = usesZernike ? (2 * idxY0) : 0;
		const T RmaxF = usesZernike ? RmaxDef : 0;
		const T iRmaxF = usesZernike ? (1.0f / RmaxF) : 0;

		const Matrix2D<T> R = createRotationMatrix<T>(angles);

		struct Program<T>::CommonKernelParameters output = {
			.idxY0 = idxY0, .idxZ0 = idxZ0, .iRmaxF = iRmaxF, .cudaClnm = transportStdVectorToGpu(clnm, stream),
			.cudaR = transportMatrix2DToGpu(R, stream),
		};

		return output;
	}

}  // namespace

template<typename PrecisionType>
Program<PrecisionType>::Program(const Program<PrecisionType>::ConstantParameters parameters)
	: cudaMV(initializePinnedMultidimArrayCuda(parameters.Vrefined(), stream)),
	  VRecMaskF(initializeMultidimArrayCuda(parameters.VRecMaskF, stream)),
	  VRecMaskB(initializeMultidimArrayCuda(parameters.VRecMaskB, stream)),
	  sigma(parameters.sigma),
	  RmaxDef(parameters.RmaxDef),
	  lastX(FINISHINGX(parameters.Vrefined())),
	  lastY(FINISHINGY(parameters.Vrefined())),
	  lastZ(FINISHINGZ(parameters.Vrefined())),
	  loopStep(parameters.loopStep),
	  cudaVL1(transportMatrix1DToGpu(parameters.vL1, stream)),
	  cudaVL2(transportMatrix1DToGpu(parameters.vL2, stream)),
	  cudaVN(transportMatrix1DToGpu(parameters.vN, stream)),
	  cudaVM(transportMatrix1DToGpu(parameters.vM, stream)),
	  blockX(std::__gcd(THREADS_IN_BLOCK, parameters.Vrefined().xdim)),
	  blockY(std::__gcd(THREADS_IN_BLOCK / blockX, parameters.Vrefined().ydim)),
	  blockZ(std::__gcd(THREADS_IN_BLOCK / (blockX * blockY), parameters.Vrefined().zdim)),
	  gridX(parameters.Vrefined().xdim / blockX),
	  gridY(parameters.Vrefined().ydim / blockY),
	  gridZ(parameters.Vrefined().zdim / blockZ)
{
	cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
}

template<typename PrecisionType>
Program<PrecisionType>::~Program()
{
	cudaFree(VRecMaskF.data);
	cudaFree(VRecMaskB.data);
	cudaFree(cudaMV.data);

	cudaFree(const_cast<int *>(cudaVL1));
	cudaFree(const_cast<int *>(cudaVL2));
	cudaFree(const_cast<int *>(cudaVN));
	cudaFree(const_cast<int *>(cudaVM));

	cudaStreamDestroy(stream);
}

template<typename PrecisionType>
template<bool usesZernike>
void Program<PrecisionType>::runForwardKernel(struct DynamicParameters &parameters)

{
	// Unique parameters
	MultidimArrayCuda<PrecisionType> *cudaP, *cudaW;
	std::vector<MultidimArrayCuda<PrecisionType>> pVector, wVector;
	std::tie(cudaP, pVector) = convertToMultidimArrayCuda(parameters.P, stream);
	std::tie(cudaW, wVector) = convertToMultidimArrayCuda(parameters.W, stream);
	auto sigma_size = sigma.size();
	auto cudaSigma = transportStdVectorToGpu(sigma, stream);
	const int step = loopStep;

	// Common parameters
	auto commonParameters = getCommonArgumentsKernel<PrecisionType>(parameters, usesZernike, RmaxDef, stream);

	forwardKernel<PrecisionType, usesZernike>
		<<<dim3(gridX, gridY, gridZ), dim3(blockX, blockY, blockZ), 0, stream>>>(cudaMV,
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

	updateVectorOfMultidimArrayWithGPUData(parameters.P, pVector, stream);
	updateVectorOfMultidimArrayWithGPUData(parameters.W, wVector, stream);

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
	auto cudaMId = initializeMultidimArrayCuda(mId, stream);
	const int step = 1;

	// Common parameters
	auto commonParameters = getCommonArgumentsKernel<PrecisionType>(parameters, usesZernike, RmaxDef, stream);

	backwardKernel<PrecisionType, usesZernike>
		<<<dim3(gridX, gridY, gridZ), dim3(blockX, blockY, blockZ), 0, stream>>>(cudaMV,
																				 cudaMId,
																				 VRecMaskB,
																				 lastZ,
																				 lastY,
																				 lastX,
																				 step,
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

	cudaFree(cudaMId.data);
	freeCommonArgumentsKernel<PrecisionType>(commonParameters);
}

template<typename PrecisionType>
void Program<PrecisionType>::recoverVolumeFromGPU(Image<PrecisionType> &Vrefined)
{
	updateMultidimArrayWithGPUData(Vrefined(), cudaMV, stream);
}

// explicit template instantiation
template class Program<float>;
template class Program<double>;
template void Program<float>::runForwardKernel<true>(struct DynamicParameters &);
template void Program<float>::runForwardKernel<false>(struct DynamicParameters &);
template void Program<double>::runForwardKernel<true>(struct DynamicParameters &);
template void Program<double>::runForwardKernel<false>(struct DynamicParameters &);
template void Program<float>::runBackwardKernel<true>(struct DynamicParameters &);
template void Program<float>::runBackwardKernel<false>(struct DynamicParameters &);
template void Program<double>::runBackwardKernel<true>(struct DynamicParameters &);
template void Program<double>::runBackwardKernel<false>(struct DynamicParameters &);
}  // namespace cuda_forward_art_zernike3D
