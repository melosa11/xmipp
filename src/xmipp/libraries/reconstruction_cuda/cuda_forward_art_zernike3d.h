#ifndef CUDA_FORWARD_ART_ZERNIKE3D_H
#define CUDA_FORWARD_ART_ZERNIKE3D_H

// Xmipp includes
#include <core/matrix1d.h>
#include <core/matrix2d.h>
#include <core/multidim_array.h>
#include <core/xmipp_image.h>
// Standard includes
#include <vector>

namespace cuda_forward_art_zernike3D {

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
	T *data;
};

struct BlockSizes {
	size_t x, y, z;
};

template<typename PrecisionType = float>
class Program {
	static_assert(std::is_floating_point<PrecisionType>::value, "Floating point type is required.");

	using PrecisionType3 = std::conditional<std::is_same<PrecisionType, float>::value, float3, double3>;

   public:
	/// Constant parameters for the computation
	struct ConstantParameters {
		MultidimArray<int> &VRecMaskF, &VRecMaskB;
		Image<PrecisionType> &Vrefined;
		Matrix1D<int> &vL1, &vN, &vL2, &vM;
		std::vector<PrecisionType> &sigma;
		int RmaxDef;
		int loopStep;
	};

	struct AngleParameters {
		PrecisionType rot, tilt, psi;
	};

	struct DynamicParameters {
		const std::vector<PrecisionType> &clnm;
		std::vector<Image<PrecisionType>> &P;
		std::vector<Image<PrecisionType>> &W;
		const Image<PrecisionType> &Idiff;
		struct AngleParameters angles;
	};

	struct CommonKernelParameters {
		size_t idxY0, idxZ0;
		PrecisionType iRmaxF;
		PrecisionType *cudaClnm, *cudaR;
	};

   public:
	template<bool usesZernike>
	void runForwardKernel(struct DynamicParameters &parameters);

	template<bool usesZernike>
	void runBackwardKernel(struct DynamicParameters &parameters);

	/// Moves Volume from GPU and writes it to Vrefined
	/// IMPORTANT: Memory heavy operation.
	void recoverVolumeFromGPU(Image<PrecisionType> &Vrefined);

	explicit Program(const ConstantParameters parameters);
	~Program();

   private:
	const MultidimArrayCuda<int> VRecMaskF, VRecMaskB;

	const MultidimArrayCuda<PrecisionType> cudaMV;

	const int lastX, lastY, lastZ;

	const int RmaxDef;

	const int loopStep;

	const int *cudaVL1, *cudaVN, *cudaVL2, *cudaVM;

	const std::vector<PrecisionType> sigma;

	const size_t blockX, blockY, blockZ, gridX, gridY, gridZ;

	const size_t blockXStep, blockYStep, blockZStep, gridXStep, gridYStep, gridZStep;
};

}  // namespace cuda_forward_art_zernike3D
#endif	// CUDA_FORWARD_ART_ZERNIKE3D_H
