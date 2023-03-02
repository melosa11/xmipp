#ifndef CUDA_FORWARD_ART_ZERNIKE3D_H
#define CUDA_FORWARD_ART_ZERNIKE3D_H

#if defined(__CUDACC__)	 // NVCC
#define MY_ALIGN(n) __align__(n)
#elif defined(__GNUC__)	 // GCC
#define MY_ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER)	 // MSVC
#define MY_ALIGN(n) __declspec(align(n))
#else
#error "Please provide a definition for MY_ALIGN macro for your host compiler!"
#endif


// Xmipp includes
#include <core/matrix1d.h>
#include <core/matrix2d.h>
#include <core/multidim_array.h>
#include <core/xmipp_image.h>
// Standard includes
#include <functional>
#include <optional>
#include <vector>

namespace cuda_forward_art_zernike3D {

template<typename T>
struct MY_ALIGN(16) MultidimArrayCuda {
	unsigned xdim;
	unsigned ydim;
	unsigned yxdim;
	int xinit;
	int yinit;
	int zinit;
	T *data;
};

struct Size3D {
	size_t x, y, z;
};

struct VolumeMask {
	std::optional<MultidimArrayCuda<int>> mask;
	Size3D blockSize;
	Size3D gridSize;
};

enum BackwardKernelMode {
	trivial,
	computeMask,
	sharedBlockMask,
};

template<typename PrecisionType = float>
class Program {
	static_assert(std::is_floating_point<PrecisionType>::value, "Floating point type is required.");

   public:
	/// Constant parameters for the computation
	struct ConstantParameters {
		MultidimArray<int> &VRecMaskF, &VRecMaskB;
		Image<PrecisionType> &Vrefined;
		Matrix1D<int> &vL1, &vN, &vL2, &vM;
		std::vector<PrecisionType> &sigma;
		int RmaxDef;
		int loopStep;
		bool computeBackwardMask;
		bool usesZernike;
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
		PrecisionType *cudaClnm;
		Matrix2D<PrecisionType> R;
	};

	struct ZernikeParameters {
		PrecisionType iRmaxF;
		size_t idxY0;
		size_t idxZ0;
		const int *cudaVL1;
		const int *cudaVN;
		const int *cudaVL2;
		const int *cudaVM;
		const PrecisionType *cudaClnm;
	};

	struct CommonBackwardKernelArguments {
		MultidimArrayCuda<PrecisionType> cudaMId;
		PrecisionType r0;
		PrecisionType r1;
		PrecisionType r2;
		PrecisionType r3;
		PrecisionType r4;
		PrecisionType r5;
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
	const BackwardKernelMode backwardMode;

	const VolumeMask backwardMask;

	const MultidimArrayCuda<PrecisionType> cudaMV;

	const int lastX, lastY, lastZ;

	const int RmaxDef;

	const int loopStep;

	const int *cudaVL1, *cudaVN, *cudaVL2, *cudaVM;

	const std::vector<PrecisionType> sigma;

	const int *cudaBlockBackwardMask;

	const PrecisionType *cudaSigma;

	size_t blockXStep, gridXStep;

	int *VRecMaskF;

	unsigned *cudaCoordinatesF;

	const unsigned xdimB, ydimB;

	size_t sizeB;

	const int xdimF, ydimF;

	size_t sizeF;
};

}  // namespace cuda_forward_art_zernike3D
#endif	// CUDA_FORWARD_ART_ZERNIKE3D_H
