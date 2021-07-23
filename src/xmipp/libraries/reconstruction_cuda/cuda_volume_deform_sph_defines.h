#ifndef CUDA_VOLUME_DEFORM_SPH_DEFINES_H
#define CUDA_VOLUME_DEFORM_SPH_DEFINES_H

#define NONE 0
#define MAXWELL 4
#define PASCAL  5
#define TURING  6
#define AMPERE  7

// Universal block sizes
// (cannot be done fully automatic because we are using BLOCK_*_DIM macros in host code)

// Set ARCH to the architecture of your GPU for better performance
#define ARCH PASCAL

#if ARCH == MAXWELL
#define BLOCK_X_DIM 16
#define BLOCK_Y_DIM 4
#define BLOCK_Z_DIM 2
#elif ARCH == PASCAL
#define BLOCK_X_DIM 16
#define BLOCK_Y_DIM 8
#define BLOCK_Z_DIM 1
#elif ARCH == TURING
#define BLOCK_X_DIM 16
#define BLOCK_Y_DIM 8
#define BLOCK_Z_DIM 1
#elif ARCH == AMPERE
#define BLOCK_X_DIM 32
#define BLOCK_Y_DIM 1
#define BLOCK_Z_DIM 4
#else
#define BLOCK_X_DIM 8
#define BLOCK_Y_DIM 4
#define BLOCK_Z_DIM 4
#endif
// Tuning parameters

#endif// CUDA_VOLUME_DEFORM_SPH_DEFINES_H
