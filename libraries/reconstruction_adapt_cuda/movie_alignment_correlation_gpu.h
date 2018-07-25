/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef MOVIE_ALIGNMENT_CORRELATION_GPU
#define MOVIE_ALIGNMENT_CORRELATION_GPU

#include "reconstruction/movie_alignment_correlation_base.h"
#include "reconstruction_cuda/cuda_gpu_movie_alignment_correlation.h"
#include "reconstruction_cuda/cuda_gpu_geo_transformer.h"
#include "data/filters.h"

template<typename T>
class ProgMovieAlignmentCorrelationGPU: public AProgMovieAlignmentCorrelation<T> {
private:
    /**
     * After running this method, all relevant images from the movie are
     * loaded in 'frameFourier' (CPU) and ready for further processing.
     * Scaling and FT is done on GPU
     * @param movie input
     * @param dark correction to be used
     * @param gain correction to be used
     * @param targetOccupancy max frequency to be preserved in FT
     * @param lpf 1D profile of the low-pass filter
     */
    void loadData(const MetaData& movie, const Image<T>& dark,
            const Image<T>& gain, T targetOccupancy,
            const MultidimArray<T>& lpf);

    /**
     * Computes shifts of all images in the 'frameFourier'
     * FT and cross-correlation is computed on GPU
     * @param N number of images to process
     * @param bX pair-wise shifts in X dimension
     * @param bY pair-wise shifts in Y dimension
     * @param A system matrix to be used
     */
    void computeShifts(size_t N, const Matrix1D<T>& bX, const Matrix1D<T>& bY,
            const Matrix2D<T>& A);

    /**
     * This method applies shifts stored in the metadata and computes 'average'
     * image
     * @param movie input
     * @param dark correction to be used
     * @param gain correction to be used
     * @param initialMic sum of the unaligned micrographs
     * @param Ninitial will store number of micrographs used for unaligned sum
     * @param averageMicrograph sum of the aligned micrographs
     * @param N will store number of micrographs used for aligned sum
     */
    void applyShiftsComputeAverage(const MetaData& movie, const Image<T>& dark,
            const Image<T>& gain, Image<T>& initialMic, size_t& Ninitial,
            Image<T>& averageMicrograph, size_t& N);

    /**
     * This method loads images, applies crop, gain and dark correction.
     * @param movie to load
     * @param noOfImgs to load from movie
     * @param dark correction
     * @param gain correction
     * @param cropInput flag stating if the images should be cropped
     * @return pointer to all loaded images
     */
    T* loadToRAM(const MetaData& movie, int noOfImgs, const Image<T>& dark,
            const Image<T>& gain, bool cropInput);

    /**
     * Method loads a single image from the movie
     * @param movie to load from
     * @param objId id of the image to load
     * @param crop flag stating if the image should be cropped
     * @param out loaded frame
     */
    void loadFrame(const MetaData& movie, size_t objId, bool crop,
            Image<T>& out);

    /**
     * Method sets sizes used for processing, e.g. proper sizes of the images,
     * batches etc.
     * @param frame reference frame
     * @param noOfImgs to process
     */
    void setSizes(Image<T> frame, int noOfImgs);

    /**
     * Estimates maximal size of the filter for given frame
     * Might be use to estimate memory requirements
     * @frame reference frame
     * @return max MB necessary for filter
     */
    int getMaxFilterSize(Image<T> &frame);

    void testFFT();
    void testFFTAndScale();
    void testScalingCpuOO();
    void testScalingCpuOE();
    void testScalingCpuEO();
    void testScalingCpuEE();

    void testScalingGpuOO();
    void testScalingGpuOE();
    void testScalingGpuEO();
    void testScalingGpuEE();

    void testFilterAndScale();

private:
    // downscaled Fourier transforms of the input images
    std::complex<T>* frameFourier;

    /**
     * Optimal sizes of the input images, i.e. images loaded from HDD
     * This might differ from actual size of the images, because these
     * values are optimized for further processing
     */
    int inputOptSizeX;
    int inputOptSizeY;
    int inputOptSizeFFTX;
    int inputOptBatchSize;

    /**
     * Optimal sizes of the down-scaled images used for e.g. cross-correlation
     *
     */
    int croppedOptSizeX;
    int croppedOptSizeY;
    int croppedOptSizeFFTX;
    int croppedOptBatchSize;

    /** Memory available for one cross-correlation batch */
    int correlationBufferSizeMB;
    /** No of images in one cross-correlation batch */
    int correlationBufferImgs;
};

#endif /* MOVIE_ALIGNMENT_CORRELATION_GPU */
