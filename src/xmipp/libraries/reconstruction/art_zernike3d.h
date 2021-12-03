/***************************************************************************
 *
 * Authors:    David Herreros Calero dherreros@cnb.csic.es
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


#ifndef _PROG_ART_ZERNIKE3D
#define _PROG_ART_ZERNIKE3D

#include <core/xmipp_metadata_program.h>
#include <core/matrix1d.h>
#include <core/xmipp_image.h>
#include <data/fourier_filter.h>
#include <data/fourier_projection.h>
#include <core/xmipp_error.h>


/** Predict Continuous Parameters. */
class ProgArtZernike3D: public XmippMetadataProgram
{
public:
     /** Filename of the reference volume */
    FileName fnVolR;
    /** Filename of the refined volume */
    FileName fnVolO;
    /// Output directory
    FileName fnOutDir;
    // Metadata with already processed images
    // FileName fnDone;
    /** Degrees of Zernike polynomials and spherical harmonics */
    int L1, L2;
    /** Zernike and SPH coefficients vectors */
    Matrix1D<int> vL1, vN, vL2, vM;
    /** Sampling rate */
    double Ts;
    /** Maximum radius */
    int RmaxDef;
    // Phase Flipped
    bool phaseFlipped;
    // Ignore CTF
    bool ignoreCTF;
    // Regularization ART
    float lambda;
    // Save each # iter
    int save_iter;
    // Correct CTF
    bool useCTF;
    // Apply Zernike
    bool useZernike;
    // Flag for enable/disabled image
    int flagEnabled;

public:
    /** Resume computations */
    bool resume;
    // Number of ART iterations
    int niter;
    // Sort last N projections
    int sort_last_N;
    // 2D and 3D masks in real space
    MultidimArray<int> mask2D;
    // Volume size
    size_t Xdim;
    // Input image
	Image<float> V, Vrefined, Ifilteredp;
    // INput image
    Image<double> I;
    // Spherical mask
    MultidimArray<int> Vmask;
	// Theoretical projection
	Image<float> P;
    // Weight Image
    Image<float> W;
    // Difference Image
    Image<float> Idiff;
    // Transformation matrix
    Matrix2D<float> A;
    // Original angles
    float rot, tilt, psi;
    // Original shift
	float shiftX, shiftY;
	// Original flip
	bool flip;
    // CTF Check
    bool hasCTF;
    // Original defocus
	float defocusU, defocusV, defocusAngle;
	// CTF
	CTFDescription ctf;
    // CTF filter
    FourierFilter FilterCTF;
	// Vector Size
	int vecSize;
	// Vector containing the degree of the spherical harmonics
	Matrix1D<float> clnm;
	// Show optimization
	bool showOptimization;
    // Row ids ordered in a orthogonal fashion
    MultidimArray<size_t> ordered_list;
    // Save iter counter
    int current_save_iter;
    // Image counter
    int num_images;
    // Current ART iteration
    int current_iter;
    // Interpolated voxel
    float voxelI;
    // Volume dimensions
    int initX, endX, initY, endY, initZ, endZ;
public:
    /// Empty constructor
	ProgArtZernike3D();

    /// Destructor
    ~ProgArtZernike3D();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /** Produce side info.
        An exception is thrown if any of the files is not found*/
    void preProcess();

    /** Create the processing working files.
     * The working files are:
     * nmaTodo.xmd for images to process (nmaTodo = mdIn - nmaDone)
     * nmaDone.xmd image already processed (could exists from a previous run)
     */
    // virtual void createWorkFiles();

    /** Predict angles and shift.
        At the input the pose parameters must have an initial guess of the
        parameters. At the output they have the estimated pose.*/
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /// Length of coefficients vector
    void numCoefficients(int l1, int l2, int &vecSize);

    /// Zernike and SPH coefficients allocation
    void fillVectorTerms(int l1, int l2, Matrix1D<int> &vL1, Matrix1D<int> &vN, 
                         Matrix1D<int> &vL2, Matrix1D<int> &vM);

    ///Deform a volumen using Zernike-Spherical harmonic basis
    void deformVol(MultidimArray<float> &mP, MultidimArray<float> &mW,
                   const MultidimArray<float> &mV,
                   float rot, float tilt, float psi);

    // void updateCTFImage(float defocusU, float defocusV, float angle);

    // ART algorithm
    template <int DIRECTION>
    void artModel();

    // Apply Zernike codeformation
    template<bool USESZERNIKE, int DIRECTION>
    void zernikeModel();

    // Interpolation weights + interpolation in 3D
    template<bool INTERPOLATE>
    void weightsInterpolation3D(float x, float y, float z, Matrix1D<float> &w);

    // Remove overdeformation from coefficients
    void removeOverdeformation();

    // virtual void checkPoint();
    
    virtual void finishProcessing();

    virtual void run();

    // Sort images in an orthogonal fashion
    void sortOrthogonal();

};
//@}
#endif