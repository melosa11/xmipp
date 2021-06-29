/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include "reconstruct_art_xray.h"
#include <core/metadata_extension.h>
#include "core/transformations.h"

//ProgReconsXrayART::ProgReconsXrayART()
//{
//
//}
//ProgReconsXrayART::~ProgReconsXrayART()
//{
//
//}
void ProgReconsXrayART::defineParams()
{
    addUsageLine("Generate 3D reconstructions from projections using ART on x-ray tomograms.");
    addUsageLine("+This program reconstructs based on irregular grids given by pseudo atomic structures. ");
    addUsageLine("+In case of regular grids, please refer to other programs as reconstruct_art ");
    addUsageLine("+or reconstruct_fourier. Optionally, a deformation file generated by Nodal Mode ");
    addUsageLine("+Alignment (NMA) can be passed with --nma parameter.");

    //    addSeeAlsoLine("volume_to_pseudoatoms, nma_alignment");

    addParamsLine("   -i <md_file>          : Metadata file with input projections");
    addParamsLine("   [-o <volume_file=\"rec_xray_art.vol\">]  : Filename for output volume.");
    addParamsLine("   --psf <psf_param_file>  : XRay-Microscope parameters file");
    addParamsLine("   [--start <basisvolume_file=\"\">]  : Start from this basis volume. The reconstruction is performed in the same grid as the one ");
    addParamsLine("                                      : in which the basis volume was stored (any -FCC or -CC or grid size value are useless)");
    addParamsLine("  [--sampling_rate <Ts=1>]  : Pixel size (Angstrom)");
    addParamsLine("  alias -s;");
    addParamsLine("  [--threshold <thr=0.05>]  : Normalized Threshold relative to maximum of PSF to reduce the volume into slabs");
    addParamsLine("  [-l <lambda=0.1>]         : Relaxation factor");
    addParamsLine("  [-n <N=1>]                : Number of iterations");
    addParamsLine("  [--thr <N=1>]             : Number of threads to use. NOTE: Not available when using MPI.");

    /// Basis parameters are forced to be voxels
    //    addParamsLine(" == Basis Parameters ==");
    //    Basis::defineParams(this);

}


void ProgReconsXrayART::readParams()
{
    fnDoc = getParam("-i");
    fnOut = getParam("-o");
    fnPSF = getParam("--psf");
    fnStart = getParam("--start");
    sampling = getDoubleParam("--sampling_rate");
    lambdaART = getDoubleParam("-l");
    Nit = getIntParam("-n");
    psfThr = getDoubleParam("--threshold");
    nThreads = getIntParam("--thr");

    // Basis parameters
    //    basis.readParams(this);

    psf.read(fnPSF);

}

void ProgReconsXrayART::preProcess(Image<double> &volVoxels)
{
    MDin.read(fnDoc);

    ImageInfo imgInfo;
    getImageInfo(MDin,imgInfo);
    projXdim = imgInfo.adim.xdim;
    projYdim = imgInfo.adim.ydim;

    psf.calculateParams(sampling*1.e-9, -1, psfThr); // sampling is read in nm
    psf.nThr = nThreads;

    /// Setting initial volume
    if ( fnStart.empty() )
    {
        volVoxels().resize(1, projXdim, projYdim, projXdim, false);
        volVoxels().initZeros();
    }
    else
        volVoxels.read(fnStart);

    volVoxels().setXmippOrigin();

    //Create threads to start working
    thMgr = new ThreadManager(nThreads);

}

double ProgReconsXrayART::singleStep(MultidimArray<double> &muVol, const Projection &projExp,
                                     double rot, double tilt, double psi)
{
    //    muVol.setXmippOrigin();

    size_t iniXdim, iniYdim, iniZdim, newXdim, newYdim, newZdim;

    iniXdim = XSIZE(muVol);
    iniYdim = YSIZE(muVol);
    iniZdim = ZSIZE(muVol);

    Projection projTheo, projNorm;
    MultidimArray<double> Idiff;
    const MultidimArray<double> &mdaProjExp = projExp();
    MultidimArray<double> &mdaProjTheo = projTheo();
    MultidimArray<double> &mdaProjNorm = projNorm();

    mdaProjTheo.initZeros(mdaProjExp);
    projTheo.setAngles(projExp.rot(), projExp.tilt(), projExp.psi());

    mdaProjNorm.initZeros(mdaProjExp);

    newXdim = iniXdim;
    newYdim = iniYdim;
    newZdim = iniZdim;

    MultidimArray<double> rotVol;

    rotVol.resizeNoCopy(newZdim, newYdim, newXdim);
    rotVol.setXmippOrigin();

    // Rotate volume ....................................................

    Matrix2D<double> R;

    Euler_angles2matrix(projExp.rot(), projExp.tilt(), projExp.psi(), R, true);

    double outside = 0; //phantom.iniVol.getPixel(0,0,0,0);

    applyGeometry(1, rotVol, muVol, R, IS_NOT_INV, DONT_WRAP, outside);

    psf.adjustParam(rotVol);

    /// Calculation of IgeoVol to optimize process
    MultidimArray<double> IgeoZb;
    IgeoVol.resize(rotVol, false);
    IgeoZb.resize(1, 1, YSIZE(rotVol), XSIZE(rotVol),false);
    IgeoZb.initConstant(1.);

//    calculateIgeo(rotVol, psf.dzo, IgeoVol, IgeoZb, psf.nThr, thMgr);

    Image<double> imTemp;
    imTemp().alias(IgeoVol);
    imTemp.write("IgeoVol.vol");

    /// Forward projection

    //    projectXrayVolume(rotVol, IgeoVol, psf, projTheo, &mdaProjNorm, thMgr);
    //
    //    projTheo.write("theproj.spi");
    //    projNorm.write("projNorm.spi");

    projectXrayGridVolume(muVol, psf, IgeoVol, projTheo, &mdaProjNorm, FORWARD, thMgr);

    projTheo.write("theproj.spi");
    projNorm.write("projNorm.spi");


    //    projectXrayGridVolume(volBasis, basis, IgeoVol, projTheo, projNorm, FORWARD, nThreads);

    double mean_error = 0;
    Idiff.initZeros(mdaProjExp);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(mdaProjExp)
    {
        // Compute difference image and error
        DIRECT_A2D_ELEM(Idiff, i, j) = DIRECT_A2D_ELEM(mdaProjExp, i, j) - DIRECT_A2D_ELEM(mdaProjTheo, i, j);
        mean_error += DIRECT_A2D_ELEM(Idiff, i, j) * DIRECT_A2D_ELEM(Idiff, i, j);

        // Compute the correction image
        DIRECT_A2D_ELEM(mdaProjNorm, i, j) = XMIPP_MAX(DIRECT_A2D_ELEM(mdaProjNorm, i, j), 1);
        DIRECT_A2D_ELEM(mdaProjNorm, i, j) =
            lambdaART * DIRECT_A2D_ELEM(Idiff, i, j) / DIRECT_A2D_ELEM(mdaProjNorm, i, j);
    }
    mean_error /= YXSIZE(mdaProjExp);

    projNorm.write("projNorm_2.spi");

    imTemp().alias(Idiff);
    imTemp.write("iDiff.spi");

    projectXrayGridVolume(muVol, psf, IgeoVol, projTheo, &mdaProjNorm, BACKWARD, thMgr);

    return mean_error;
}

void ProgReconsXrayART::run()
{
    //    show();

    Image<double> volVoxels;


    preProcess(volVoxels);
    Projection projExp;
    for (int it=0; it<Nit; it++)
    {
        //        preIteration();
        double itError=0;
        FOR_ALL_OBJECTS_IN_METADATA(MDin)
        {
            FileName fnExp;
            MDin.getValue( MDL_IMAGE, fnExp,__iter.objId);
            double rot;
            MDin.getValue( MDL_ANGLE_ROT, rot,__iter.objId);
            double tilt;
            MDin.getValue( MDL_ANGLE_TILT, tilt,__iter.objId);
            double psi;
            MDin.getValue( MDL_ANGLE_PSI, psi,__iter.objId);
            double shiftX;
            MDin.getValue( MDL_SHIFT_X, shiftX,__iter.objId);
            double shiftY;
            MDin.getValue( MDL_SHIFT_Y, shiftY,__iter.objId);

            projExp.read(fnExp);
            projExp().setXmippOrigin();
            projExp.setAngles(rot, tilt, psi);

            itError += singleStep(MULTIDIM_ARRAY(volVoxels), projExp, rot, tilt, psi);
        }
        //        postIteration();

        if (MDin.size()>0)
            itError/=MDin.size();
        std::cerr << "Error at iteration " << it << " = " << itError << std::endl;
    }

    volVoxels.write(fnOut);
}
