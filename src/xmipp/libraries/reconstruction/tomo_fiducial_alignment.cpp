/***************************************************************************
 *
 * Authors:    Federico de Isidro Gomez          fdeisidro@cnb.csic.es (2021)
 * Authors:    Jose Luis Vilas                   jlvilas@cnb.csic.es (2021)
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

#include "tomo_fiducial_alignment.h"
#include "core/xmipp_image_generic.h"


void ProgTomoFidAlign::readParams()
{
	fnImg = getParam("--img");
    fidSize = getDoubleParam("--fiducialSize");
    sampling = getDoubleParam("--sampling");
	fnOut = getParam("-o");
}


void ProgTomoFidAlign::defineParams()
{
	addUsageLine("Projections ...");
	addParamsLine("  --img <img_file=\"\">   			: Input Image");
	addParamsLine("  --fiducialSize <s=10>              : Diameter of the fiducials");
    addParamsLine("  --sampling <s=1>                   : Sampling rate in A");
	addParamsLine("  -o <output=\"projection.mrc\">	    : Projections");
}


void ProgTomoFidAlign::lowPassFilter(
		MultidimArray<std::complex<double>> &fftImg,
		MultidimArray<double> &imgTofilter)
{
	// Filter frequencies
	double highFreqFilt = sampling/fidSize;
	double tail = highFreqFilt + 0.02;

    double lowFreqFilt = sampling/fidSize;

	double idelta = PI/(highFreqFilt-tail);

    double uy, ux, u, uy2;
    size_t ydimImg = YSIZE(imgTofilter);
    size_t xdimImg = XSIZE(imgTofilter);

	long n=0;
	for(size_t i=0; i<YSIZE(fftImg); ++i)
	{
		FFT_IDX2DIGFREQ(i, ydimImg, uy);
		uy2=uy*uy;
		for(size_t j=0; j<XSIZE(fftImg); ++j)
		{
			FFT_IDX2DIGFREQ(j, xdimImg, ux);
			u=sqrt(uy2+ux*ux);
            if (u>=highFreqFilt && u<=lowFreqFilt)
            {
                //double H=0.5*(1+cos((un-w1)*ideltal));
                DIRECT_MULTIDIM_ELEM(fftImg, n) *= 0.5*(1+cos((u-highFreqFilt)*idelta));//H;
            }
            else if (u>tail)
            {
                DIRECT_MULTIDIM_ELEM(fftImg, n) = 0;
            }
			++n;
		}
	}

	FourierTransformer transformer_inv;
	transformer_inv.inverseFourierTransform(fftImg, imgTofilter);

}

void ProgTomoFidAlign::run()
{
    std::cout << "Starting..." << std::endl;
	size_t Xdim, Ydim;
	
	MetaData tiltseriesmd;

    if (fnImg.isMetaData())
    {
        tiltseriesmd.read(fnImg);
    }
    else
    {
        ImageGeneric tiltSeriesImages;
        tiltSeriesImages.read(fnImg, HEADER);

        size_t Zdim, Ndim;
        tiltSeriesImages.getDimensions(Xdim, Ydim, Zdim, Ndim);

        if (fnImg.getExtension() == "mrc" and Ndim == 1)
            Ndim = Zdim;

        size_t id;
        FileName fn;
        for (size_t i = 0; i < Ndim; i++) 
        {
            id = tiltseriesmd.addObject();
            fn.compose(i + FIRST_IMAGE, fnImg);
            tiltseriesmd.setValue(MDL_IMAGE, fn, id);
        }
    }
    tiltseriesmd.write("mtes.xmd");


	long  n = 1;

	FileName fnTSimg;
	size_t objId, objId_ts;
	Image<double> imgTS;

	MultidimArray<double> &ptrImg = imgTS();
    MultidimArray<double> projImgTS;
    MultidimArray<double> filteredImg;
    MultidimArray<std::complex<double>> fftImg;
    MultidimArray<double> freqMap;

	projImgTS.initZeros(Ydim, Xdim);

    FOR_ALL_OBJECTS_IN_METADATA(tiltseriesmd)
	{
		objId = __iter.objId;
		tiltseriesmd.getValue(MDL_IMAGE, fnTSimg, objId);

        std::cout << fnTSimg << std::endl;
        imgTS.read(fnTSimg);

        FourierTransformer transformer1(FFTW_BACKWARD);
        transformer1.FourierTransform(ptrImg, fftImg);

        lowPassFilter(fftImg, ptrImg);//), freqMap);

        for (size_t i =0; i<Ydim;++i)
        {
            for (size_t j =0; j<Xdim;++j)
            {
                DIRECT_A2D_ELEM(projImgTS, i, j) += DIRECT_A2D_ELEM(ptrImg, i, j);
            }
        }
	}

    imgTS() = projImgTS;
	imgTS.write(fnOut);
	
}


