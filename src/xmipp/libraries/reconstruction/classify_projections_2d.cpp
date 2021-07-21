/***************************************************************************
 *
 * Authors:    Federico P. de Isidro Gomez       fdeisidro@cnb.csic.es (2021)
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

#include "classify_projections_2d.h"


void ProgClassifyProjection2D::readParams()
{
	fnImg = getParam("--img");
    // sampling = getDoubleParam("--sampling");
    angStep = getDoubleParam("--angularStep");
	fnOut = getParam("-o");
}


void ProgClassifyProjection2D::defineParams()
{
	addUsageLine("Projections ...");
	addParamsLine("  --img <img_file=\"\">   			: Input Image");
    // addParamsLine("  --sampling <sampling=1.0>   			: Sampling rate (A)");
    addParamsLine("  --angularStep <angStep=1.0>   			: Angular step for projecting (in degrees)");
	addParamsLine("  -o <output=\"projection.mrc\">	    : Projections");
}


void ProgClassifyProjection2D::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;

}

    template<typename T>
std::vector<MultidimArray<double>> ProgClassifyProjection2D::projectImage2D(MultidimArray<T> &img, double angularStep, int angularOffset)
{
    size_t xdim, ydim;

    xdim = XSIZE(img);
    ydim = YSIZE(img);

    size_t boxsize = xdim;

    size_t semiboxsize = boxsize/2;

    // angle_proj = 90:1:269;
    // Nangles = length(angle_proj);
    // projImg = zeros(Nangles, xdim);

    size_t n = (size_t)180/angularStep;

    MultidimArray<double> uniDimProj;
    std::vector<MultidimArray<double>> imgProjSet(n);

    // Recycle variable as counter
    n = 0;

    for (int a = angularOffset; a < 180 + angularOffset; a += angularStep)
    {
        double ang = a * PI / 180.0;
        // std::cout << "angle in rad = " << ang << std::endl;
        double c = cos(ang);
        double s = sin(ang);

        if (xdim % 2 == 0)
        {
            uniDimProj = evenCase(a, c, s, semiboxsize, xdim, ydim, img);
            // std::cout << ";" << std::endl;

            imgProjSet[n] = uniDimProj;
            // projImg(a,:) = uniDimProj;
        }
        // else
        // {
        //     uniDimProj = oddCase(a, c, s, semiboxsize, xdim, ydim, img);
        //     projImg(a,:) = uniDimProj;
        // }

        n++;
    }

    return imgProjSet;
}


MultidimArray<double> ProgClassifyProjection2D::evenCase(double a, double c, double s, int semiboxsize, 
                                        size_t xdim, size_t ydim, MultidimArray<double> &img)
{
    MultidimArray<double> projProfile;
    projProfile.initZeros(xdim);

    int pxidx;

    for (size_t row = 0; row<xdim; row++)
    {
        for (size_t col = 0; col<ydim; col++)
        {
            double p;
            p = sqrt( (row-semiboxsize-0.5)*(row-semiboxsize-0.5) + (col-semiboxsize-0.5)*(col-semiboxsize-0.5));
            
            double cosAngle;
            cosAngle = ( (row-0.5-semiboxsize)*c + (col-0.5-semiboxsize)*s )/p;

            double px = p*cosAngle;

            if (abs(px)<=semiboxsize)
            {
                int signOfpx;
                signOfpx = (-px);

                if (px >= 0)
                {
                    signOfpx = 1;
                }
                else
                {
                    signOfpx = -1;
                }

                pxidx = signOfpx * ceil(abs(px)) + semiboxsize;

                if (signOfpx <=0)
                {
                    pxidx += 1;
                }

                DIRECT_MULTIDIM_ELEM(projProfile, pxidx-1) += DIRECT_A2D_ELEM(img, row, col);
            }
            
        }
    }

    for (size_t i = 0; i<xdim; i++ )
        // std::cout << DIRECT_MULTIDIM_ELEM(projProfile, i) << " ";

    return projProfile;
}


// void ProgClassifyProjection2D::oddCase(double a, double c, double s, 
                                        // int semiboxsize, size_t xdim, size_t ydim, size_t zdim)
// {
    
// }

void ProgClassifyProjection2D::correlateProjections(MultidimArray<double> &xProjection,
                                                    MultidimArray<double> &yProjection,
                                                    std::vector<MultidimArray<double>> &referenceProjections)
{
    std::cout << "Correlatng projections" << std::endl;

    size_t numberReferenceProjections = referenceProjections.size();
    // *** pensar sobre /2
    size_t stepReferenceProjections = numberReferenceProjections / 2;
    
    std::cout << "numberReferenceProjections " << numberReferenceProjections<< std::endl;

    MultidimArray<double> correlationVectorX, correlationVectorY;
    std::vector<double> maximumCorrVector(numberReferenceProjections);
    std::vector<double> maximumCorrVectorX(numberReferenceProjections);
    std::vector<double> maximumCorrVectorY(numberReferenceProjections);
    std::vector<std::vector<int>> shiftMaximumCorrVector(numberReferenceProjections);

    for(int n = 0; n < numberReferenceProjections; n++)
    { 
        #ifdef DEBUG  
        std::cout << "nxyzSize xProj " << NZYXSIZE(xProjection) << std::endl;
        std::cout << "xSize xProj " << XSIZE(xProjection) << std::endl;
        std::cout << "ySize xProj " << YSIZE(xProjection) << std::endl;
        std::cout << "zSize xProj " << ZSIZE(xProjection) << std::endl;

        std::cout << "nxyzSize yProj " << NZYXSIZE(yProjection) << std::endl;
        std::cout << "xSize yProj " << XSIZE(yProjection) << std::endl;
        std::cout << "ySize yProj " << YSIZE(yProjection) << std::endl;
        std::cout << "zSize yProj " << ZSIZE(yProjection) << std::endl;

        std::cout << "xSize referenceProj " << XSIZE(referenceProjections[n]) << std::endl;
        std::cout << "ySize referenceProj " << YSIZE(referenceProjections[n]) << std::endl;
        std::cout << "zSize referenceProj " << ZSIZE(referenceProjections[n]) << std::endl;
        std::cout << "nSize referenceProj " << referenceProjections.size() << std::endl;

        std::cout << "xSize referenceProj " << XSIZE(referenceProjections[(stepReferenceProjections + n) % numberReferenceProjections]) << std::endl;
        std::cout << "ySize referenceProj " << YSIZE(referenceProjections[(stepReferenceProjections + n) % numberReferenceProjections]) << std::endl;
        std::cout << "zSize referenceProj " << ZSIZE(referenceProjections[(stepReferenceProjections + n) % numberReferenceProjections]) << std::endl;
        std::cout << "nSize referenceProj " << referenceProjections.size() << std::endl;
        #endif
        
        correlateProjectionsVectors(xProjection, referenceProjections[n], correlationVectorX);   
        size_t n90 = (stepReferenceProjections + n) % numberReferenceProjections;
        correlateProjectionsVectors(yProjection, referenceProjections[n], correlationVectorY);

        // double correlationIndex(xProjection, referenceProjections[n]);     
        // double correlationIndex(yProjection, referenceProjections[n]);

        #ifdef DEBUG
        std::cout << "ANGLES------------------>" << std::endl;
        std::cout << "xProj angle " << n * angStep << std::endl;
        std::cout << "xProj n " << n << std::endl;
        std::cout << "yProj angle " << n90 * angStep << std::endl;
        std::cout << "yProj n " << n90 << std::endl;
        std::cout << "---------------------------------------------------------------------------" << std::endl;
        #endif

        #ifdef DEBUG
        if(n==90)
        {
            std::cout << "correlationVectorX = [" << std::endl;
            for(size_t i = 0; i < XSIZE(correlationVectorX); i++)
            {
                std::cout << correlationVectorX[i] << " ";
            }
            std::cout << "];" << std::endl;   
            
            std::cout << "correlationVectorY = [" << std::endl;
            for(size_t i = 0; i < XSIZE(correlationVectorY); i++)
            {
                std::cout << correlationVectorY[i] << " ";
            }
            std::cout << "];" << std::endl;
        }
        #endif


        double maximumCorrX = -1e304;
        double maximumCorrY = -1e304;

        int shiftX, shiftY;

        for(size_t i = 0; i < XSIZE(correlationVectorX); i++)
        {
            if(correlationVectorX[i] > maximumCorrX)
            {
                maximumCorrX = correlationVectorX[i];
                shiftX = i;
            }

            if(correlationVectorY[i] > maximumCorrY)
            {
                maximumCorrY = correlationVectorY[i];
                shiftY = i;
            }
        }

        // TODO: check for negative values
        maximumCorrVectorX[n] = maximumCorrX;
        maximumCorrVectorY[n] = maximumCorrY;
        maximumCorrVector[n] = (maximumCorrX + maximumCorrY);
        
        std::vector<int> shiftsVector {shiftX, shiftY};
        shiftMaximumCorrVector[n] = shiftsVector;
    }


    std::cout << "both = [" << std::endl;

    for(size_t i = 0; i < maximumCorrVector.size(); i++)
    {
        std::cout << maximumCorrVector[i] << " ";
    }
    std::cout << "];" << std::endl;


    std::cout << "x = [" << std::endl;

    for(size_t i = 0; i < maximumCorrVectorX.size(); i++)
    {
        std::cout << maximumCorrVectorX[i] << " ";
    }
    std::cout << "];" << std::endl;


    std::cout << "y = [" << std::endl;

    for(size_t i = 0; i < maximumCorrVectorY.size(); i++)
    {
        std::cout << maximumCorrVectorY[i] << " ";
    }
    std::cout << "];" << std::endl;

    
    double maxCorr = MINDOUBLE;
    std::vector<int> maxShift;
    float angle;

    for(size_t n = 0; n < maximumCorrVector.size(); n++) 
    {
        if(maximumCorrVector[n] > maxCorr)
        {
            maxCorr = maximumCorrVector[n];
            maxShift = shiftMaximumCorrVector[n];
            angle = n * angStep;
        }
    }

    size_t vectorSize = XSIZE(xProjection);
    int shiftX = (maxShift[0] - vectorSize / 2) / 2;
    int shiftY = (maxShift[1] - vectorSize / 2) / 2;


    std::cout << MINDOUBLE << std::endl;
    std::cout << "maxCorr " << maxCorr << std::endl;
    std::cout << "angle " << angle << std::endl;
    std::cout << "shiftX " << shiftX << std::endl;
    std::cout << "shiftY " << shiftY << std::endl;

    // return angle, shX, shY;
}

template <typename T>
void ProgClassifyProjection2D::correlateProjectionsVectors(const MultidimArray< T > & m1,
                        const MultidimArray< T > & m2,
                        MultidimArray< double >& R)
{
    // m1.checkDimension(2);
    // m2.checkDimension(2);

    double m1Norm = 0, m2Norm = 0;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(m1)
    {
        m1Norm += m1(i);
    }
    
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(m2)
    {
        m2Norm += m2(i);
    }

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(m1)
    {
        m1(i) /= m1Norm;
    }
    
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(m2)
    {
        m2(i) /= m2Norm;
    }


    // Compute the Fourier Transforms
    MultidimArray< std::complex< double > > FFT1, FFT2;
    FourierTransformer transformer1, transformer2;
    R=m1;
    transformer1.FourierTransform(R, FFT1, false);
    transformer2.FourierTransform((MultidimArray<T> &)m2, FFT2, false);

    // Multiply FFT1 * FFT2'
    double dSize=XSIZE(m1);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(FFT1)
    FFT1(i) *= dSize * conj(FFT2(i));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}



void ProgClassifyProjection2D::run()
{
    Image<double> inputImg;
	inputImg.read(fnImg);
    
    std::cout << "starting" << std::endl;
    auto img= inputImg();

    Image<double> inputTestImage;
    inputTestImage.read("/home/fede/AA_2dClassProy/elipse2d_rot90.mrc");
    auto testImage = inputTestImage();

    std::cout << "continue" << std::endl;

    std::vector<MultidimArray<double>> normalProjections, imgProjSet;

    imgProjSet = projectImage2D(img, angStep, 0);

    normalProjections = projectImage2D(testImage, 90, 45);

    MultidimArray<double> xProjection, yProjection;

    std::cout << normalProjections.size() << std::endl;

    xProjection = normalProjections[0];
    yProjection = normalProjections[1];

    std::cout << "xProj = [" << std::endl;
    for(size_t i = 0; i < XSIZE(xProjection); i++)
    {
        std::cout << xProjection[i] << " ";
    }
    std::cout << "];" << std::endl;   
    
    std::cout << "yProj = [" << std::endl;
    for(size_t i = 0; i < XSIZE(yProjection); i++)
    {
        std::cout << yProjection[i] << " ";
    }
    std::cout << "];" << std::endl;


    std::cout << "reference0 = [" << std::endl;
    for(size_t i = 0; i < XSIZE(imgProjSet[0]); i++)
    {
        std::cout << imgProjSet[0][i] << " ";
    }
    std::cout << "];" << std::endl;   
    
    std::cout << "reference90 = [" << std::endl;
    for(size_t i = 0; i < XSIZE(imgProjSet[90]); i++)
    {
        std::cout << imgProjSet[90][i] << " ";
    }
    std::cout << "];" << std::endl;

    correlateProjections(xProjection, yProjection, imgProjSet);
}


