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

    // Recycle variable as coounter
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
            std::cout << ";" << std::endl;

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
        std::cout << DIRECT_MULTIDIM_ELEM(projProfile, i) << " ";

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
    size_t vectorSize = XSIZE(xProjection);
    size_t numberReferenceProjections = referenceProjections.size();
    size_t stepReferenceProjections = numberReferenceProjections / 2;

    float angleStepRefenceProjections = 180 / numberReferenceProjections;

    MultidimArray<double> correlationVectorX, correlationVectorY;
    std::vector<double> maximumCorrVector(numberReferenceProjections);
    std::vector<std::vector<size_t>> shiftMaximumCorrVector(numberReferenceProjections);

    for(int n = 0; n < numberReferenceProjections; n++)
    {   
        correlation_vector(xProjection, referenceProjections[n], correlationVectorX);
        correlation_vector(yProjection, referenceProjections[(stepReferenceProjections + n) % numberReferenceProjections], correlationVectorY);

        double maximumCorrX = MINDOUBLE;
        double maximumCorrY = MINDOUBLE;

        size_t shiftX, shiftY;

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
        maximumCorrVector[n] = (maximumCorrX + maximumCorrY);
        std::vector<size_t> shiftsVector {shiftX, shiftY};
        shiftMaximumCorrVector[n] = shiftsVector;
    }
    
    double maxCorr = MINDOUBLE;
    std::vector<size_t> maxShift;
    float angle;

    for(size_t n = 0; n < maximumCorrVector.size(); n++) 
    {
        if(maximumCorrVector[n] > maxCorr)
        {
            maxCorr = maximumCorrVector[n];
            maxShift = shiftMaximumCorrVector[n];
            angle = n * angleStepRefenceProjections;
        }
    }

    size_t shiftX = (maxShift[0] - vectorSize / 2) / 2;
    size_t shiftY = (maxShift[1] - vectorSize / 2) / 2;

    // return angle, shX, shY;
}


void ProgClassifyProjection2D::run()
{
    Image<double> inputImg;
	inputImg.read(fnImg);
    
    std::cout << "starting" << std::endl;
    auto img= inputImg();

    std::cout << "continue" << std::endl;

    std::vector<MultidimArray<double>> normalProjections, imgProjSet;

    imgProjSet = projectImage2D(img, angStep, 0);

    normalProjections = projectImage2D(img, 90, 45);

    MultidimArray<double> xProjection, yProjection;

    xProjection = normalProjections[0];
    yProjection = normalProjections[1];

    correlateProjections(xProjection, yProjection, imgProjSet);
}


