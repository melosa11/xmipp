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

#include "classify_projections_2d.h"


void ProgClassifyProjection2D::readParams()
{
	fnImg = getParam("--img");
    sampling = getDoubleParam("--sampling");
	fnOut = getParam("-o");
}


void ProgClassifyProjection2D::defineParams()
{
	addUsageLine("Projections ...");
	addParamsLine("  --img <img_file=\"\">   			: Input Image");
    addParamsLine("  --sampling <sampling=1.0>   			: Sampling rate (A)");
	addParamsLine("  -o <output=\"projection.mrc\">	    : Projections");
}


void ProgClassifyProjection2D::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;

}

template<typename T>
void ProgClassifyProjection2D::projectImage2D(MultidimArray<T> &img)
{
    size_t xdim = XSIZE(img), ydim = YSIZE(img);

    //TODO: check the limits
    size_t semiboxsize = xdim/2;
    size_t lambdalimit = xdim/2;//floor(xdim*sqrt(2)/2);

    std::cout << "boxsize= " << boxsize << "  " << "xdim = " << xdim << "  " << "ydim = " << ydim << std::endl;

    MultidimArray<T> projImg;

    size_t Nprojections = 1;

    //projImg is the radon Tranform. The first row is the projection at 
    // 0 degrees, the second row is the projection at 1 degree and so on
    projImg.initZeros(xdim);

    //TODO: check if the number of pixels is odd or even

    // for (int angle = 0; angle<180; angle++)
    int angle = 0.0;
    double ang= angle*PI/180.0;

    std::vector<double> galleryProjections[Nprojections];

    for (int n = -semiboxsize; n<semiboxsize; n++)
    {
        for (int lambda = -lambdalimit; lambda<lambdalimit; lambda++)
        {
            double c = cos(ang);
            size_t j = round(lambda*c) + semiboxsize;
            size_t i = round(lambda*sin(ang)) + semiboxsize + ((double) n )/c;
            
            std::cout << "i = " << i << "   j = " << j << std::endl;

            if ((j<=xdim) && (i<=ydim))
            {
                A2D_ELEM(projImg, n) += A2D_ELEM(img, i, j);
            }

        }
    }

    galleryProjections[projIdx] = A2D_ELEM(projImg, Nprojections, n);


    Image<double> saveImg;
    saveImg() = projImg;
    saveImg.write("kkk.mrc");
}


void ProgClassifyProjection2D::run()
{
    Image<double> inputImg;
	inputImg.read(fnImg);
    
    std::cout << "starting" << std::endl;
    auto img= inputImg();

    std::cout << "continue" << std::endl;

    projectImage2D(img);



	
}


