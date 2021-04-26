/***************************************************************************
 *
 * Authors:    FFederico P. de Isidro Gómez		  fp.deisidro@cnb.csic.es
 * 			   Jose Luis Vilas, 				  jlvilas@cnb.csic.es (2021)
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

#include "resolution_tomo_radial.h"


void ProgResTomoRad::readParams()
{
	fnVol = getParam("--vol");
	fnOut = getParam("-o");
	aroundcenter = checkParam("--aroundCenter");
	thr = getDoubleParam("--thr");

	mask = checkParam("-mask")

	if(mask)
	{
		fnMask = getParam("-mask");
	}
}


void ProgResTomoRad::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file=\"\">                   : Input volume");
	addParamsLine("  -o <output=\"MGresolution.vol\">        : Local resolution volume (in Angstroms)");
	addParamsLine("	 [-mask <vol_mask=\"\"> ] 				 : Mask of regions of interest where resolution values must be considered in the radial profile average");
	addParamsLine("  [--aroundCenter] 					     : Radial average around the center. If this flag is not set, then the radial average is computer around the axis");
	addParamsLine("  [--thr <thr=0.75>]                		 : Threshold (A/px)");
}


void ProgResTomoRad::produceSideInfo()
{

	std::cout << "Starting..." << std::endl;

	Image<double> V;
	V.read(fnVol);

	//V().setXmippOrigin();

	MultidimArray<double> &locresmap=V();
	std::vector<double> resVector;
	
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(locresmap)
		resVector.push_back(DIRECT_MULTIDIM_ELEM(locresmap, n));
	
	std::sort(resVector.begin(),resVector.end());

	double thresholdResolution = resVector[size_t(resVector.size()*thr)];

	std::cout << "resolution threshold = " << thresholdResolution << std::endl;

	size_t xdim, ydim, zdim, ndim;
	locresmap.getDimensions(xdim, ydim, zdim, ndim);

	std::cout << "Resmap dimensions:" << std::endl;
	std::cout << xdim << std::endl;
	std::cout << ydim << std::endl;
	std::cout << zdim << std::endl;
	std::cout << ndim << std::endl;

	xdim = xdim/2;	
	ydim = ydim/2;
	zdim = zdim/2;

	MultidimArray<double> radAvg(xdim), counter(xdim);
	radAvg.initZeros();
	counter.initZeros();

	if(mask){
		Image<int> I;
		I.read(fnMask);
		MultidimArray<int> &maskMap=I();

		size_t xdimM, ydimM, zdimM, ndimM;
		maskMap.getDimensions(xdimM, ydimM, zdimM, ndimM);

		std::cout << "Mask dimensions:" << std::endl;
		std::cout << xdimM << std::endl;
		std::cout << ydimM << std::endl;
		std::cout << zdimM << std::endl;
		std::cout << ndimM << std::endl;
	}

	// MultidimArray<int> maskMap;
	// if(!fnMask.isEmpty())
	// {
	// 	I.read(fnMask);
	// 	MultidimArray<int> &maskMap=I();
	// }



	if (aroundcenter)
	{
		// std::cout << "This is millestone 1" << std::endl;

		FOR_ALL_ELEMENTS_IN_ARRAY3D(locresmap)
		{
			// std::cout << "This is millestone 2" << std::endl;

			double res = A3D_ELEM(locresmap, k, i, j);
			int radius = floor(sqrt((i-ydim)*(i-ydim) + (j-xdim)*(j-xdim) + (k-zdim)*(k-zdim)));

			// std::cout << "This is millestone 3" << std::endl;
			
			int maskValue = A3D_ELEM(maskMap, k, i, j); 
			
			// std::cout << "This is millestone 4" << std::endl;
			// std::cout << maskValue << std::endl;


			if ((res<=thresholdResolution) && (radius<xdim))
			{
				if(mask)
				{
					if(maskValue != 0)
					{
						std::cout << k << std::endl;
						std::cout << i << std::endl;
						std::cout << j << std::endl;
						//std::cout << "i " << i << " j " << j << " k" << k << " " << radius << std::endl;
						DIRECT_MULTIDIM_ELEM(radAvg, radius) +=res;
						DIRECT_MULTIDIM_ELEM(counter, radius) +=1;
					}

				}
				else
				{	
					//std::cout << "i " << i << " j " << j << " k" << k << " " << radius << std::endl;
					DIRECT_MULTIDIM_ELEM(radAvg, radius) +=res;
					DIRECT_MULTIDIM_ELEM(counter, radius) +=1;			

				}
			}
		}
	}
	else
	{
		FOR_ALL_ELEMENTS_IN_ARRAY3D(locresmap)
		{
			double res = A3D_ELEM(locresmap, k, i, j);
			int radius = floor(sqrt((j-xdim)*(j-xdim)));

			if ((res<=thresholdResolution) && (radius<xdim))
			{
				if(mask)
				{
					//std::cout << "i " << i << " j " << j << " k" << k << " " << radius << std::endl;
					DIRECT_MULTIDIM_ELEM(radAvg, radius) +=res;
					DIRECT_MULTIDIM_ELEM(counter, radius) +=1;
				}
				else
				{
					if(A3D_ELEM(maskMap, k, i, j) != 0)
					{
						//std::cout << "i " << i << " j " << j << " k" << k << " " << radius << std::endl;
						DIRECT_MULTIDIM_ELEM(radAvg, radius) +=res;
						DIRECT_MULTIDIM_ELEM(counter, radius) +=1;
					}
				}
			}
		}
	}


	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(radAvg)
		DIRECT_MULTIDIM_ELEM(radAvg, n) /= DIRECT_MULTIDIM_ELEM(counter, n);

	std::cout << radAvg << std::endl;
	
}


void ProgResTomoRad ::testradialAvg()
{
	MultidimArray<double> ball;
	size_t xdim = 1024, ydim= 1440, zdim= 300;
	ball.initZeros(xdim, ydim,zdim);

	xdim = xdim/2;
	ydim = ydim/2;
	zdim = zdim/2;

	FOR_ALL_ELEMENTS_IN_ARRAY3D(ball)
	{
		int radius = floor(sqrt((i-ydim)*(i-ydim) + (j-xdim)*(j-xdim) + (k-zdim)*(k-zdim)));

		A3D_ELEM(ball, k, i, j) = radius;
	}

	FileName fn_ball;
	fn_ball = "mibola.mrc";

	Image<double> imgBall;
	imgBall() = ball;
	imgBall.write(fn_ball);
}


void ProgResTomoRad::run()
{
	produceSideInfo();
	testradialAvg();
}
