/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2018)
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
}


void ProgResTomoRad::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file=\"\">                   : Input volume");
	addParamsLine("  -o <output=\"MGresolution.vol\">        : Local resolution volume (in Angstroms)");
	addParamsLine("  [--aroundCenter] 					     : Radial average around the center, if this flag is not set, then the radial average is computer around the axis");
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

	xdim = xdim/2;	
	ydim = ydim/2;
	zdim = zdim/2;

	MultidimArray<double> radAvg(xdim), counter(xdim);
	radAvg.initZeros();
	counter.initZeros();

	if (aroundcenter)
	{
		FOR_ALL_ELEMENTS_IN_ARRAY3D(locresmap)
		{
			double res = A3D_ELEM(locresmap, k, i, j);
			int radius = floor(sqrt((i-ydim)*(i-ydim) + (j-xdim)*(j-xdim) + (k-zdim)*(k-zdim)));

			if ((res<=thresholdResolution) && (radius<xdim))
			{
				//std::cout << "i " << i << " j " << j << " k" << k << " " << radius << std::endl;
				DIRECT_MULTIDIM_ELEM(radAvg, radius) +=res;
				DIRECT_MULTIDIM_ELEM(counter, radius) +=1;
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
				//std::cout << "i " << i << " j " << j << " k" << k << " " << radius << std::endl;
				DIRECT_MULTIDIM_ELEM(radAvg, radius) +=res;
				DIRECT_MULTIDIM_ELEM(counter, radius) +=1;
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
	size_t xdim = 100, ydim= 100, zdim= 100;
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
