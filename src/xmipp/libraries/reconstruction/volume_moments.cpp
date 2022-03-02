/***************************************************************************
 *
 * Authors:    Oier Lauzirika Zarrabeitia (oierlauzi@bizkaia.eu)
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

#include "core/xmipp_image_generic.h"
#include "volume_moments.h"

// Read arguments ==========================================================
void ProgVolumeMoments::readParams()
{
	// Initialize optional params
	degGeometric = -1;
	degCentralGeometric = -1;

    fnVol = getParam("-i");
    fnOutputRoot = getParam("--oroot");
	if(checkParam("--degGeom"))
		degGeometric = getIntParam("--degGeom");
	if(checkParam("--degCentralGeom"))
		degCentralGeometric = getIntParam("--degCentralGeom");
    if(checkParam("--mask"))
    	mask.readParams(this);
}

// Show ====================================================================
void ProgVolumeMoments::show() const
{
	if (verbose==0)
		return;
    std::cout << "Input volume:  " 						<< fnVol     			<< std::endl;
    std::cout << "Output root:  " 						<< fnOutputRoot			<< std::endl;
    std::cout << "Degree of geometric moments:  " 		<< degGeometric			<< std::endl;
    std::cout << "Degree of central geometric moments: "<< degCentralGeometric	<< std::endl;
    if (mask.type!=NO_MASK)
    	mask.show();
}

// usage ===================================================================
void ProgVolumeMoments::defineParams()
{
	addUsageLine("Compute the specified moments, within a mask (optional).");
    addParamsLine("   -i <volume>             	: Input volume");
    addParamsLine("   -oroot <dir>             	: Output directory");
    addParamsLine("  [--degGeom <N>]            : Degree of geometric moments");
    addParamsLine("  [--degCentralGeom <N>]     : Degree of central geometric moments");
    mask.defineParams(this,INT_MASK);
}

// Produce side information ================================================
void ProgVolumeMoments::produce_side_info()
{
	/*mdVols.read(fnVols);
	mdVols.removeDisabled();

	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(mdVols, Xdim, Ydim, Zdim, Ndim);
	if (mask.type!=NO_MASK)
		mask.generate_mask(Zdim,Ydim,Xdim);
	else
	{
		mask.imask.resizeNoCopy(Zdim,Ydim,Xdim);
		mask.imask.initConstant(1);
	}*/
}

void ProgVolumeMoments::run()
{
    show();
    produce_side_info();

	Image<double> img;
	MultidimArray<double> result;

	// Calculate geometric moments if requested
	if(degGeometric > 0)
	{
		if(mask.type == NO_MASK)
			calcGeometricMoments(img, result);
		else
			calcGeometricMoments(img, mask, result);
	}
	
	// Calculate central geometric moments if requested
	if(degCentralGeometric > 0)
	{
		if(mask.type == NO_MASK)
			calcGeometricMoments(img, result);
		else
			calcGeometricMoments(img, mask, result);

	}

	// TODO
}

void ProgVolumeMoments::calcGeometricMoments(  	const Image<double>& img,
												MultidimArray<double>& result)
{
	double x0, y0, z0;
	calcCenter(img, x0, y0, z0);

	FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(result) {
		DIRECT_NZYX_ELEM(result, l, k, i, j) = calcGeometricMoment(img, j, i, k, x0, y0, z0);
	}
}

void ProgVolumeMoments::calcGeometricMoments(  	const Image<double>& img,
												const Mask& mask,
												MultidimArray<double>& result)
{
	double x0, y0, z0;
	calcCenter(img, x0, y0, z0);

	FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(result) {
		DIRECT_NZYX_ELEM(result, l, k, i, j) = calcGeometricMoment(img, mask, j, i, k, x0, y0, z0);
	}
}

void ProgVolumeMoments::calcCentralGeometricMoments(const Image<double>& img,
													MultidimArray<double>& result)
{
	double x0, y0, z0;
	calcCentroid(img, x0, y0, z0);

	FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(result) {
		DIRECT_NZYX_ELEM(result, l, k, i, j) = calcGeometricMoment(img, j, i, k, x0, y0, z0);
	}
}

void ProgVolumeMoments::calcCentralGeometricMoments(const Image<double>& img,
													const Mask& mask,
													MultidimArray<double>& result)
{
	double x0, y0, z0;
	calcCentroid(img, mask, x0, y0, z0);

	FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(result) {
		DIRECT_NZYX_ELEM(result, l, k, i, j) = calcGeometricMoment(img, mask, j, i, k, x0, y0, z0);
	}
}

double ProgVolumeMoments::calcGeometricMoment(  const Image<double>& img, 
												uint p, uint q, uint r, 
												double x0, double y0, double z0 )
{
	double result = 0;

	FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img()) {
		const auto sample = DIRECT_NZYX_ELEM(img(), l, k, i, j);
		double xmx0 = j - x0;
		double ymy0 = i - y0;
		double zmz0 = k - z0;
		result += sample*pow(xmx0, p)*pow(ymy0, q)*pow(zmz0, r);
	}

	return result;
}

double ProgVolumeMoments::calcGeometricMoment(  const Image<double>& img, 
												const Mask& mask,
												uint p, uint q, uint r, 
												double x0, double y0, double z0 )
{
	double result = 0;

	FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img()) {
		if(DIRECT_NZYX_ELEM(mask.imask, l, k, i, j)) {
			const auto sample = DIRECT_NZYX_ELEM(img(), l, k, i, j);
			double xmx0 = j - x0;
			double ymy0 = i - y0;
			double zmz0 = k - z0;
			result += sample*pow(xmx0, p)*pow(ymy0, q)*pow(zmz0, r);
		}
	}

	return result;
}

void ProgVolumeMoments::calcCenter( const Image<double>& img, 
									double& x0, double& y0, double& z0 )
{
	x0 = (double)(img().xdim) / 2.0;
	y0 = (double)(img().ydim) / 2.0;
	z0 = (double)(img().zdim) / 2.0;
}

void ProgVolumeMoments::calcCentroid(  	const Image<double>& img, 
										double& x0, double& y0, double& z0 )
{
	double xSum = 0.0, ySum = 0.0, zSum = 0.0, sum = 0.0;
	
	FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img()) {
		const auto sample = DIRECT_NZYX_ELEM(img(), l, k, i, j);
		xSum += j*sample;
		ySum += i*sample;
		zSum += k*sample;
		sum += sample;
	}

	x0 = xSum / sum;
	y0 = ySum / sum;
	z0 = zSum / sum;
}

void ProgVolumeMoments::calcCentroid(  	const Image<double>& img, 
										const Mask& mask,
										double& x0, double& y0, double& z0 )
{
	double xSum = 0.0, ySum = 0.0, zSum = 0.0, sum = 0.0;
	
	FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img()) {
		if(DIRECT_NZYX_ELEM(mask.imask, l, k, i, j)) {
			const auto sample = DIRECT_NZYX_ELEM(img(), l, k, i, j);
			xSum += j*sample;
			ySum += i*sample;
			zSum += k*sample;
			sum += sample;
		}
	}

	x0 = xSum / sum;
	y0 = ySum / sum;
	z0 = zSum / sum;
}
