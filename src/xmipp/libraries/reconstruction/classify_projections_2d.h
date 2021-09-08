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

#ifndef _PROG_CLASSIFY_PROJECTIONS_2D
#define _PROG_CLASSIFY_PROJECTIONS_2D

#include <iostream>
#include <core/xmipp_program.h>
#include <core/xmipp_image.h>
#include <core/metadata.h>
#include <core/xmipp_fft.h>
#include <core/xmipp_fftw.h>
#include <math.h>
#include <limits>
#include <complex>
#include <string>


/**@defgroup Classify Image 2D
   @ingroup ReconsLibrary */
//@{


class ProgClassifyProjection2D : public XmippProgram
{
public:
	 /** Filenames */
	FileName fnOut, fnImg, fnImgTest;

	/** sampling rate */
	double sampling, angStep;

public:

    void defineParams();
    void readParams();
    void produceSideInfo();
    
    template<typename T>
    std::vector<MultidimArray<double>> projectImage2D(MultidimArray<T> &img, double angularStep, int angularOffset);
    MultidimArray<double> evenCase(double a, double c, double s, int semiboxsize, 
                                        size_t xdim, size_t ydim, MultidimArray<double> &img);

    void correlateProjections(MultidimArray<double> &xProjection,
                              MultidimArray<double> &yProjection,
                              std::vector<MultidimArray<double>> &referenceProjections);
    template<typename T>
    void correlateProjectionsVectors(const MultidimArray< T > & m1,
                        const MultidimArray< T > & m2,
                        MultidimArray< double >& R);
    void run();

};
//@}
#endif
