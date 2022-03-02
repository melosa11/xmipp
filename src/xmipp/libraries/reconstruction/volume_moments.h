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
#ifndef _PROG_VOLUME_MOMENTS
#define _PROG_VOLUME_MOMENTS

#include "core/xmipp_image.h"
#include "core/xmipp_program.h"
#include "data/mask.h"

///@defgroup VolumeMoments Volume Moments
///@ingroup ReconsLibrary
//@{
/** Volume Moments parameters. */
class ProgVolumeMoments: public XmippProgram
{
public:
    /// Input volume
    FileName fnVol;
    /// Input mask
    Mask mask;
    /// Output root dir
    FileName fnOutputRoot;
    /// Maximum degree of geometric moments
    int degGeometric;
    /// Maximum degree of central geometric moments
    int degCentralGeometric;

public:
public:
    /// Read arguments
    void readParams();

    /// Show
    void show() const;

    /// Define parameters
    void defineParams();

    /** Produce side info.*/
    void produce_side_info();

    /** Run */
    void run();

    /** Calculate geometric moments */
    static void calcGeometricMoments(   const Image<double>& img,
                                        MultidimArray<double>& result);

    /** Calculate geometric moments with mask */
    static void calcGeometricMoments(   const Image<double>& img,
                                        const Mask& mask,
                                        MultidimArray<double>& result);

    /** Calculate central geometric moments */
    static void calcCentralGeometricMoments(const Image<double>& img,
                                            MultidimArray<double>& result);

    /** Calculate central geometric moments with mask */
    static void calcCentralGeometricMoments(const Image<double>& img,
                                            const Mask& mask,
                                            MultidimArray<double>& result);

    /** Calculate central geometric moment */
    static double calcGeometricMoment(  const Image<double>& img, 
                                        uint p, uint q, uint r, 
                                        double x0, double y0, double z0 );
    
    /** Calculate central geometric moment with mask*/
    static double calcGeometricMoment(  const Image<double>& img, 
                                        const Mask& mask,
                                        uint p, uint q, uint r, 
                                        double x0, double y0, double z0 );

    /** Calculate centre */
    static void calcCenter( const Image<double>& img, 
                            double& x0, double& y0, double& z0 );

    /** Calculate centroid */
    static void calcCentroid(   const Image<double>& img, 
                                double& x0, double& y0, double& z0 );
    
    /** Calculate centroid with mask */
    static void calcCentroid(   const Image<double>& img, 
                                const Mask& mask,
                                double& x0, double& y0, double& z0 );
};
//@}
#endif
