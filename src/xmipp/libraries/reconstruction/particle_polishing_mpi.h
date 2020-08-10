/***************************************************************************
 *
 * Authors: Amaya Jimenez      (ajimenez@cnb.csic.es)
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

#ifndef LIBRARIES_RECONSTRUCTION_PARTICLE_POLISHING_MPI_H_
#define LIBRARIES_RECONSTRUCTION_PARTICLE_POLISHING_MPI_H_

#include <core/xmipp_program.h>
#include <core/xmipp_image.h>
#include <core/metadata.h>
#include <data/fourier_projection.h>
#include <data/fourier_filter.h>
#include <data/filters.h>
#include <data/numerical_tools.h>


class ProgParticlePolishingMpi: public XmippMetadataProgram
{
protected:
	FileName fnMdMov; // File with the input movie particles metadata
	FileName fnMdPart; // File with the input particles metadata
	FileName fnVol;
	FileName fnOut;
    MetaData mdPart;
	Image<double> V;
	FourierProjector *projectorV;
	int nFrames, nMics, nFilters, xmov, ymov;
	double samplingRate;

	std::vector<int> frIds, mvIds, partIds, xcoors, ycoors, enables, objIds;
	std::vector<double> rots, tilts, psis, xs, ys;
	std::vector<bool> flips;
	std::vector<std::string> fnParts;
	CTFDescription *ctfs;
	std::vector<double> *stkIdxs;
	size_t mdPartSize;
	std::vector<std::vector<double>> resultShiftX;
	std::vector<std::vector<double>> resultShiftY;
	bool fixedBW;
	std::vector<int> mvIdsAux;
	std::vector<double> slopes;
	std::vector<double> intercepts;
	std::vector<int> processParts;


protected:

    void defineParams();
    void readParams();
    void show();
    void startProcessing();
    //void preProcess();
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    void similarity (const MultidimArray<double> &I1, const MultidimArray<double> &I2, double &corrN, double &corrM, double &corrW, double &imed, const double &meanF);
    void averagingAll(MultidimArray<double> &Iout, const std::vector<double> stks, FileName myfn, bool noCurrent, bool applyAlign, const std::vector<double> shiftX, const std::vector<double> shiftY, int Xdim, int Ydim);
    void calculateCurve_1(const MultidimArray<double> &Iavg, const MultidimArray<double> &Iproj, MultidimArray<double> &vectorAvg, int nStep, double step, double offset, double Dmin, double Dmax);
    void calculateCurve_2(const MultidimArray<double> &Iproj, MultidimArray<double> &vectorAvg, int nStep, double &slope, double &intercept, double Dmin, double Dmax);
    void writingOutput(size_t xdim, size_t ydim);
    void calculateWeightedFrequency(MultidimArray<double> &Ipart, int Nsteps, const double *frC, const std::vector<double> weights);
    void averagingWindow(MultidimArray<double> &Iout, const std::vector<double> stks, FileName myfn,
    		bool applyAlign, const std::vector<double> shiftX, const std::vector<double> shiftY, int Xdim, int Ydim, int window);

};

#endif /* LIBRARIES_RECONSTRUCTION_PARTICLE_POLISHING_MPI_H_ */
