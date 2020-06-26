/***************************************************************************
 *
 * Authors:    Amaya Jimenez    (ajimenez@cnb.csic.es)
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

#include <iostream>
#include <random>
#include "particle_polishing_calculate_curve.h"
#include <iterator>
#include <core/xmipp_fft.h>


void ProgParticlePolishingCurve::defineParams()
{
    addUsageLine("Particle polishing from a stack of movie particles");
    addParamsLine(" -i <movie>: Input movie particle metadata");
    //addParamsLine(" -iPart <part>: Input particles metadata");
    addParamsLine(" -vol <volume>: Input volume to generate the reference projections");
    addParamsLine(" --s <samplingRate=1>: Sampling rate");
    addParamsLine(" --nFrames <nFrames>: Number of frames");
    addParamsLine(" --nMics <nMics>: Number of micrographs");
    addParamsLine(" --filter <nFilter=1>: The number of filters to apply");
    addParamsLine(" --movxdim <xmov> : Movie size in x dimension");
    addParamsLine(" --movydim <ymov> : Movie size in y dimension");
    addParamsLine(" [-o <fnOut=\"out.xmd\">]: Output metadata with weighted particles");
    addParamsLine(" [--fixedBW]      : fixed bandwith for the filters. If this flag does not appear, the bandwith will be lower in low frequencies");

}

void ProgParticlePolishingCurve::readParams()
{
	fnMdMov=getParam("-i");
	//fnMdPart=getParam("-iPart");
	fnVol = getParam("-vol");
	fnOut=getParam("-o");
	nFrames=getIntParam("--nFrames");
	nMics=getIntParam("--nMics");
	nFilters=getIntParam("--filter");
	samplingRate=getDoubleParam("--s");
	xmov = getIntParam("--movxdim");
	ymov = getIntParam("--movydim");
        fixedBW  = checkParam("--fixedBW");

}


void ProgParticlePolishingCurve::show()
{
	if (verbose==0)
		return;
	std::cout
	<< "Input movie particle metadata:     " << fnMdMov << std::endl
	<< "Input volume to generate the reference projections:     " << fnVol << std::endl
	;
}



void ProgParticlePolishingCurve::averagingAll(MultidimArray<double> &Iout, const std::vector<double> stks, FileName myfn,
		bool noCurrent, bool applyAlign, const std::vector<double> shiftX, const std::vector<double> shiftY){

	FileName fnPart;
	Image<double> Ipart;
	double count=0.0;

	Iout.initZeros(Ydim, Xdim);
	Iout.setXmippOrigin();
	for (int a=0; a<nFrames; a++){

		int id = (int)stks[a];
		fnPart.compose(id, myfn.removePrefixNumber());

		if (noCurrent){
			std::string fn1 = fnPart.getString().c_str();
			std::string fn2 = myfn.getString().c_str();
			if (fn1.compare(fn2)==0)
				continue;
		}
		Ipart.read(fnPart);
		Ipart().setXmippOrigin();
		if(applyAlign){
			selfTranslate(NEAREST, Ipart(), vectorR2(shiftX[a], shiftY[a]), DONT_WRAP, 0.0);
		}
		Iout+=Ipart();
		count+=1.0;
	}
	Iout/=count;

}



void ProgParticlePolishingCurve::calculateCurve_1(const MultidimArray<double> &Iavg, const MultidimArray<double> &Iproj, MultidimArray<double> &vectorAvg, int nStep, double step, double offset, double Dmin, double Dmax){

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iproj){
		if(DIRECT_MULTIDIM_ELEM(Iproj, n)<Dmin)
			DIRECT_MULTIDIM_ELEM(Iproj, n)=Dmin;
		if(DIRECT_MULTIDIM_ELEM(Iproj, n)>Dmax)
			DIRECT_MULTIDIM_ELEM(Iproj, n)=Dmax;
		int pos = int(floor((DIRECT_MULTIDIM_ELEM(Iproj, n)/step)+offset));
		if(pos>=nStep)
			pos=nStep-1;
		DIRECT_A2D_ELEM(vectorAvg, 0, pos)+=1;
		DIRECT_A2D_ELEM(vectorAvg, 1, pos)+=DIRECT_MULTIDIM_ELEM(Iavg, n);
	}

}


void ProgParticlePolishingCurve::calculateCurve_2(const MultidimArray<double> &Iproj, MultidimArray<double> &vectorAvg, int nStep, double &slope, double &intercept, double Dmin, double Dmax){

	std::vector<double> x, y;
	double aux;
	for(int n=0; n<nStep; n++){
		if (DIRECT_A2D_ELEM(vectorAvg, 0, n)!=0){
			DIRECT_A2D_ELEM(vectorAvg, 1, n)/=DIRECT_A2D_ELEM(vectorAvg, 0, n);
			aux = Dmin + n*((Dmax-Dmin)/(nStep-1));
			x.push_back(double(aux));
			y.push_back(DIRECT_A2D_ELEM(vectorAvg, 1, n));
		}
	}

	double xSum=0, ySum=0, xxSum=0, xySum=0;
	double n = y.size();
	for (int i = 0; i < n; i++){
		xSum += x[i];
		ySum += y[i];
		xxSum += x[i] * x[i];
		xySum += x[i] * y[i];
	}
	slope = (n * xySum - xSum * ySum) / (n * xxSum - xSum * xSum);
	intercept = (ySum - slope * xSum) / n;
}



void ProgParticlePolishingCurve::run()
{

	//TO ESTIMATE THE CURVE FOR EVERY MIC
	MetaData mdPartPrev, mdPart, mdOut;
	MDRow currentRow, outRow;
	mdPartPrev.read(fnMdMov,NULL);
	mdPartSize = mdPartPrev.size();
	mdPart.sort(mdPartPrev, MDL_PARTICLE_ID); //or sort by MDL_MICROGRAPH_ID ¿?

	MDIterator *iterPart = new MDIterator();

	FileName fnPart;

	Image<double> projV, Iavg, I;
	Matrix2D<double> A;
	Projection PV;

	int dataInMovie;
	double slope=0., intercept=0.;
	int nStep=30;
	MultidimArray<double> vectorAvg;

	double Dmin=0., Dmax=0.;
	double stepCurve;
	double offset;

	double rot, tilt, psi, x, y;
	bool flip;
	size_t frId, mvId, partId;
	int enabled;
	std::vector<double> stks;
	CTFDescription ctf;
	size_t objId;

	std::vector<int> mvIdsAux;
	std::vector<double> slopes;
	std::vector<double> intercepts;

	FileName fnAux=fnOut.insertBeforeExtension("_aux");
	fnAux.removeAllExtensions();
	fnAux.addExtension("stk");
	if(fnAux.exists())
		fnAux.deleteFile();

	//PROJECT THE INPUT VOLUME
	V.read(fnVol);
    V().setXmippOrigin();
	int xdimVol = (int)XSIZE(V());
	int ydimVol = (int)YSIZE(V());
	projectorV = new FourierProjector(V(),2,0.5,BSPLINE3);

	for(int m=1; m<=nMics; m++){

		dataInMovie = 0;
		vectorAvg.initZeros(2, nStep);
		slope=0.;
		intercept=0.;

		iterPart->init(mdPart);
		for(int i=0; i<mdPartSize; i++){

			mdPart.getRow(currentRow, iterPart->objId);
			currentRow.getValue(MDL_IMAGE,fnPart);
			if(i==0 && m==1){
				I.read(fnPart, HEADER);
				Ydim=YSIZE(I());
				Xdim=XSIZE(I());
			}
			currentRow.getValue(MDL_OBJID, objId);
			currentRow.getValue(MDL_PARTICLE_ID,partId);
			currentRow.getValue(MDL_MICROGRAPH_ID,mvId);
			currentRow.getValue(MDL_FRAME_ID,frId);
			currentRow.getValue(MDL_ENABLED,enabled);

			if(enabled!=1){
				if(iterPart->hasNext())
					iterPart->moveNext();
				continue;
			}
			if(mvId!=m){
				if(iterPart->hasNext())
					iterPart->moveNext();
				continue;
			}
			dataInMovie++;

			currentRow.getValue(MDL_ANGLE_ROT,rot);
			currentRow.getValue(MDL_ANGLE_TILT,tilt);
			currentRow.getValue(MDL_ANGLE_PSI,psi);
			currentRow.getValue(MDL_SHIFT_X,x);
			currentRow.getValue(MDL_SHIFT_Y,y);
			currentRow.getValue(MDL_FLIP,flip);
			currentRow.getValue(MDL_DM3_VALUE,stks);
			ctf.readFromMdRow(currentRow);
			ctf.produceSideInfo();

			A.initIdentity(3);
			MAT_ELEM(A,0,2)=x;
			MAT_ELEM(A,1,2)=y;
			if (flip)
			{
				MAT_ELEM(A,0,0)*=-1;
				MAT_ELEM(A,0,1)*=-1;
				MAT_ELEM(A,0,2)*=-1;
			}

			if ((xdimVol != Xdim) || (ydimVol != Ydim)){
				std::cerr << "Error: The input particles and volume have different sizes" << std::endl;
				exit(1);
			}

			std::cerr << "Part id " << partId << " " << (int)objId << std::endl;

			projectVolume(*projectorV, PV, xdimVol, xdimVol,  rot, tilt, psi);
			applyGeometry(LINEAR,projV(),PV(),A,IS_INV,DONT_WRAP,0.);
			projV().setXmippOrigin();

			Image<double> aux=PV();
			aux.write(fnAux,objId,true,WRITE_APPEND);
			currentRow.setValue(MDL_POLISHING_PROJ_STK, formatString("%06d@%s", objId, fnAux));

			//filtering the projections with the ctf
			ctf.applyCTF(projV(), samplingRate, false);

			//to obtain the points of the curve (intensity in the projection) vs (counted electrons)
			//the movie particles are averaged (all frames) to compare every pixel value
			bool noCurrent=false;
			bool applyAlign=false;
			std::vector<double> vectorX, vectorY;
			averagingAll(Iavg(), stks, fnPart, noCurrent, applyAlign, vectorX, vectorY);

			/*/DEBUG
			projV.write(formatString("Testprojection_%i_%i.tif", frId, partId));
			Ipart.write(formatString("particle_%i_%i.tif", frId, partId));
			Iavg.write(formatString("Testaverage_%i_%i.tif", frId, partId));
			//END DEBUG*/

			//With Iavg and projV, we calculate the curve (intensity in the projection) vs (counted electrons)
			if(Dmin==0. && Dmax==0.){
				projV().computeDoubleMinMax(Dmin, Dmax);
				Dmin=Dmin-0.2*Dmin;
				Dmax=Dmax+0.2*Dmax;
				stepCurve = (Dmax-Dmin)/double(nStep);
				offset = -Dmin/double(stepCurve);
								//std::cout << Dmin << " " << Dmax << " " << stepCurve << " " << offset << " " << std::endl;
			}
			calculateCurve_1(Iavg(), projV(), vectorAvg, nStep, stepCurve, offset, Dmin, Dmax);

			if(iterPart->hasNext())
				iterPart->moveNext();

		}//end particles loop


		if(dataInMovie>0){
			calculateCurve_2(projV(), vectorAvg, nStep, slope, intercept, Dmin, Dmax);
			//Vectors to store some results
			slopes.push_back(slope);
			intercepts.push_back(intercept);
			mvIdsAux.push_back(m);
			std::cerr << "Estimated curve for movie " << m << ". Slope: " << slope << ". Intercept: " << intercept << std::endl;

		}


	}//end loop mics


	iterPart->init(mdPart);
	for(int i=0; i<mdPartSize; i++){

		//To store the obtained slope and intercept
		mdPart.getRow(currentRow, iterPart->objId);
		currentRow.getValue(MDL_MICROGRAPH_ID,mvId);
		outRow=currentRow;
		for (int hh=0; hh<mvIdsAux.size(); hh++){
			if(mvIdsAux[hh]==mvId){
				outRow.setValue(MDL_POLISHING_SLOPE_CURVE, slopes[hh]);
				outRow.setValue(MDL_POLISHING_INTERCEPT_CURVE, intercepts[hh]);
				break;
			}
		}
		//mdPart.addRow(currentRow);
		mdOut.addRow(outRow);

		if(iterPart->hasNext())
			iterPart->moveNext();
	}//end loop particles

	mdOut.write(fnOut, MD_OVERWRITE);

}


