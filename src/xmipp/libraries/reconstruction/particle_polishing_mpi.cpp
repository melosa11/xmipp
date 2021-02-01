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

#include "particle_polishing_mpi.h"

#include <iostream>
#include <random>
#include <iterator>
#include <fstream>


void ProgParticlePolishingMpi::defineParams()
{
    each_image_produces_an_output = true;
    allow_apply_geo = false;
    save_metadata_stack = true;
    keep_input_columns = true;
    addUsageLine("Particle polishing from a stack of movie particles");
	XmippMetadataProgram::defineParams();
    addParamsLine(" -vol <volume>: Input volume to generate the reference projections");
    addParamsLine(" --s <samplingRate=1>: Sampling rate");
    addParamsLine(" --nFrames <nFrames>: Number of frames");
    addParamsLine(" --nMics <nMics>: Number of micrographs");
    //addParamsLine(" --filter <nFilter=1>: The number of filters to apply");
    addParamsLine(" --movxdim <xmov> : Movie size in x dimension");
    addParamsLine(" --movydim <ymov> : Movie size in y dimension");
    //addParamsLine(" [--fixedBW]      : Fixed bandwith for the filters. If this flag does not appear, the bandwith will be lower in low frequencies");

}

void ProgParticlePolishingMpi::readParams()
{

	XmippMetadataProgram::readParams();
	fnVol = getParam("-vol");
	nFrames=getIntParam("--nFrames");
	nMics=getIntParam("--nMics");
	//nFilters=getIntParam("--filter");
	samplingRate=getDoubleParam("--s");
	xmov = getIntParam("--movxdim");
	ymov = getIntParam("--movydim");
    //fixedBW  = checkParam("--fixedBW");

}


void ProgParticlePolishingMpi::show()
{
	XmippMetadataProgram::show();
}


void ProgParticlePolishingMpi::similarity (const MultidimArray<double> &I, const MultidimArray<double> &Iexp, double &corrN,
		double &corrM, double &corrW, double &imed, const double &meanF=0.){


	MultidimArray<double> Idiff;
	Idiff=I;
	Idiff-=Iexp;
	double meanD, stdD;
	Idiff.computeAvgStdev(meanD,stdD);
	Idiff.selfABS();
	double thD=stdD;

	/*/DEBUG
	Image<double> Idiff2, I2, Iexp2;
	Idiff2() = Idiff;
	Idiff2.write("diff.mrc");
	I2()=I;
	I2.write("projection.mrc");
	Iexp2()=Iexp;
	Iexp2.write("experimental.mrc");
	//END DEBUG/*/

	double meanI, stdI, thI;
	I.computeAvgStdev(meanI,stdI);
	thI = stdI;

	//std::cerr << "- THI: " << thI << ",  THD: " << thD << std::endl;

	double NI=0;
	double ND=0;
	double sumMI=0, sumMIexp=0;
	double sumI=0, sumIexp=0;
	double sumWI=0, sumWIexp=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Idiff)
	{
		double p=DIRECT_MULTIDIM_ELEM(I,n);
		double pexp=DIRECT_MULTIDIM_ELEM(Iexp,n);
		sumI+=p;
		sumIexp+=pexp;
		if (DIRECT_MULTIDIM_ELEM(Idiff,n)>thD)
		{
			sumWI+=p;
			sumWIexp+=pexp;
			ND+=1.0;
		}
		if(p>thI){
			sumMI+=p;
			sumMIexp+=pexp;
			NI+=1.0;
		}
	}

	double avgI, avgIexp, avgMI, avgMIexp, avgWI, avgWIexp, iNI, iND;
	double isize=1.0/MULTIDIM_SIZE(Idiff);
	avgI=sumI*isize;
	avgIexp=sumIexp*isize;
	if(meanF!=0.){
		//printf("Changing the mean of the experimental images %lf %lf \n", avgIexp, meanF);
		avgIexp=meanF;
		//avgI=100;
	}//else{
		//printf("NO changing the mean of the experimental images %lf \n", avgIexp);
	//}

	if (NI>0){
		iNI=1.0/NI;
		avgMI=sumMI*iNI;
		avgMIexp=sumMIexp*iNI;
	}
	if(ND>0){
		iND=1.0/ND;
		avgWI=sumWI*iND;
		avgWIexp=sumWIexp*iND;
	}

	double sumIIexp=0.0, sumII=0.0, sumIexpIexp=0.0;
	double sumMIIexp=0.0, sumMII=0.0, sumMIexpIexp=0.0;
	double sumWIIexp=0.0, sumWII=0.0, sumWIexpIexp=0.0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Idiff)
	{
		double p=DIRECT_MULTIDIM_ELEM(I,n);
		double pexp=DIRECT_MULTIDIM_ELEM(Iexp,n);
		double pIa=p-avgI;
		double pIexpa=pexp-avgIexp;
		sumIIexp+=(pIa*pIexpa);
		sumII +=(pIa*pIa);
		sumIexpIexp +=(pIexpa*pIexpa);
		if (p>thI){
			pIa=p-avgMI;
			pIexpa=pexp-avgMIexp;
			sumMIIexp+=pIa*pIexpa;
			sumMII +=pIa*pIa;
			sumMIexpIexp +=pIexpa*pIexpa;
		}
		if (DIRECT_MULTIDIM_ELEM(Idiff,n)>thD)
		{
			double w=DIRECT_MULTIDIM_ELEM(Idiff,n);
			pIa=p-avgWI;
			pIexpa=pexp-avgMIexp;
			sumWIIexp+=w*pIa*pIexpa;
			sumWII +=w*pIa*pIa;
			sumWIexpIexp +=w*pIexpa*pIexpa;
		}
	}

	//printf("Some values %lf, %lf, ", sumIIexp, isize);
	sumIIexp*=isize;
	//printf(" %lf, ", sumIIexp);
	sumII*=isize;
	sumIexpIexp*=isize;

	sumMIIexp*=iNI;
	sumMII*=iNI;
	sumMIexpIexp*=iNI;

	sumWIIexp*=iND;
	sumWII*=iND;
	sumWIexpIexp*=iND;

	corrN=sumIIexp/sqrt(sumII*sumIexpIexp);
	corrM=sumMIIexp/sqrt(sumMII*sumMIexpIexp);
	corrW=sumWIIexp/sqrt(sumWII*sumWIexpIexp);
	//corrN=sumIIexp;
	//corrM=sumMIIexp;
	//corrW=sumWIIexp;
	if(std::isnan(corrN))
		corrN=-1.0;
	if(std::isnan(corrM))
		corrM=-1.0;
	if(std::isnan(corrW))
		corrW=-1.0;
	imed=imedDistance(I, Iexp);

}



void ProgParticlePolishingMpi::averagingWindow(MultidimArray<double> &Iout, const std::vector<double> stks, FileName myfn,
		bool applyAlign, const std::vector<double> shiftX, const std::vector<double> shiftY, int Xdim, int Ydim, int window){


	FileName fnPart;
	Image<double> Ipart;
	double count=0.0;
	if(window%2==1)
		window+=1;
	int w=window/2;
	int posi;

	Iout.initZeros(Ydim, Xdim);
	Iout.setXmippOrigin();
	for (int a=0; a<nFrames; a++){
		int id = (int)stks[a];
		fnPart.compose(id, myfn.removePrefixNumber());
		std::string fn1 = fnPart.getString().c_str();
		std::string fn2 = myfn.getString().c_str();
		if (fn1.compare(fn2)==0){
			posi=a;
			break;
		}
	}

	int limI, limS;
	limI=posi-w;
	limS=posi+w;
	//if(limI<0){
	//	limS+=(-limI);
	//	limI=0;
	//}
	//if(limS>nFrames-1){
	//	limI-=(limS-(nFrames-1));
	//	limS=nFrames-1;
	//}
	for (int b=limI; b<=limS; b++){
		if(b>=0 && b<nFrames){
			int id = (int)stks[b];
			fnPart.compose(id, myfn.removePrefixNumber());
			Ipart.read(fnPart);
			Ipart().setXmippOrigin();
			if(applyAlign){
				selfTranslate(NEAREST, Ipart(), vectorR2(shiftX[b], shiftY[b]), DONT_WRAP, 0.0);
			}
			Iout+=Ipart();
			count+=1.0;
		}
	}
	Iout/=count;

}





void ProgParticlePolishingMpi::averagingAll(MultidimArray<double> &Iout, const std::vector<double> stks, FileName myfn,
		bool noCurrent, bool applyAlign, const std::vector<double> shiftX, const std::vector<double> shiftY, int Xdim, int Ydim){

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



void ProgParticlePolishingMpi::calculateCurve_1(const MultidimArray<double> &Iavg, const MultidimArray<double> &Iproj, MultidimArray<double> &vectorAvg, int nStep, double step, double offset, double Dmin, double Dmax){

	//int nStep=30;
	//MultidimArray<double> vectorAvg;
	//vectorAvg.initZeros(2, nStep);
	//double Dmin, Dmax;
	//Iproj.computeDoubleMinMax(Dmin, Dmax);
	//double step = (Dmax-Dmin)/double(nStep);
	//double offset = -Dmin/double(step);
	//std::cerr << "In calculateCurve_1 variables: " << step << ", " << offset << ", " << Dmin << ", " << Dmax <<std::endl;
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
		//if (n==20544){
		//if (pos==0){
		//	std::cerr << n << ", " << pos << ", " << DIRECT_MULTIDIM_ELEM(Iproj, n) << ", " << DIRECT_MULTIDIM_ELEM(Iavg, n) << ", " << DIRECT_A2D_ELEM(vectorAvg, 1, pos) << ", " << DIRECT_A2D_ELEM(vectorAvg, 0, pos) << std::endl;
		//}
	}

}


void ProgParticlePolishingMpi::calculateCurve_2(const MultidimArray<double> &Iproj, MultidimArray<double> &vectorAvg, int nStep, double &slope, double &intercept, double Dmin, double Dmax){

	//std::cerr << "vectorAvg calculado: " << vectorAvg << std::endl;
	//double Dmin, Dmax;
	//Iproj.computeDoubleMinMax(Dmin, Dmax);
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

	/*std::cout << "Vector X " << std::endl;
	for (std::vector<double>::const_iterator it = x.begin(); it != x.end(); ++it)
	{
	    std::cout << *it << std::endl;
	}
	std::cout << "Vector Y " << std::endl;
	for (std::vector<double>::const_iterator it = y.begin(); it != y.end(); ++it)
	{
	    std::cout << *it << std::endl;
	}*/

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


void ProgParticlePolishingMpi::writingOutput(size_t xdim, size_t ydim){

	//MOVIE PARTICLES IMAGES
	MetaData mdPartPrev, mdPart;
	MDRow currentRow;
	mdPartPrev.read(fnMdMov,NULL);
	mdPartSize = mdPartPrev.size();
	mdPart.sort(mdPartPrev, MDL_PARTICLE_ID);

	MDIterator *iterPart = new MDIterator();
	std::vector<int>::iterator it2;

	double rot, tilt, psi, x, y;
	bool flip;
	size_t frId, mvId, partId;
	int xcoor, ycoor;
	int enabled;
	FileName fnPart;
	String aux;

	MetaData SFq2;
	FileName fnRoot2=fnMdMov.insertBeforeExtension("_out_particles_new");
	FileName fnStackOut2=formatString("%s.stk",fnRoot2.c_str());
	if(fnStackOut2.exists())
		fnStackOut2.deleteFile();

	iterPart->init(mdPart);

	int partIdPrev=-1;
	Image<double> Ipart, Ifinal;
	int countForPart=0;
	int countOutMd=0;
	for(int i=0; i<mdPartSize; i++){

		//Project the volume with the parameters in the image
		mdPart.getRow(currentRow, iterPart->objId);
		currentRow.getValue(MDL_PARTICLE_ID,partId);
		currentRow.getValue(MDL_ENABLED,enabled);

		if(enabled==-1){
			if(iterPart->hasNext())
				iterPart->moveNext();
			continue;
		}

		if(partIdPrev!=partId){
			countForPart=0;
			Ifinal().initZeros(ydim, xdim);
			Ifinal().setXmippOrigin();
			partIdPrev=partId;
		}
		currentRow.getValue(MDL_IMAGE,fnPart);
		Ipart.read(fnPart);
		Ipart().setXmippOrigin();

		it2 = find(partIds.begin(), partIds.end(), (int)partId);
		int index = std::distance(partIds.begin(), it2);
		currentRow.getValue(MDL_FRAME_ID,frId);
		//printf("Particle %d frame %d shiftX %d shiftY %d \n", partId, frId, resultShiftX[index+countForPart], resultShiftY[index+countForPart]);

		//AJ TODO
		//selfTranslate(NEAREST, Ipart(), vectorR2((double)resultShiftX[index+frId-1], (double)resultShiftY[index+frId-1]), DONT_WRAP, 0.0);
		Ifinal()+=Ipart();
		countForPart++;

		if (countForPart==nFrames-1){
			Ifinal.write(fnStackOut2,countOutMd,true,WRITE_APPEND);
			countOutMd++;
			FileName fnToSave2;
			fnToSave2.compose(countOutMd, fnStackOut2);
			//size_t id = SFq.addObject();
			//SFq.setValue(MDL_IMAGE, fnToSave, id);
			currentRow.setValue(MDL_IMAGE, fnToSave2);
			SFq2.addRow(currentRow);
		}

		if(iterPart->hasNext())
			iterPart->moveNext();

	} //end mdPartSize loop

	FileName fnOut2;
	fnOut2 = fnMdMov.insertBeforeExtension("_out_new");
	printf("Writing output metadata\n");
	printf("%s \n ", fnOut2.getString().c_str());
	SFq2.write(fnOut2);

}

void ProgParticlePolishingMpi::calculateWeightedFrequency(MultidimArray<double> &Ipart, int Nsteps, const double *frC, const std::vector<double> weights){

	MultidimArray< std::complex<double> > fftIpart;
	MultidimArray<double> wfftIpart;
	FourierTransformer transformer;
	// Auxiliary vector for representing frequency values
	Matrix1D<double> fidx;

	transformer.FourierTransform(Ipart, fftIpart, false);

    //MultidimArray<double> vMag;
    //FFT_magnitude(fftIpart, vMag);
	//Image<double> imageA;
	//imageA()=vMag;
    //imageA.write(formatString("amplitude1.tif"));

	//wfftIpart.initZeros(fftIpart);
	fidx.resizeNoCopy(3);
    for (size_t k=0; k<ZSIZE(fftIpart); k++)
    {
        FFT_IDX2DIGFREQ(k,ZSIZE(Ipart),ZZ(fidx));
        for (size_t i=0; i<YSIZE(fftIpart); i++)
        {
            FFT_IDX2DIGFREQ(i,YSIZE(Ipart),YY(fidx));
            for (size_t j=0; j<XSIZE(fftIpart); j++)
            {
                FFT_IDX2DIGFREQ(j,XSIZE(Ipart),XX(fidx));
                double absw = fidx.module();
                for (int n=0; n<Nsteps; n++){
                	if(absw<=frC[n] && n==0){
                    	//DIRECT_A3D_ELEM(wfftIpart,k,i,j)=weights[0];
                		DIRECT_A3D_ELEM(fftIpart,k,i,j)*= weights[0];
                    	break;
                	}else if((absw<=frC[n] && n!=0)){
                		//DIRECT_A3D_ELEM(wfftIpart,k,i,j)=((weights[n]-weights[n-1])/((n+1)*bandSize - (n)*bandSize)*(absw - (n)*bandSize)) + weights[n-1];
                		DIRECT_A3D_ELEM(fftIpart,k,i,j)*= (((weights[n]-weights[n-1])/(frC[n] - frC[n-1]))*(absw - frC[n-1])) + weights[n-1];
                		break;
                	}else if((absw>=frC[n] && n==Nsteps-1)){
                    	//DIRECT_A3D_ELEM(wfftIpart,k,i,j)=weights[Nsteps-1];
                		DIRECT_A3D_ELEM(fftIpart,k,i,j)*= weights[Nsteps-1];
                    	break;
                	}
                }
            }
        }
    }
	//Image<double> imageW;
	//imageW()=wfftIpart;
    //imageW.write(formatString("pesos.tif"));

    transformer.inverseFourierTransform();

}



void ProgParticlePolishingMpi::calculateWeightedFrequency2(MultidimArray<double> &Ipart, const double aa, const double bb, const double cc){

	MultidimArray< std::complex<double> > fftIpart;
	MultidimArray<double> wfftIpart;
	FourierTransformer transformer;
	// Auxiliary vector for representing frequency values
	Matrix1D<double> fidx;

	transformer.FourierTransform(Ipart, fftIpart, false);

    //MultidimArray<double> vMag;
    //FFT_magnitude(fftIpart, vMag);
	//Image<double> imageA;
	//imageA()=vMag;
    //imageA.write(formatString("amplitude1.tif"));

	//wfftIpart.initZeros(fftIpart);
	fidx.resizeNoCopy(3);
    for (size_t k=0; k<ZSIZE(fftIpart); k++)
    {
        FFT_IDX2DIGFREQ(k,ZSIZE(Ipart),ZZ(fidx));
        for (size_t i=0; i<YSIZE(fftIpart); i++)
        {
            FFT_IDX2DIGFREQ(i,YSIZE(Ipart),YY(fidx));
            for (size_t j=0; j<XSIZE(fftIpart); j++)
            {
                FFT_IDX2DIGFREQ(j,XSIZE(Ipart),XX(fidx));
                double absw = fidx.module();
                DIRECT_A3D_ELEM(fftIpart,k,i,j)*=(aa*exp(-bb*absw)+cc);
            }
        }
    }
	//Image<double> imageW;
	//imageW()=wfftIpart;
    //imageW.write(formatString("pesos.tif"));

    transformer.inverseFourierTransform();

}




void ProgParticlePolishingMpi::calculateWeightedFrequency3(MultidimArray<double> &Ipart, const std::vector<double> aa,
		const std::vector<double> bb, const std::vector<double> cc, const int nFrames, const int posi){

	MultidimArray< std::complex<double> > fftIpart;
	MultidimArray<double> wfftIpart;
	FourierTransformer transformer;
	// Auxiliary vector for representing frequency values
	Matrix1D<double> fidx;

	transformer.FourierTransform(Ipart, fftIpart, false);

    //MultidimArray<double> vMag;
    //FFT_magnitude(fftIpart, vMag);
	//Image<double> imageA;
	//imageA()=vMag;
    //imageA.write(formatString("amplitude1.tif"));

	//wfftIpart.initZeros(fftIpart);
	fidx.resizeNoCopy(3);
    for (size_t k=0; k<ZSIZE(fftIpart); k++)
    {
        FFT_IDX2DIGFREQ(k,ZSIZE(Ipart),ZZ(fidx));
        for (size_t i=0; i<YSIZE(fftIpart); i++)
        {
            FFT_IDX2DIGFREQ(i,YSIZE(Ipart),YY(fidx));
            for (size_t j=0; j<XSIZE(fftIpart); j++)
            {
                FFT_IDX2DIGFREQ(j,XSIZE(Ipart),XX(fidx));
                double absw = fidx.module();
                double myw=0;
                for (int nn=0; nn<nFrames; nn++){
                	myw+=aa[nn]*exp(-bb[nn]*absw)+cc[nn];
                }
                DIRECT_A3D_ELEM(fftIpart,k,i,j)=DIRECT_A3D_ELEM(fftIpart,k,i,j)*(aa[posi]*exp(-bb[posi]*absw)+cc[posi])/myw;
                std::cout << "My w: " << (aa[posi]*exp(-bb[posi]*absw)+cc[posi])/myw << std::endl;
            }
        }
    }
	//Image<double> imageW;
	//imageW()=wfftIpart;
    //imageW.write(formatString("pesos.tif"));

    transformer.inverseFourierTransform();

}



void ProgParticlePolishingMpi::startProcessing()
{

	FileName fnMovieData("moviedata.txt");
	if (fnMovieData.exists())
		std::cout<< "Detected a previous estimation of movie params"<<std::endl;

	//TO ESTIMATE THE CURVE FOR EVERY MIC
	MetaData mdPartPrev, mdPart, mdOut;
	MDRow currentRow, outRow;
	mdPartPrev.read(fn_in,NULL);
	mdPartSize = mdPartPrev.size();
	mdPart.sort(mdPartPrev, MDL_PARTICLE_ID); //or sort by MDL_MICROGRAPH_ID ¿?

	MDIterator *iterPart = new MDIterator();

	FileName fnPart;
	ctfs = new CTFDescription[mdPartSize];
	stkIdxs = new std::vector<double>[mdPartSize];

    Image<double> I;
	Image<double> Ipart, projV, Iavg, Iproj, IpartOut, Iout, Ifinal, IfinalAux, projAux;
	Image<double> Ipartaux, projVaux;
	Matrix2D<double> A, Aout;
	MultidimArray<double> dataArray, softArray;
	Projection PV;

	int dataInMovie;
	double slope=0., intercept=0.;
	int nStep=30;
	MultidimArray<double> vectorAvg;

	double Dmin=0., Dmax=0.;
	double stepCurve;
	double offset;
	int count=0;

	double rot, tilt, psi, x, y;
	bool flip;
	size_t frId, mvId, partId;
	int enabled;
	std::vector<double> stks;
	CTFDescription ctf;
	size_t objId;
	int Xdim, Ydim;
	size_t newId=0;


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

			projectVolume(*projectorV, PV, xdimVol, xdimVol,  rot, tilt, psi);
			applyGeometry(LINEAR,projV(),PV(),A,IS_INV,DONT_WRAP,0.);
			projV().setXmippOrigin();

			//filtering the projections with the ctf
			ctf.applyCTF(projV(), samplingRate, false);

			FileName fnAux;
			newId++;
			fnAux.compose(newId, fnPart.removePrefixNumber());
			//fnAux = fnPart.insertBeforeExtension("proj"); ANTES
			Image<double> aux=projV(); //AJ antes era PV()
			aux.write(fnAux);
		    std::ofstream outproj;
		    outproj.open(fnPart.getString().c_str()); // opens the file
		    if( !outproj ) { // file couldn't be opened
		 	   std::cerr << "Error: file could not be opened" << std::endl;
		 	   exit(1);
		    }
			outproj << fnAux.getString().c_str() << std::endl;
			outproj.close();

			//to obtain the points of the curve (intensity in the projection) vs (counted electrons)
			//the movie particles are averaged (all frames) to compare every pixel value

           //if (fnMovieData.exists()==false){

				bool noCurrent=false;
				bool applyAlign=false;
				std::vector<double> vectorX, vectorY;
				averagingAll(Iavg(), stks, fnPart, noCurrent, applyAlign, vectorX, vectorY, Xdim, Ydim);

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

            //} //fin fnMovieData.exists()==false

			if(iterPart->hasNext())
				iterPart->moveNext();

		}//end particles loop

       //if (fnMovieData.exists()==false){

			if(dataInMovie>0){
				calculateCurve_2(projV(), vectorAvg, nStep, slope, intercept, Dmin, Dmax);
				//Vectors to store some results
				slopes.push_back(slope);
				intercepts.push_back(intercept);
				mvIdsAux.push_back(m);
				std::cerr << "Estimated curve for movie " << m << ". Slope: " << slope << ". Intercept: " << intercept << std::endl;

			}

        //}// end if (fnMovieData.exists()==false){


	}//end loop mics

   //if (fnMovieData.exists()==false){

	   std::ofstream outdata;
	   outdata.open("moviedata.txt"); // opens the file
	   if( !outdata ) { // file couldn't be opened
		   std::cerr << "Error: file could not be opened" << std::endl;
		   exit(1);
	   }
	   for (int a=0; a<mvIdsAux.size(); a++)
		   outdata << mvIdsAux[a] << " " << slopes[a] << " " << intercepts[a] << std::endl;
	   outdata.close();

    //} //end if (fnMovieData.exists()==false){



	/*iterPart->init(mdPart);
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

	mdOut.write(fn_in, MD_OVERWRITE);*/

}



/*
void ProgParticlePolishingMpi::finishProcessing()
{

	//TO ESTIMATE THE CURVE FOR EVERY MIC
	MetaData mdPart;
	MDRow currentRow;
	mdPart.read(fn_in,NULL);
	mdPartSize = mdPart.size();
	MDIterator *iterPart = new MDIterator();
	FileName fnPart, fnAux;

        iterPart->init(mdPart);
	for(int i=0; i<mdPartSize; i++){
		mdPart.getRow(currentRow, iterPart->objId);
		currentRow.getValue(MDL_IMAGE,fnPart);
		fnAux = fnPart.insertBeforeExtension("proj");
		fnAux.deleteFile();
		if(iterPart->hasNext())
			iterPart->moveNext();

	}


}
*/



double expAdjusted_L1(double *x, void *data)
{
	double a=x[1];
	double b=x[2];
	double c=x[3];
	double frC[15]={0.0078, 0.01, 0.0156, 0.0313, 0.05, 0.0625, 0.10, 0.1250, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45};

	double retval=0;
	MultidimArray<double> *auxp = (MultidimArray<double> *) data;
	MultidimArray<double> myData=*(auxp);
	int count=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myData){
		retval+=fabs(DIRECT_MULTIDIM_ELEM(myData,n)-(a*exp(-b*(frC[n]))+c));
		count++;
	}
	return retval;
}



void ProgParticlePolishingMpi::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{

	FileName fnPart, fnAux;
	Image<double> Ipart, projV, Iavg, Iproj, IpartOut, Iout, Ifinal, IfinalAux, projAux;
	Image<double> Ipartaux, projVaux;
	Matrix2D<double> A, Aout;
	MultidimArray<double> dataArray, softArray;
	Projection PV;
	CTFDescription ctf;
	std::vector<double> vectorX, vectorY;

	////////////////////////////////////
	/// FIRST PART - SHIFT ESTIMATION
	////////////////////////////////////

	//ML estimation of shifts
	//-11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
	//-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6
	int shiftX[]={-3, -2, -1, 0, 1, 2, 3};
	int shiftY[]={-3, -2, -1, 0, 1, 2, 3};
    int myL = (int)(sizeof(shiftX)/sizeof(*shiftX));
	MultidimArray<double> lkresults;
	double maxShift = XSIZE(Iavg())/2-10;
	double maxShift2 = maxShift*maxShift;
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	double finalPosX[nFrames];
	double finalPosY[nFrames];

	int countForPart=0;
	int countDisabled=0;
	size_t frId, mvId, partId;
	int enabled;
	double slope=0., intercept=0.;
	std::vector<int>::iterator it;
	std::vector<double> stks;

	//Project the volume with the parameters in the image
	double rot, tilt, psi, xValue, yValue;
	bool flip;
	int xcoor, ycoor;
	int objId;
	int estX;
	int estY;

	rowOut = rowIn;

	rowIn.getValue(MDL_PARTICLE_ID,partId);
	rowIn.getValue(MDL_MICROGRAPH_ID,mvId);
	rowIn.getValue(MDL_FRAME_ID,frId);
	rowIn.getValue(MDL_ENABLED,enabled);

    Image<double> I;
    I.read(fnImg, HEADER);
    int Ydim=YSIZE(I());
    int Xdim=XSIZE(I());

	if(enabled==1){

		rowIn.getValue(MDL_IMAGE,fnAux);
		rowIn.getValue(MDL_ANGLE_ROT,rot);
		rowIn.getValue(MDL_ANGLE_TILT,tilt);
		rowIn.getValue(MDL_ANGLE_PSI,psi);
		rowIn.getValue(MDL_SHIFT_X,xValue);
		rowIn.getValue(MDL_SHIFT_Y,yValue);
		rowIn.getValue(MDL_FLIP,flip);
		rowIn.getValue(MDL_DM3_VALUE,stks);
		//ctf.readFromMdRow(rowIn);
		//ctf.produceSideInfo();

		A.initIdentity(3);
		MAT_ELEM(A,0,2)=xValue;
		MAT_ELEM(A,1,2)=yValue;
		if (flip)
		{
			MAT_ELEM(A,0,0)*=-1;
			MAT_ELEM(A,0,1)*=-1;
			MAT_ELEM(A,0,2)*=-1;
		}

		//creating projection
		//projectVolume(*projectorV, PV, Xdim, Ydim,  rot, tilt, psi);
		//applyGeometry(LINEAR,projVaux(),PV(),A,IS_INV,DONT_WRAP,0.);
		Image<double> projVread;
		//FileName fnAux2=fnAux.insertBeforeExtension("proj"); ANTES
		std::ifstream infile2(fnAux.getString().c_str());
		std::string mystr;
		infile2 >> mystr;
		FileName fnAux2(mystr);
		infile2.close();
		projVaux.read(fnAux2);
		projVaux().setXmippOrigin();

		//AJ puede que esto ya no lo necesite si guardo en proj la proyeccion con el applyGeometry y ctf hechos
		/*projVread.read(fnAux2);
		projVread().setXmippOrigin();
		projVaux().initZeros(Xdim, Ydim);
		applyGeometry(LINEAR,projVaux(),projVread(),A,IS_INV,DONT_WRAP,0.);
		projVaux().setXmippOrigin();
		//apply ctf to projection
		ctf.applyCTF(projVaux(), samplingRate, false);
		*/


		//DEBUG
		//projVaux.write(formatString("myProjFiltered.mrc"));
		//END DEBUG

		//for (int hh=0; hh<mvIdsAux.size(); hh++){
		//	if(mvIdsAux[hh]==mvId){
		//		slope=slopes[hh];
		//		intercept=intercepts[hh];
		//		break;
		//	}
		//}

		std::ifstream infile("moviedata.txt");
		int m;
		double s, i;
		while (infile >> m >> s >> i){
			if(m==mvId){
				slope=s;
				intercept=i;
			}
		}
		infile.close();

		//std::cerr << " FOR Movie ID: " << mvId << ". Working with Slope: " << slope << ". Intercept: " << intercept << std::endl;

		for (int a=0; a<nFrames; a++){

			int id = (int)stks[a];
			fnPart.compose(id, fnImg.removePrefixNumber());

			//std::cerr << fnPart << " " << fnAux2 << std::endl;

			//reading movie particle
			Ipart.read(fnPart);
			Ipart().setXmippOrigin();

			lkresults.initZeros(myL, myL);
			for(int jj=0; jj<myL; jj++){
				for(int hh=0; hh<myL; hh++){

					projV().initZeros(projVaux());
					projV().setXmippOrigin();

					if (shiftX[jj]!=0 || shiftY[hh]!=0){
						translate(LINEAR, projV(), projVaux(), vectorR2((double)shiftX[jj], (double)shiftY[hh]), DONT_WRAP, 0.0); //translation of projV to avoid interpolation problem moving the frame particle
					}else{
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(projVaux()){
							DIRECT_MULTIDIM_ELEM(projV(),n) = DIRECT_MULTIDIM_ELEM(projVaux(),n);
						}
					}
					double likelihood=0.;
					double lambda, fact;
					for(int n=0; n<YSIZE(Ipart()); n++){
						for(int m=0; m<XSIZE(Ipart()); m++){
							if ((n-round(YSIZE(Ipart())/2))*(n-round(YSIZE(Ipart())/2))+(m-round(XSIZE(Ipart())/2))*(m-round(XSIZE(Ipart())/2))>maxShift2) // continue if the Euclidean distance is too far
								continue;
							fact=1.;
							lambda = slope*DIRECT_A2D_ELEM(projV(), n, m)+intercept;
							if (DIRECT_A2D_ELEM(Ipart(), n, m)>0){
								for(int aa=1; aa<=DIRECT_A2D_ELEM(Ipart(), n, m); aa++)
									fact = fact*aa;
							}
							likelihood += -1.*lambda + DIRECT_A2D_ELEM(Ipart(), n, m)*log(lambda) - log(fact);
						}
					}
					DIRECT_A2D_ELEM(lkresults, hh, jj) = likelihood;
				}
			}


			double estXAux, estYAux;
			bestShift(lkresults, estXAux, estYAux, NULL, -1);
			int bestPosX = (int)round(estXAux);
			if(bestPosX>myL-1)
				bestPosX=myL-1;
			else if(bestPosX<0)
				bestPosX=0;
			int bestPosY = (int)round(estYAux);
			if(bestPosY>myL-1)
				bestPosY=myL-1;
			else if(bestPosY<0)
				bestPosY=0;

			//negative because these values are calculated with a displacement in the projection, but then they will be applied to the frames
			finalPosX[a]=-shiftX[bestPosX];
			finalPosY[a]=-shiftY[bestPosY];

			//printf(". BEST POS for particle %d. Shift %lf, %lf \n", a, finalPosX[a], finalPosY[a]);

		} // end loop stks


		PseudoInverseHelper pseudoInverter;
		Matrix2D<double> &X = pseudoInverter.A;
		Matrix1D<double> &y = pseudoInverter.b;
		Matrix1D<double> resultX(nFrames);
		Matrix1D<double> resultY(nFrames);
		X.resizeNoCopy(nFrames,3);
		y.resizeNoCopy(nFrames);
		for(int jj=0; jj<nFrames; jj++){
			X(jj,0)=1;
			X(jj,1)=jj+1;
			X(jj,2)=(jj+1)*(jj+1);
			//X(ii,3)=(ii+1)*(ii+1)*(ii+1);
			y(jj)=finalPosX[jj];
		}
		Matrix1D<double> alpha(3);
		solveLinearSystem(pseudoInverter, alpha);
		//printf("SOLVING LINEAR SYSTEM FOR X \n");
		//printf("alpha(0)=%lf \n", alpha(0));
		//printf("alpha(1)=%lf \n", alpha(1));
		//printf("alpha(2)=%lf \n", alpha(2));
		//printf("alpha(3)=%lf \n", alpha(3));
		matrixOperation_Ax(X, alpha, resultX);

		y.resizeNoCopy(nFrames);
		for(int jj=0; jj<nFrames; jj++){
			y(jj)=finalPosY[jj];
		}
		solveLinearSystem(pseudoInverter, alpha);
		//printf("SOLVING LINEAR SYSTEM FOR Y \n");
		//printf("alpha(0)=%lf \n", alpha(0));
		//printf("alpha(1)=%lf \n", alpha(1));
		//printf("alpha(2)=%lf \n", alpha(2));
		//printf("alpha(3)=%lf \n", alpha(3));
		matrixOperation_Ax(X, alpha, resultY);


		for(int a=0; a<nFrames; a++){
			vectorX.push_back(round(resultX(a)));
			vectorY.push_back(round(resultY(a)));
			//printf("Particle with id %d in frame %d: shiftX %lf shiftY %lf \n", (int)partId, a, resultX(a), resultY(a));
		}
		//resultShiftX.push_back(vectorX);
		//resultShiftY.push_back(vectorY);
		rowOut.setValue(MDL_POLISHING_X, vectorX);
		rowOut.setValue(MDL_POLISHING_Y, vectorY);

		//Writing intermediate output - remove in some moment
		//writingOutput(XSIZE(Ipart()), YSIZE(Ipart()));
		//Writing the ouput
		Ifinal().initZeros(projV());
		Ifinal().setXmippOrigin();
		for(int a=0; a<nFrames; a++){

			int id = (int)stks[a];
			fnPart.compose(id, fnImg.removePrefixNumber());
			Ipartaux.read(fnPart);
			Ipartaux().setXmippOrigin();
			selfTranslate(NEAREST, Ipartaux(), vectorR2(vectorX[a], vectorY[a]), DONT_WRAP, 0.0);
			Ifinal()+=Ipartaux();

		}
		FileName myFnOut = (String)fnImgOut.getString();
		myFnOut = myFnOut.insertBeforeExtension("_shift");
		//printf("Writing output in %s, %s", fnImgOut.getString().c_str(), myFnOut.getString().c_str());
		Ifinal.write(myFnOut);



		////////////////////////////////////
		/// SECOND PART - WEIGHT ESTIMATION
		////////////////////////////////////

		MultidimArray<double> matrixWeights, maxvalues, matrixWeightsPart;
		FourierFilter FilterBP, FilterLP;
		FilterBP.FilterBand=BANDPASS;
		FilterBP.FilterShape=RAISED_COSINE;
		FilterLP.FilterBand=LOWPASS;
		FilterLP.FilterShape=RAISED_COSINE;
		double cutfreq;
		//double bandSize=0.5/(double)nFilters;
		nFilters=15;
		double frC[15]={0.0078, 0.01, 0.0156, 0.0313, 0.05, 0.0625, 0.10, 0.1250, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45};

		//if(fixedBW){
		//	for (int n=0; n<nFilters; n++)
		//		frC[n]=(n+1)*bandSize;
		//}else{
		//	for (int n=nFilters; n>0; n--){
		//		if(n==nFilters)
		//			frC[n-1]=0.5;
		//		else
		//			frC[n-1]=frC[n]/2;
		//	}
		//}

		//printf("FREQ: \n");
		//for (int n=0; n<nFilters; n++)
		//	printf(" %lf ",frC[n]);
		//printf("\n");

		matrixWeightsPart.initZeros(nFrames, nFilters);

		//Creating the projection image
		//projectVolume(*projectorV, PV, Xdim, Ydim,  rot, tilt, psi);
		//applyGeometry(LINEAR,projVaux(),PV(),A,IS_INV,DONT_WRAP,0.);
		/*fnAux2=fnAux.insertBeforeExtension("proj");
		projVread.read(fnAux2);
		projVread().setXmippOrigin();
		projVaux().initZeros(Xdim, Ydim);
		applyGeometry(LINEAR,projVaux(),projVread(),A,IS_INV,DONT_WRAP,0.);
		projVaux().setXmippOrigin();*/

		//To invert contrast in the projections
		double Dmin, Dmax, irange, val;
		projVaux().computeDoubleMinMax(Dmin, Dmax);
		irange=1.0/(Dmax - Dmin);

		for(int a; a<nFrames; a++){

			//Reading the frame
			//double aux=0.;
			int id = (int)stks[a];
			fnPart.compose(id, fnImg.removePrefixNumber());
			//Ipartaux.read(fnPart);
			//Ipartaux().setXmippOrigin();

			for(int n=0; n<nFilters; n++){

				projV().initZeros(projVaux());
				projV().setXmippOrigin();
				//To invert contrast in the projections
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(projVaux()){
					val=DIRECT_MULTIDIM_ELEM(projVaux(),n);
					DIRECT_MULTIDIM_ELEM(projV(),n) = (Dmax - val) * irange;
				}
				//AJ creo que no hace falta porque ya se la aplique antes
				//filtering the projections with the ctf
				//ctf.applyCTF(projV(), samplingRate, false);

				//averaging movie particle image with the ones in all the frames but without the current one (the first true), and applying the align to all of them (last true)
				//bool noCurrent=true;
				//bool applyAlign=true;
				//averagingAll(Iavg(), stks, fnPart, noCurrent, applyAlign, vectorX, vectorY, Xdim, Ydim);

				bool applyAlign=true;
				averagingWindow(Iavg(), stks, fnPart, applyAlign, vectorX, vectorY, Xdim, Ydim, 4);

				//filtering the projected particles with the lowpass filter
				FilterLP.w1=frC[n]; //(n+1)*bandSize;
				FilterLP.generateMask(projV());
				FilterLP.applyMaskSpace(projV());

				//filtering the averaged movie particle (leaving out the j-frame)
				FilterLP.generateMask(Iavg());
				FilterLP.applyMaskSpace(Iavg());

				//DEBUG
				//Iavg.write(formatString("myAvgFiltered.mrc"));
				//projV.write(formatString("myProjFiltered.mrc"));
				//END DEBUG

				//calculate similarity measures between averaged movie particles and filtered projection
				double corrN, corrM, corrW, imed;
				similarity(projV(), Iavg(), corrN, corrM, corrW, imed); //with or without meanD (as last parameter) ¿?
				DIRECT_A2D_ELEM(matrixWeightsPart, a, n) = corrN;
				//aux += DIRECT_A2D_ELEM(matrixWeightsPart, a, n);

			} //end frequencies loop

		} //end loop frames


		std::vector<double> vectorAux;
		FOR_ALL_ELEMENTS_IN_ARRAY2D(matrixWeightsPart)
			vectorAux.push_back(A2D_ELEM(matrixWeightsPart,i,j));
		rowOut.setValue(MDL_POLISHING_FREQ_COEFFS_BEFORE, vectorAux);

		//To normalize the weights
		//for(int nn=0; nn<nFilters; nn++){
		//	double aux=0;
		//	for(int mm=0; mm<nFrames; mm++)
		//		aux += DIRECT_A2D_ELEM(matrixWeightsPart, mm, nn);
		//	for(int mm=0; mm<nFrames; mm++)
		//		DIRECT_A2D_ELEM(matrixWeightsPart, mm, nn) /= aux;
		//}

		std::vector<double> vectorAux2;
		FOR_ALL_ELEMENTS_IN_ARRAY2D(matrixWeightsPart)
			vectorAux2.push_back(A2D_ELEM(matrixWeightsPart,i,j));
		rowOut.setValue(MDL_POLISHING_FREQ_COEFFS_AFTER_NORM, vectorAux2);

		//Adjust to a negative exponential function
		MultidimArray<double> data;
		std::vector <double> aa,bb,cc;
		printf("Data to fit \n");
		for(int mm=0; mm<nFrames; mm++){
			data.initZeros(nFilters);
			for(int nn=0; nn<nFilters; nn++){
				DIRECT_A1D_ELEM(data, nn) = DIRECT_A2D_ELEM(matrixWeightsPart, mm, nn);
				printf("%lf ", DIRECT_A1D_ELEM(data, nn));
			}
			printf("\n");
		    Matrix1D<double> p(3), steps(3);
		    p(0)=1; // a in a*exp(-b*x)+c
		    p(1)=1; // b in a*exp(-b*x)+c
		    p(2)=0; // c in a*exp(-b*x)+c
		    steps.initConstant(1);
		    double cost;
		    int iter;
			powellOptimizer(p, 1, 3, &expAdjusted_L1, &data, 0.001, cost, iter, steps, false);
			aa.push_back(p(0));
			bb.push_back(p(1));
			cc.push_back(p(2));
			printf("Obtained params for frame %d: %lf, %lf, %lf with cost %lf and iters %d \n", mm, p(0), p(1), p(2), cost, iter);
		}

		std::vector<double> vectorAux3;
		for(int a=0; a<nFrames; a++){
			vectorAux3.push_back(aa[a]);
			vectorAux3.push_back(bb[a]);
			vectorAux3.push_back(cc[a]);
		}
		rowOut.setValue(MDL_POLISHING_FREQ_COEFFS_AFTER_SVD, vectorAux3);

		//Writing the ouput
		Ifinal().initZeros(projV());
		Ifinal().setXmippOrigin();
		for(int a=0; a<nFrames; a++){

			int id = (int)stks[a];
			fnPart.compose(id, fnImg.removePrefixNumber());
			Ipartaux.read(fnPart);
			Ipartaux().setXmippOrigin();

			selfTranslate(NEAREST, Ipartaux(), vectorR2(vectorX[a], vectorY[a]), DONT_WRAP, 0.0);

			//calculateWeightedFrequency2(Ipartaux(), aa[a], bb[a], cc[a]);
			calculateWeightedFrequency3(Ipartaux(), aa, bb, cc, nFrames, a);
			Ifinal()+=Ipartaux();

		}

		rowOut.setValue(MDL_IMAGE, fnImgOut);
		Ifinal.write(fnImgOut);
		printf("Particle %s finished \n", fnImgOut.getString().c_str());



		/*
		std::vector<double> vectorAux;
		FOR_ALL_ELEMENTS_IN_ARRAY2D(matrixWeightsPart)
			vectorAux.push_back(A2D_ELEM(matrixWeightsPart,i,j));
		rowOut.setValue(MDL_POLISHING_FREQ_COEFFS_BEFORE, vectorAux);

		//to normalize the weights
		for(int n=0; n<nFilters; n++){
			double aux=0, aux2=0;
			for(int mm=0; mm<nFrames; mm++){
				aux += DIRECT_A2D_ELEM(matrixWeightsPart, mm, n);
			}
			for(int mm=0; mm<nFrames; mm++){
				DIRECT_A2D_ELEM(matrixWeightsPart, mm, n) = (2*(aux/nFrames) - DIRECT_A2D_ELEM(matrixWeightsPart, mm, n)); //aux;
				if (DIRECT_A2D_ELEM(matrixWeightsPart, mm, n)<0)
					DIRECT_A2D_ELEM(matrixWeightsPart, mm, n)=0.;
				aux2 += DIRECT_A2D_ELEM(matrixWeightsPart, mm, n);
			}
			for(int mm=0; mm<nFrames; mm++){
				if (aux2!=0)
					DIRECT_A2D_ELEM(matrixWeightsPart, mm, n) /= aux2;
				else
					DIRECT_A2D_ELEM(matrixWeightsPart, mm, n) = 0.;
			}
		}

		std::vector<double> vectorAux2;
		FOR_ALL_ELEMENTS_IN_ARRAY2D(matrixWeightsPart)
			vectorAux2.push_back(A2D_ELEM(matrixWeightsPart,i,j));
		rowOut.setValue(MDL_POLISHING_FREQ_COEFFS_AFTER_NORM, vectorAux2);

		Matrix2D<double> matrixWeightsMat;
		matrixWeightsPart.copy(matrixWeightsMat);
		Matrix2D<double> U,V;
		Matrix1D<double> S;
		svdcmp(matrixWeightsMat,U,S,V);
		Matrix2D<double> Smat;
		Smat.initZeros(S.size(),S.size());
		for (int h=0; h<S.size(); h++){
			if (h<2)
				Smat(h,h)=S(h);
		}
		Matrix2D<double> result1;
		Matrix2D<double> resultWeights;
		matrixOperation_AB(U, Smat, result1);
		matrixOperation_ABt(result1, V, resultWeights);
		//std::cerr << "For particle " << fnImg << std::endl;
		//std::cerr << "- SVD recons: " << resultWeights << std::endl;

		std::vector<double> vectorAux3;
		for(int a=0; a<nFrames; a++)
			for(int nn=0; nn<nFilters; nn++)
				vectorAux3.push_back(resultWeights(a,nn));
		rowOut.setValue(MDL_POLISHING_FREQ_COEFFS_AFTER_SVD, vectorAux3);
		*/

		/*
		//Writing the ouput
		FileName fnTest;
		Ifinal().initZeros(projV());
		Ifinal().setXmippOrigin();
		for(int a=0; a<nFrames; a++){

			std::vector<double> myweights;
			int id = (int)stks[a];
			fnPart.compose(id, fnImg.removePrefixNumber());
			Ipartaux.read(fnPart);
			Ipartaux().setXmippOrigin();

			selfTranslate(NEAREST, Ipartaux(), vectorR2(vectorX[a], vectorY[a]), DONT_WRAP, 0.0);

			for(int nn=0; nn<nFilters; nn++)
				myweights.push_back(resultWeights(a, nn));

			calculateWeightedFrequency(Ipartaux(), nFilters, frC, myweights);
			Ifinal()+=Ipartaux();

		}

		//FileName myFnOut = (String)fnImgOut.getString();
		//myFnOut = myFnOut.removeAllExtensions();
		//myFnOut = myFnOut.addExtension("stk");
		//printf("Writing output in %s, %s", fnImgOut.getString().c_str(), myFnOut.getString().c_str());
		Ifinal.write(fnImgOut);
		printf("Particle %s finished \n", fnImgOut.getString().c_str());
		*/

	} //end if enabled

}



