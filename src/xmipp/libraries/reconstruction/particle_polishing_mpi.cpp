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
    addParamsLine(" --movxdim <xmov> : Movie size in x dimension");
    addParamsLine(" --movydim <ymov> : Movie size in y dimension");

}

void ProgParticlePolishingMpi::readParams()
{

	XmippMetadataProgram::readParams();
	fnVol = getParam("-vol");
	nFrames=getIntParam("--nFrames");
	nMics=getIntParam("--nMics");
	samplingRate=getDoubleParam("--s");
	xmov = getIntParam("--movxdim");
	ymov = getIntParam("--movydim");

}


void ProgParticlePolishingMpi::show()
{
	XmippMetadataProgram::show();
}


void ProgParticlePolishingMpi::similarity (const MultidimArray<double> &I, const MultidimArray<double> &Iexp, double &corrN,
		double &corrM, double &corrW, double &imed, const double &meanF=0.){

	//en esta funcion se calculan varios tipos de correlaciones, por si alguna medida mas fuera util
	//pero ahora mismo solo se usa una correlacion normal y creo que es la unica que tiene sentido

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
		avgIexp=meanF;
	}

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

	sumIIexp*=isize;
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

	//En esta funcion se recorren los frames de una particula para promediarlos, pero solo aquellos que esten dentro de la ventana considerada
	//se puede hacer con o sin shifts utilizando applyAlign, shiftX y shiftY

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

	//En esta funcion se recorren los frames de una particula para promediarlos
	//se puede hacer con o sin shifts utilizando applyAlign, shiftX y shiftY

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

	//vectorAvg contendra el numero de pixeles asignados a cierto rango de amplitudes de la proyeccion DIRECT_A2D_ELEM(vectorAvg, 0, pos)+=1;
	//y la suma total de electrones asignados a ese rango DIRECT_A2D_ELEM(vectorAvg, 1, pos)+=DIRECT_MULTIDIM_ELEM(Iavg, n);
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


void ProgParticlePolishingMpi::calculateCurve_2(const MultidimArray<double> &Iproj, MultidimArray<double> &vectorAvg, int nStep, double &slope, double &intercept, double Dmin, double Dmax){

	std::vector<double> x, y;
	double aux;
	for(int n=0; n<nStep; n++){
		if (DIRECT_A2D_ELEM(vectorAvg, 0, n)!=0){
			DIRECT_A2D_ELEM(vectorAvg, 1, n)/=DIRECT_A2D_ELEM(vectorAvg, 0, n);
			aux = Dmin + n*((Dmax-Dmin)/(nStep-1));
			//En x e y tenemos valor de proyeccion, numero de electrones promedio
			x.push_back(double(aux));
			y.push_back(DIRECT_A2D_ELEM(vectorAvg, 1, n));
		}
	}

	//Usamos x e y para calcular pendiente y punto de intercepcion
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


/*
void ProgParticlePolishingMpi::calculateWeightedFrequency2(MultidimArray<double> &Ipart, const double aa, const double bb, const double cc){

	MultidimArray< std::complex<double> > fftIpart;
	MultidimArray<double> wfftIpart;
	FourierTransformer transformer;
	// Auxiliary vector for representing frequency values
	Matrix1D<double> fidx;

	transformer.FourierTransform(Ipart, fftIpart, false);

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
    transformer.inverseFourierTransform();

}
*/




void ProgParticlePolishingMpi::calculateWeightedFrequency3(MultidimArray<double> &Ipart, const std::vector<double> aa,
		const std::vector<double> bb, const std::vector<double> cc, const int nFrames, const int posi){

	MultidimArray< std::complex<double> > fftIpart;
	FourierTransformer transformer;
	// Auxiliary vector for representing frequency values
	Matrix1D<double> fidx;

	//en espacio de fourier aplicamos los pesos, que deben estar normalizados
	//se normalizan dividiendo el peso para un frame y frecuencia, entre la suma de los pesos para todos los frames en esa frecuencia

	transformer.FourierTransform(Ipart, fftIpart, false);

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
                	myw+=aa[nn]*exp(-bb[nn]*absw)+cc[nn]; //suma de pesos para todos los frames para normalizar
                }
                DIRECT_A3D_ELEM(fftIpart,k,i,j)=DIRECT_A3D_ELEM(fftIpart,k,i,j)*(aa[posi]*exp(-bb[posi]*absw)+cc[posi])/myw;
            }
        }
    }
    transformer.inverseFourierTransform();

}



void ProgParticlePolishingMpi::startProcessing()
{

	//TO ESTIMATE THE CURVE FOR EVERY MIC
	MetaData mdPartPrev, mdPart;
	MDRow currentRow;
	mdPartPrev.read(fn_in,NULL);
	mdPartSize = mdPartPrev.size();
	mdPart.sort(mdPartPrev, MDL_PARTICLE_ID);

	MDIterator *iterPart = new MDIterator();

	FileName fnPart;
    Image<double> I, projV, Iavg;
	Matrix2D<double> A;
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
	//size_t newId=0;


	//PROJECT THE INPUT VOLUME
	V.read(fnVol);
    V().setXmippOrigin();
	int xdimVol = (int)XSIZE(V());
	int ydimVol = (int)YSIZE(V());
	projectorV = new FourierProjector(V(),2,0.5,BSPLINE3);


	//Vamos recorriendo las micrografias, para cada una de ellas calcularemos la curva amplitud vs electrones
	for(int m=1; m<=nMics; m++){

		dataInMovie = 0;
		vectorAvg.initZeros(2, nStep);
		slope=0.;
		intercept=0.;

		iterPart->init(mdPart);
		//recorremos todo el metadata con las particulas de entrada
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
			//nos quedamos con las particulas que corresponden a la micrografia que estamos estudiando
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
			//MDL_DM3_VALUE lo he utilizado en un nuevo protocolo para extraer particulas a nivel de frame
			//este vector de doubles contiene, para cada particula,
			//los indices en el stack file que contienen cada uno de sus frames
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

			//proyectamos el volumen en la direccion indicada por la particula
			projectVolume(*projectorV, PV, xdimVol, xdimVol,  rot, tilt, psi);
			applyGeometry(LINEAR,projV(),PV(),A,IS_INV,DONT_WRAP,0.);
			projV().setXmippOrigin();

			//filtering the projections with the ctf
			ctf.applyCTF(projV(), samplingRate, false);

			//guardamos la proyeccion filtrada con la ctf porque me sera util mas adelante
			//asi evitamos proyectar el volumen varias veces por particula, que es muy lento
			FileName fnAux;
			fnAux = fnPart.insertBeforeExtension("proj");
			Image<double> aux=projV();
			aux.write(fnAux);

			//to obtain the points of the curve (intensity in the projection) vs (counted electrons)
			//the movie particles are averaged (all frames) to compare every pixel value
			bool noCurrent=false;
			bool applyAlign=false;
			std::vector<double> vectorX, vectorY;
			averagingAll(Iavg(), stks, fnPart, noCurrent, applyAlign, vectorX, vectorY, Xdim, Ydim);

			/*/DEBUG
			projV.write(formatString("Testprojection_%i_%i.tif", frId, partId));
			I.write(formatString("Testparticle_%i_%i.tif", frId, partId));
			Iavg.write(formatString("Testaverage_%i_%i.tif", frId, partId));
			//END DEBUG*/

			//With Iavg and projV, we calculate the curve (intensity in the projection) vs (counted electrons)
			//con la primera particula establezco los valores maximos y minimos que vamos a permitir
			//no se me ocurrio otra forma de hacerlo
			if(Dmin==0. && Dmax==0.){
				projV().computeDoubleMinMax(Dmin, Dmax);
				Dmin=Dmin-0.2*Dmin;
				Dmax=Dmax+0.2*Dmax;
				stepCurve = (Dmax-Dmin)/double(nStep);
				offset = -Dmin/double(stepCurve);
			}
			//basicamente esta funcion va construyendo un histograma
			//eje x-> valores de amplitud en la proyeccion (desde Dmin hasta Dmax hay nStep posibles niveles)
			//eje y-> numero de electrones promedio, que es lo mismo que los valores de amplitud en la particula promediada usando los frames
			calculateCurve_1(Iavg(), projV(), vectorAvg, nStep, stepCurve, offset, Dmin, Dmax);

			if(iterPart->hasNext())
				iterPart->moveNext();


		}//end particles loop


		//Una vez recorridas las particulas de esa micrografia calculamos los parametros de la curva
		//en realidad es una recta, solo hay pendiente y punto de intercepcion
		if(dataInMovie>0){
			calculateCurve_2(projV(), vectorAvg, nStep, slope, intercept, Dmin, Dmax);
			//Vectors to store some results for every mic
			slopes.push_back(slope);
			intercepts.push_back(intercept);
			mvIdsAux.push_back(m);
			std::cerr << "Estimated curve for movie " << m << ". Slope: " << slope << ". Intercept: " << intercept << std::endl;
		}

	}//end loop mics

	//Guardo los datos en un archivo que me servira para leer los datos facilmente luego al procesar cada particula
	//este archivo podria utilizarse tambien para, comproabr si existe, y de esa forma en un hipotetico re-run, no repetir esta parte (bastante lenta)
	std::ofstream outdata;
	outdata.open("moviedata.txt"); // opens the file
    if( !outdata ) { // file couldn't be opened
    	std::cerr << "Error: file could not be opened" << std::endl;
	    exit(1);
    }
    for (int a=0; a<mvIdsAux.size(); a++)
    	outdata << mvIdsAux[a] << " " << slopes[a] << " " << intercepts[a] << std::endl;
    outdata.close();

}



/*
//Just to remove intermediate files
void ProgParticlePolishingMpi::finishProcessing()
{

    std::cout << "Finish processing" << std::endl;
    XmippMetadataProgram::finishProcessing();

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


//funcion para ajustar mediante powell a una exponencial negativa de la forma a*exp(-b*x)+c
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
	Image<double> Ipart, projV, Iavg, Ifinal, Ipartaux, projVaux;
	std::vector<double> vectorX, vectorY;

	////////////////////////////////////
	/// FIRST PART - SHIFT ESTIMATION
	////////////////////////////////////

	//ML estimation of shifts
	//Possible shift values to check
	//-11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
	//-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6
	int shiftX[]={-3, -2, -1, 0, 1, 2, 3};
	int shiftY[]={-3, -2, -1, 0, 1, 2, 3};
    int myL = (int)(sizeof(shiftX)/sizeof(*shiftX));
	MultidimArray<double> lkresults;
	double maxShift, maxShift2;
	double finalPosX[nFrames];
	double finalPosY[nFrames];
	size_t frId, mvId, partId;
	int enabled;
	double slope=0., intercept=0.;
	std::vector<double> stks;

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

		//Necesitamos los indices del stack que contienen cada frame de esta particula para luego recorrerlos y calcular su shift
		rowIn.getValue(MDL_DM3_VALUE,stks);

		//reading projection
		//esta proyeccion ya esta filtrada con la ctf
		Image<double> projVread;
		FileName fnAux2=fnAux.insertBeforeExtension("proj");
		projVaux.read(fnAux2);
		projVaux().setXmippOrigin();

		//leemos la pendiente y el punto de intercepcion que corresponden con esta mic
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

		//Recorremos los frames de esta particula
		for (int a=0; a<nFrames; a++){

			//reading movie particle
			int id = (int)stks[a];
			fnPart.compose(id, fnImg.removePrefixNumber());
			Ipart.read(fnPart);
			Ipart().setXmippOrigin();
			maxShift = XSIZE(Ipart())/2-10;
			maxShift2 = maxShift*maxShift;

			//bucles para aplicar los shifts en las direcciones x e y
			lkresults.initZeros(myL, myL);
			for(int jj=0; jj<myL; jj++){
				for(int hh=0; hh<myL; hh++){

					projV().initZeros(projVaux());
					projV().setXmippOrigin();

					//aplicamos el shift, si ambos son cero no hace falta usar translate
					//projV contiene la proyeccion desplazada con los shifts correspondientes
					if (shiftX[jj]!=0 || shiftY[hh]!=0){
						translate(LINEAR, projV(), projVaux(), vectorR2((double)shiftX[jj], (double)shiftY[hh]), DONT_WRAP, 0.0); //translation of projV to avoid interpolation problem moving the frame particle
					}else{
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(projVaux()){
							DIRECT_MULTIDIM_ELEM(projV(),n) = DIRECT_MULTIDIM_ELEM(projVaux(),n);
						}
					}
					//aplicamos ML, es decir, lambda nos dice el valor esperado dada cierta amplitud en un pixel de la proyeccion (aplicando la recta antes calculada)
					//Vemos como de probable es obtener el numero de pixeles que nos dice Ipart en cada posicion dado ese lambda
					double likelihood=0.;
					double lambda, fact;
					for(int n=0; n<YSIZE(Ipart()); n++){
						for(int m=0; m<XSIZE(Ipart()); m++){
							//nos fijamos solo en la parte central de la particula
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
					//la verosimilitud se va almacenando en la matriz lkresults
					DIRECT_A2D_ELEM(lkresults, hh, jj) = likelihood;
				}
			}


			//utilizamos la funcion bestshift que nos dice donde se encuentra el maximo de lkresults
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

			//almacenamos resultados en finalPosX y finalPosY
			//negative because these values are calculated with a displacement in the projection,
			//but then they will be applied to the frames
			finalPosX[a]=-shiftX[bestPosX];
			finalPosY[a]=-shiftY[bestPosY];

			//printf(". BEST POS for particle %d. Shift %lf, %lf \n", a, finalPosX[a], finalPosY[a]);

		} // end loop frames

		//suavizamos las estimaciones anteriores porque son muy ruidosas ajustando a un polinomio de grado 2
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
			y(jj)=finalPosX[jj];
		}
		Matrix1D<double> alpha(3);
		solveLinearSystem(pseudoInverter, alpha);
		matrixOperation_Ax(X, alpha, resultX);
		//para la estimacion del shift en x

		y.resizeNoCopy(nFrames);
		for(int jj=0; jj<nFrames; jj++){
			y(jj)=finalPosY[jj];
		}
		solveLinearSystem(pseudoInverter, alpha);
		matrixOperation_Ax(X, alpha, resultY);
		//para la estimacion del shift en y

		//guardamos resultados, las variables de metadata MDL_POLISHING_X y MDL_POLISHING_Y las utilizo para luego poder analizar resultados
		for(int a=0; a<nFrames; a++){
			vectorX.push_back(round(resultX(a)));
			vectorY.push_back(round(resultY(a)));
		}
		rowOut.setValue(MDL_POLISHING_X, vectorX);
		rowOut.setValue(MDL_POLISHING_Y, vectorY);

		//Writing intermediate output - remove in some moment
		Ifinal().initZeros(projVaux());
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
		Ifinal.write(myFnOut);


		////////////////////////////////////
		/// SECOND PART - WEIGHT ESTIMATION
		////////////////////////////////////

		MultidimArray<double> matrixWeightsPart;
		FourierFilter FilterLP;
		FilterLP.FilterBand=LOWPASS;
		FilterLP.FilterShape=RAISED_COSINE;
		double cutfreq;
		nFilters=15;
		//frecuencias que vamos a usar
		double frC[15]={0.0078, 0.01, 0.0156, 0.0313, 0.05, 0.0625, 0.10, 0.1250, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45};

		matrixWeightsPart.initZeros(nFrames, nFilters);

		//To invert contrast in the projections
		double Dmin, Dmax, irange, val;
		projVaux().computeDoubleMinMax(Dmin, Dmax);
		irange=1.0/(Dmax - Dmin);

		//recorremos los frames de esta particula
		for(int a; a<nFrames; a++){

			//Reading the frame
			int id = (int)stks[a];
			fnPart.compose(id, fnImg.removePrefixNumber());

			for(int n=0; n<nFilters; n++){

				projV().initZeros(projVaux());
				projV().setXmippOrigin();
				//To invert contrast in the projections
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(projVaux()){
					val=DIRECT_MULTIDIM_ELEM(projVaux(),n);
					DIRECT_MULTIDIM_ELEM(projV(),n) = (Dmax - val) * irange;
				}

				//averaging movie particle image with the frames inside a window of size 4, and applying the align to all of them
				bool applyAlign=true;
				averagingWindow(Iavg(), stks, fnPart, applyAlign, vectorX, vectorY, Xdim, Ydim, 4);

				//filtering the projected particles with the lowpass filter
				FilterLP.w1=frC[n];
				FilterLP.generateMask(projV());
				FilterLP.applyMaskSpace(projV());

				//filtering the averaged movie particle
				FilterLP.generateMask(Iavg());
				FilterLP.applyMaskSpace(Iavg());

				//DEBUG
				//Iavg.write(formatString("myAvgFiltered.mrc"));
				//projV.write(formatString("myProjFiltered.mrc"));
				//END DEBUG

				//calculate similarity measures between averaged movie particles and filtered projection
				//basicamente calculamos la correlacion entre ambas y se alamcena ese valor como peso
				double corrN, corrM, corrW, imed;
				similarity(projV(), Iavg(), corrN, corrM, corrW, imed);
				DIRECT_A2D_ELEM(matrixWeightsPart, a, n) = corrN;

			} //end frequencies loop

		} //end loop frames


		//Guardamos en el metadata MDL_POLISHING_FREQ_COEFFS_BEFORE los valores de correlacion obtenidos para poder analizarlos facilmente
		std::vector<double> vectorAux;
		FOR_ALL_ELEMENTS_IN_ARRAY2D(matrixWeightsPart)
			vectorAux.push_back(A2D_ELEM(matrixWeightsPart,i,j));
		rowOut.setValue(MDL_POLISHING_FREQ_COEFFS_BEFORE, vectorAux);

		//Adjust to a negative exponential function
		//utilizando el optmizador powell
		MultidimArray<double> data;
		std::vector <double> aa,bb,cc;

		for(int mm=0; mm<nFrames; mm++){
			data.initZeros(nFilters);
			for(int nn=0; nn<nFilters; nn++){
				DIRECT_A1D_ELEM(data, nn) = DIRECT_A2D_ELEM(matrixWeightsPart, mm, nn);
			}
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
			//printf("Obtained params for frame %d: %lf, %lf, %lf with cost %lf and iters %d \n", mm, p(0), p(1), p(2), cost, iter);
		}

		std::vector<double> vectorAux3;
		for(int a=0; a<nFrames; a++){
			vectorAux3.push_back(aa[a]);
			vectorAux3.push_back(bb[a]);
			vectorAux3.push_back(cc[a]);
		}
		//guardamos los pesos de la funcion exponencial obtenidos para poder analizarlos en la variable MDL_POLISHING_FREQ_COEFFS_AFTER_SVD
		rowOut.setValue(MDL_POLISHING_FREQ_COEFFS_AFTER_SVD, vectorAux3);

		//Writing the ouput
		//para guardar la salida debemos aplicar el shift previamente calculado
		//luego la funcion exponencial nos da el peso a aplicar en cada frecuencia, que se hace en la funcion calculateWeightedFrequency3
		//se suman todos los frames y tenemos la particula final
		Ifinal().initZeros(projV());
		Ifinal().setXmippOrigin();
		for(int a=0; a<nFrames; a++){

			int id = (int)stks[a];
			fnPart.compose(id, fnImg.removePrefixNumber());
			Ipartaux.read(fnPart);
			Ipartaux().setXmippOrigin();
			selfTranslate(NEAREST, Ipartaux(), vectorR2(vectorX[a], vectorY[a]), DONT_WRAP, 0.0);
			calculateWeightedFrequency3(Ipartaux(), aa, bb, cc, nFrames, a);
			Ifinal()+=Ipartaux();

		}

		rowOut.setValue(MDL_IMAGE, fnImgOut);
		Ifinal.write(fnImgOut);
		printf("Particle %s finished \n", fnImgOut.getString().c_str());

	} //end if enabled


}



