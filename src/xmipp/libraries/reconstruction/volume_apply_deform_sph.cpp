/***************************************************************************
 *
 * Authors:    David Herreros Calero             dherreros@cnb.csic.es
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/numerical_tools.h>
#include <data/basis.h>
#include "volume_apply_deform_sph.h"
#include <data/numerical_tools.h>

void ProgApplyVolDeformSph::defineParams()
{
	addUsageLine("Deform a PDB according to a list of SPH deformation coefficients");
	addParamsLine("-i <volume>             : Volume to deform");
	addParamsLine("--clnm <metadata_file>  : List of deformation coefficients");
	addParamsLine("-o <volume>             : Deformed volume");
	addExampleLine("xmipp_apply_deform_sph -i input.vol -o volume_deformed.vol --clnm coefficients.txt");
}

void ProgApplyVolDeformSph::readParams()
{
	fn_vol=getParam("-i");
	fn_sph=getParam("--clnm");
	fn_out=getParam("-o");
}

void ProgApplyVolDeformSph::show()
{
	if (verbose==0)
		return;
	std::cout
	<< "Volume:               " << fn_vol << std::endl
	<< "Coefficient list:     " << fn_sph << std::endl
	<< "Output:               " << fn_out << std::endl
	;
}

void ProgApplyVolDeformSph::run()
{
	Image<double> VI, VO;
	VI.read(fn_vol);
	VI().setXmippOrigin();
	VO().initZeros(VI());
	VO().setXmippOrigin();
	int nCoeff = numberCoefficients();
	clnm.resize(nCoeff);
	clnm.read(fn_sph);
	int l,n,m;
	size_t idxY0=VEC_XSIZE(clnm)/4;
	size_t idxZ0=2*idxY0;
	size_t idxR=3*idxY0;
	const MultidimArray<double> &mVI=VI();
	double voxelI;
	for (int k=STARTINGZ(mVI); k<=FINISHINGZ(mVI); k++)
	{
		for (int i=STARTINGY(mVI); i<=FINISHINGY(mVI); i++)
		{
			for (int j=STARTINGX(mVI); j<=FINISHINGX(mVI); j++)
			{
				double gx=0.0, gy=0.0, gz=0.0;
				for (size_t idx=0; idx<idxY0; idx++)
				{
					double Rmax=VEC_ELEM(clnm,idx+idxR);
					double Rmax2=Rmax*Rmax;
					double iRmax=1.0/Rmax;
					double k2=k*k;
					double kr=k*iRmax;
					double k2i2=k2+i*i;
					double ir=i*iRmax;
					double r2=k2i2+j*j;
					double jr=j*iRmax;
					double rr=std::sqrt(r2)*iRmax;
					double zsph=0.0;
					if (r2<Rmax2)
					{
						spherical_index2lnm(idx,l,n,m);
						zsph=ZernikeSphericalHarmonics(l,n,m,jr,ir,kr,rr);
					}
					if (rr>0 || l==0)
					{
						gx += VEC_ELEM(clnm,idx)        *(zsph);
						gy += VEC_ELEM(clnm,idx+idxY0)  *(zsph);
						gz += VEC_ELEM(clnm,idx+idxZ0)  *(zsph);
					}
				}
				voxelI=mVI.interpolatedElement3D(j+gx,i+gy,k+gz);
				VO(k,i,j)=voxelI;
			}
		}
	}
	VO.write(fn_out);
}

int ProgApplyVolDeformSph::numberCoefficients()
{
	int nCoeff = 0;
	std::ifstream coeff_file;
	coeff_file.open(fn_sph.getString());
	float i;
	while (coeff_file >> i)
	{
		nCoeff++;
	}
	return nCoeff;
}