/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
 *
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

#ifndef MPI_ANGULAR_ACCURACY_PCA_H__
#define MPI_ANGULAR_ACCURACY_PCA_H_

#include <reconstruction/angular_accuracy_pca.h>
#include "parallel/xmipp_mpi.h"

/**@defgroup MpiProgAngularAccuracyPCA Determine the angular determination accuracy of a set of particles and a 3D reconstruction (MPI)
   @ingroup ParallelLibrary */
//@{

class MpiProgAngularAccuracyPCA: public ProgAngularAccuracyPCA
{
public:
	MpiNode *node=nullptr;
public:
	// Destructor
	MpiProgAngularAccuracyPCA() {}
	~MpiProgAngularAccuracyPCA();

	MpiProgAngularAccuracyPCA(const MpiProgAngularAccuracyPCA &)=delete;
	MpiProgAngularAccuracyPCA(const MpiProgAngularAccuracyPCA &&)=delete;
	MpiProgAngularAccuracyPCA & operator= (const MpiProgAngularAccuracyPCA &)=delete;
	MpiProgAngularAccuracyPCA & operator= (const MpiProgAngularAccuracyPCA &&)=delete;

	// Redefine how to read the command line
	void read(int argc, char** argv);

	// Redefine how to synchronize
	void synchronize();

	// Redefine how to gather the alignment
    void gatherResults();
};
//@}
#endif /*ANGULAR_ACCURACY_PCA_H_ */
