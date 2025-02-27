/***************************************************************************
 *
 * Authors:  David Herreros Calero (dherreros@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <mpi.h>
#include <parallel/xmipp_mpi.h>
#include <reconstruction/forward_zernike_subtomos.h>


class MpiProgForwardZernikeSubtomos: public ProgForwardZernikeSubtomos, public MpiMetadataProgram
{

//AJ new
private:
	MpiFileMutex *fileMutex;
//END AJ

public:
    void defineParams()
    {
        ProgForwardZernikeSubtomos::defineParams();
        MpiMetadataProgram::defineParams();
    }
    void readParams()
    {
        MpiMetadataProgram::readParams();
        ProgForwardZernikeSubtomos::readParams();
    }
    void read(int argc, char **argv, bool reportErrors = true)
    {
    	//AJ new
    	fileMutex = new MpiFileMutex(node);
    	//END AJ
        MpiMetadataProgram::read(argc,argv);
    }
    void showProgress()
    {
        if (node->isMaster())
        {
            time_bar_done=first+1;
            ProgForwardZernikeSubtomos::showProgress();
        }
    }
    void preProcess()
    {
        //Master node should prepare some stuff before start working
        MetaData &mdIn = *getInputMd(); //get a reference to input metadata

        if (node->isMaster())
        {
        	ProgForwardZernikeSubtomos::createWorkFiles();
            // mdIn.write(fnOutDir + "/sphTodo.xmd");
        }
        node->barrierWait();//Sync all before start working
        ProgForwardZernikeSubtomos::preProcess();
        // mdIn.read(fnOutDir + "/sphTodo.xmd");
        mdIn.findObjects(imgsId);//get objects ids
        distributor = new MpiTaskDistributor(mdIn.size(), 1, node);
    }
    void startProcessing()
    {
        if (node->rank==1)
        {
        	verbose=1;
        	ProgForwardZernikeSubtomos::startProcessing();
        }
        node->barrierWait();
    }
    virtual bool getImageToProcess(size_t &objId, size_t &objIndex) override
    {
        //return getTaskToProcess(objId, objIndex);
        size_t first, last;
        bool moreTasks = distributor->getTasks(first, last);

        if (moreTasks)
        {
            time_bar_done = first + 1;
            objIndex = first;
            objId = imgsId[first];
            return true;
        }
        time_bar_done = getInputMd()->size();
        objId = BAD_OBJID;
        objIndex = BAD_INDEX;
        return false;
    }
    void finishProcessing()
    {
    	distributor->wait();

        //All nodes wait for each other
        node->barrierWait();
        if (node->isMaster())
        	ProgForwardZernikeSubtomos::finishProcessing();
        node->barrierWait();
    }
    //AJ new
    void checkPoint()
    {
        fileMutex->lock();
        ProgForwardZernikeSubtomos::checkPoint();
        fileMutex->unlock();
    }
    //END AJ
    void wait()
    {
		distributor->wait();
    }
};