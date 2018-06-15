#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:     I. Foche Perez (ifoche@cnb.csic.es)
# *              J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import join
from glob import glob
from datetime import datetime

Import('env')

AddOption('--no-opencv', dest='opencv', action='store_false', default=True,
          help='Avoid compilation of opencv programs')


# Define some variables used by Scons. Note that some of
# the variables will be passed by Scipion in the environment (env).

env['CUDA_SDK_PATH'] = os.environ.get('CUDA_SDK_PATH', '')
env['CUDA_LIB_PATH'] = os.environ.get('CUDA_LIB_PATH', '')

get = lambda x: os.environ.get(x, '0').lower() in ['true', 'yes', 'y', '1']

gtest = get('GTEST')
cuda = get('CUDA')
debug = get('DEBUG')
matlab = get('MATLAB')
opencv = env.GetOption('opencv') and get('OPENCV')

if opencv:
    opencvLibs = ['opencv_core',
                  'opencv_imgproc',
                  'opencv_video',
                  'libopencv_calib3d']
                  
else:
    opencvLibs = []



if 'MATLAB' in os.environ:
    # Must be removed to avoid problems in Matlab compilation.
    del os.environ['MATLAB']
    # Yeah, nice

# Read some flags
CYGWIN = env['PLATFORM'] == 'cygwin'
MACOSX = env['PLATFORM'] == 'darwin'
MINGW = env['PLATFORM'] == 'win32'

XMIPP_PATH = Dir('.').abspath
XMIPP_BUNDLE = Dir('..').abspath


#  ***********************************************************************
#  *                      Xmipp C++ Libraries                            *
#  ***********************************************************************

ALL_LIBS = {'fftw3', 'tiff', 'jpeg', 'sqlite3', 'hdf5', 'hdf5_cpp','XmippCore'}

# Create a shortcut and customized function
# to add the Xmipp CPP libraries
def addLib(name, **kwargs):
    ALL_LIBS.add(name)
    # Install all libraries in scipion/software/lib
    # COSS kwargs['installDir'] = '#software/lib'
    # Add always the xmipp path as -I for include and also xmipp/libraries
    incs = kwargs.get('incs', []) + [join(XMIPP_PATH, 'external'),
                                     join(XMIPP_PATH, 'libraries'), join(XMIPP_BUNDLE,'xmippCore')]
    kwargs['incs'] = incs

    deps = kwargs.get('deps', [])
    kwargs['deps'] = deps

    # Add libraries in libs as deps if not present
    libs = kwargs.get('libs', [])
    for lib in libs:
        if lib in ALL_LIBS and lib not in deps:
            deps.append(lib)

    # If pattern not provided use *.cpp as default
    patterns = kwargs.get('patterns', '*.cpp')
    kwargs['patterns'] = patterns
    
    libpath = kwargs.get('libpath', [])
    kwargs['libpath'] = libpath+[join(XMIPP_BUNDLE,'xmippCore','lib'),join(XMIPP_BUNDLE,'xmipp','lib')]

    if 'cuda' in kwargs and kwargs['cuda']:
    	lib = env.AddCppLibraryCuda(name, **kwargs)
    else:    	
    	lib = env.AddCppLibrary(name, **kwargs)
    	
    env.Alias('xmipp-libs', lib)

    return lib


# Add first external libraries (alglib, bilib, condor)
#NOTE: for alglib and condor, the dir can not be where the source
# code is because try to use .h files to make the final .so library

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]
env['PYTHONINCFLAGS'] = os.environ.get('PYTHONINCFLAGS', '').split()
if len(env["PYTHONINCFLAGS"])>0:
    python_incdirs = [remove_prefix(os.path.expandvars(x),"-I") for x in env["PYTHONINCFLAGS"]]
else:
    python_incdirs = []

# Basic libraries
dirs = ['external','external','external','libraries','libraries','libraries','libraries','libraries']
patterns=['condor/*.cpp','delaunay/*.cpp','gtest/*.cc','data/*.cpp','reconstruction/*.cpp',
		'classification/*.cpp','dimred/*.cpp','interface/*.cpp']
if cuda:
       dirs+=['libraries']
       patterns+=['reconstruction_adapt_cuda/*.cpp']
addLib('Xmipp', dirs=dirs, patterns=patterns, incs=python_incdirs, libs=['pthread','python2.7'])


# CUDA
if cuda:
    addLib('XmippReconsCuda',
       dirs=['libraries'],
       patterns=['reconstruction_cuda/*.cpp'],
       cuda=True)

# MPI
dirs = ['libraries']
patterns=['parallel/*.cpp']
if cuda:
	dirs+=['libraries']
	patterns+=['parallel_adapt_cuda/*.cpp']
addLib('XmippParallel',dirs=dirs,patterns=patterns, libs=['pthread', 'fftw3_threads','Xmipp'], mpi=True)


#  ***********************************************************************
#  *                      Xmipp Programs and Tests                       *
#  ***********************************************************************

XMIPP_LIBS = ['Xmipp']
PROG_DEPS = XMIPP_LIBS

PROG_LIBS = XMIPP_LIBS + ['hdf5','hdf5_cpp']

def addRunTest(testName, prog):
    """ Add a Scons target for running xmipp tests. """
    xmippTestName = 'xmipp_' + testName
    xmlFileName = join(XMIPP_PATH, 'applications', 'tests', 'OUTPUT',
                       xmippTestName+".xml")
    if os.path.exists(xmlFileName):
        os.remove(xmlFileName)
    testCase = env.Command(
        xmlFileName,
        join(XMIPP_PATH, 'bin/%s' % xmippTestName),
        "%s/scipion run $SOURCE --gtest_output=xml:$TARGET" % os.environ['SCIPION_HOME'])
    env.Alias('run_' + testName, testCase)
    env.Depends(testCase, prog)
    env.Alias('xmipp-runtests', testCase)

    AlwaysBuild(testCase)

    return testCase


# Shortcut function to add the Xmipp programs.
def addProg(progName, **kwargs):
    """ Shortcut to add the Xmipp programs easily.
    Params:
        progName: the name of the program without xmipp_ prefix that will be added.
        if 'src' not in kwargs, add: 'applications/programs/progName' by default.
        if progName starts with 'mpi_' then mpi will be set to True.
    """
    isTest = progName.startswith('test_')

    progsFolder = 'tests' if isTest else 'programs'

    src = kwargs.get('src', [join('applications', progsFolder, progName)])

    kwargs['src'] = src
    # Add all xmipp libraries just in case
    kwargs['libs'] = kwargs.get('libs', []) + PROG_LIBS + ['XmippCore']
    kwargs['deps'] = PROG_DEPS

    # Add always the xmipp path as -I for include and also xmipp/libraries
    incs = kwargs.get('incs', []) + [join(XMIPP_PATH, 'external'),
                                     join(XMIPP_PATH, 'libraries'), join(XMIPP_BUNDLE,'xmippCore')]
    kwargs['incs'] = incs
    if 'libPaths' in kwargs:
        kwargs['libPaths'] +=[join(XMIPP_BUNDLE,'xmippCore','lib'),join(XMIPP_BUNDLE,'xmipp','lib')]
    else:
        kwargs['libPaths'] =[join(XMIPP_BUNDLE,'xmippCore','lib'),join(XMIPP_BUNDLE,'xmipp','lib')]

    if progName.startswith('mpi_'):
        kwargs['mpi'] = True
        kwargs['libs'] += ['XmippParallel']

    if progName.startswith('cuda_'):
        kwargs['cuda'] = True
        #kwargs['nvcc'] = True
        kwargs['libs'] += ['XmippReconsCuda']

    if progName.startswith('mpi_cuda_'):
        kwargs['cuda'] = True
        #kwargs['nvcc'] = True
        kwargs['libs'] += ['XmippReconsCuda','XmippParallelAdaptCuda']

    xmippProgName = 'xmipp_%s' % progName

    if progName.startswith('test_'):
        env.Alias('xmipp-tests', xmippProgName)
        addRunTest(progName, xmippProgName)

    prog = env.AddProgram(xmippProgName, **kwargs)

    # Add some aliases before return
    env.Alias(xmippProgName, prog)
    env.Alias('xmipp-programs', prog)

    return prog

# Define the list of all programs
PROGRAMS_WITH_PYTHON = ['volume_align_prog','mpi_classify_CLTomo_prog']
PROGRAMS_WITH_OPENCV = ['movie_optical_alignment_cpu','mpi_volume_homogenizer']
PROGRAMS_WITH_OPENCV_AND_CUDA = ['movie_optical_alignment_gpu']
for p in glob(os.path.join(XMIPP_PATH,'applications','programs','*')):
	pname = os.path.basename(p)
	if pname in PROGRAMS_WITH_PYTHON:
		addProg(pname,incs=python_incdirs,libs=['python2.7'])
	elif pname in PROGRAMS_WITH_OPENCV:
		addProg(pname,libs=opencvLibs,deps=['opencv'])
	elif pname in PROGRAMS_WITH_OPENCV_AND_CUDA:
		if cuda:
			addProg(pname,libs=opencvLibs,deps=['opencv'],cuda=True)
	elif  'cuda' in pname:
		if cuda:
			addProg(pname)
	else:
		addProg(pname)

addLib('xmipp.so',
       dirs=['bindings'],
       patterns=['python/*.cpp'],
       incs=python_incdirs,
       libs=['python2.7', 'XmippCore', 'Xmipp'],
       prefix='', target='xmipp')
       
#  ***********************************************************************
#  *                      Xmipp Scripts                                  *
#  ***********************************************************************
def addBatch(batchName, script, scriptFolder='applications/scripts'):
    """ Add a link to xmipp/bin folder prepending xmipp_ prefix.
    The script should be located in from xmipp root,
    by default in 'applications/scripts/'
    """
    xmippBatchName = 'xmipp_%s' % batchName
    batchLink = env.SymLink(join(XMIPP_PATH, 'bin', xmippBatchName),
                            join(XMIPP_PATH, scriptFolder, script))
    env.Alias('xmipp-batchs', batchLink)

    return batchLink


# Batches (apps)
for scriptName in glob(os.path.join(XMIPP_PATH,'applications','scripts','*.[ps]*')):
	dirName = os.path.basename(os.path.dirname(scriptName))
	addBatch(dirName,scriptName)

# # Python tests
# testPythonInterface = env.SymLink('bin/xmipp_test_pythoninterface', 'applications/tests/test_pythoninterface/batch_test_pythoninterface.py')
# Depends(testPythonInterface, packageDeps)
# AddXmippTest('test_pythoninterface', testPythonInterface, "$SOURCE $TARGET")
#
# testPySqlite = env.SymLink('bin/xmipp_test_pysqlite', 'applications/tests/test_pysqlite/batch_test_pysqlite.py')
# Depends(testPySqlite, packageDeps)
# AddXmippTest('test_pysqlite', testPySqlite, "$SOURCE $TARGET")
#
# testEMX = env.SymLink('bin/xmipp_test_emx', 'applications/tests/test_emx/batch_test_emx.py')
# Depends(testEMX, packageDeps)
# AddXmippTest('test_emx', testEMX, "$SOURCE $TARGET")


def compileMatlabBinding(target, source, env):
    matlabDir = join(XMIPP_PATH, 'bindings', 'matlab')

    incStr = ' '.join('-I%s' % p for p in [os.path.join(XMIPP_BUNDLE, 'xmippCore'),
                                            os.path.join(XMIPP_BUNDLE, 'xmipp','libraries')]+env['EXTERNAL_INCDIRS'] )
    libStr = ' '.join('-L%s' % p for p in [os.path.join(XMIPP_BUNDLE, 'xmippCore','lib'),
                                            os.path.join(XMIPP_BUNDLE, 'xmipp','lib')]+env['EXTERNAL_LIBDIRS'])

    libs = ' '.join('-l%s' % lib for lib in ['Xmipp','XmippCore'])

    mex = join(env['MATLAB_DIR'], 'bin', 'mex')
    command = '%s -O CFLAGS="\$CFLAGS -std=c99 -fpermissive" -outdir %s %s %s %s %s ' % (mex, matlabDir, incStr, libStr, libs, source[0])
    print command
    os.system(command)

# Matlab programs
def addMatlabBinding(name):
    """ Add options to compile xmipp-Matlab bindings. """
    matlabDir = join(XMIPP_PATH, 'bindings', 'matlab')
    source = join(matlabDir, name)
    target = join(matlabDir, name.replace(".cpp",".mexa64"))

    cmdTarget = env.Command(target, source, compileMatlabBinding)
    env.Alias('xmipp-matlab', cmdTarget)

    return cmdTarget

if matlab:
	for pname in glob(os.path.join(XMIPP_PATH,'bindings','matlab','*.cpp')):
		addMatlabBinding(pname)
	env.Default('xmipp-matlab')
	env.Alias('xmipp', 'xmipp-matlab')

XmippAlias = env.Alias('xmipp', ['xmipp-libs',
                                 'xmipp-programs',
                                 'xmipp-batchs'])

Return('XmippAlias')
