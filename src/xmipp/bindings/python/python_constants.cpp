/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
#include "xmippmodule.h"
#include "core/transformations_defines.h"
#include "core/metadata_sql.h"
#include "core/metadata_writemode.h"
#include "core/xmipp_datatype.h"
#include "core/xmipp_color.h"
#include "core/xmipp_image_base.h"
#include "core/axis_view.h"
#include "core/xmipp_random_mode.h"

/***************************************************************/
/*                    MDLabels constants                       */
/**************************************************************/

void addIntConstant(PyObject * dict, const char * name, const long &value)
{
    PyObject * pyValue = PyLong_FromLong(value);
    PyDict_SetItemString(dict, name, pyValue);
    Py_DECREF(pyValue);
}

//Macro to add constants to xmipp module with the same name
#define ADD_CONST(label) addIntConstant(dict, #label, (long)label)
// and with different name
#define ADD_CONST2(name, label) addIntConstant(dict, name, (long)label)

void addLabels(PyObject * dict)
{

    //Add constants
    ADD_CONST(xmipp_transformation::NEAREST);
    ADD_CONST(xmipp_transformation::LINEAR);
    ADD_CONST(xmipp_transformation::BSPLINE2);
    ADD_CONST(xmipp_transformation::BSPLINE3);
    ADD_CONST(xmipp_transformation::BSPLINE4);


    ADD_CONST(AGGR_COUNT);
    ADD_CONST(AGGR_MAX);
    ADD_CONST(AGGR_SUM);
    ADD_CONST(AGGR_AVG);
    ADD_CONST(AGGR_MIN);

    ADD_CONST(UNION);
    ADD_CONST(UNION_DISTINCT);
    ADD_CONST(INTERSECTION);
    ADD_CONST(SUBSTRACTION);
    ADD_CONST(INNER_JOIN);
    ADD_CONST(LEFT_JOIN);
    ADD_CONST(NATURAL_JOIN);
    ADD_CONST(OUTER_JOIN);
    ADD_CONST(INNER);
    ADD_CONST(LEFT);
    ADD_CONST(OUTER);
    ADD_CONST(NATURAL);
    ADD_CONST(EQ);
    ADD_CONST(NE);
    ADD_CONST(GT);
    ADD_CONST(LT);
    ADD_CONST(GE);
    ADD_CONST(LE);
    ADD_CONST(MD_OVERWRITE);
    ADD_CONST(MD_APPEND);
    ADD_CONST(MDL_UNDEFINED);
    ADD_CONST(MDL_FIRST_LABEL);

    //Metadata Labels
    ADD_CONST(MDL_OBJID);
    ADD_CONST(MDL_ANGLE_PSI2);
    ADD_CONST(MDL_ANGLE_PSI);
    ADD_CONST(MDL_ANGLE_PSI_DIFF);
    ADD_CONST(MDL_ANGLE_ROT2);
    ADD_CONST(MDL_ANGLE_ROT);
    ADD_CONST(MDL_ANGLE_ROT_DIFF);
    ADD_CONST(MDL_ANGLE_TILT2);
    ADD_CONST(MDL_ANGLE_TILT);
    ADD_CONST(MDL_ANGLE_TILT_DIFF);
    ADD_CONST(MDL_ANGLE_DIFF0);
    ADD_CONST(MDL_ANGLE_DIFF);
    ADD_CONST(MDL_ANGLE_DIFF2);
    ADD_CONST(MDL_ANGLE_Y);
    ADD_CONST(MDL_ANGLE_Y2);
    ADD_CONST(MDL_ANGLE_TEMPERATURE);
    ADD_CONST(MDL_ASSIGNED_DIR_REF_CC);
    ADD_CONST(MDL_AVG);
    ADD_CONST(MDL_AVG_CHANGES_ORIENTATIONS);
    ADD_CONST(MDL_AVG_CHANGES_OFFSETS);
    ADD_CONST(MDL_AVG_CHANGES_CLASSES);

    ADD_CONST(MDL_BGMEAN);
    ADD_CONST(MDL_BFACTOR);
    ADD_CONST(MDL_BLOCK_NUMBER);

    ADD_CONST(MDL_CLASS_COUNT);
    ADD_CONST(MDL_CLASS_PERCENTAGE);
    ADD_CONST(MDL_CLASSIFICATION_DATA);
    ADD_CONST(MDL_CLASSIFICATION_DATA_SIZE);
    ADD_CONST(MDL_CLASSIFICATION_DPR_05);
    ADD_CONST(MDL_CLASSIFICATION_FRC_05);
    ADD_CONST(MDL_CLASSIFICATION_INTRACLASS_DISTANCE);

    ADD_CONST(MDL_COLOR);
    ADD_CONST(MDL_COMMENT);
    ADD_CONST(MDL_COST);
    ADD_CONST(MDL_COST_PERCENTILE);
    ADD_CONST(MDL_COUNT);
    ADD_CONST(MDL_COUNT2);
    ADD_CONST(MDL_CONTINUOUS_GRAY_A);
    ADD_CONST(MDL_CONTINUOUS_GRAY_B);
    ADD_CONST(MDL_CONTINUOUS_X);
    ADD_CONST(MDL_CONTINUOUS_Y);
    ADD_CONST(MDL_CONTINUOUS_SCALE_X);
    ADD_CONST(MDL_CONTINUOUS_SCALE_Y);
    ADD_CONST(MDL_CONTINUOUS_FLIP);
    ADD_CONST(MDL_CORR_DENOISED_PROJECTION);
    ADD_CONST(MDL_CORR_DENOISED_NOISY);

    ADD_CONST(MDL_CORRELATION_IDX);
    ADD_CONST(MDL_CORRELATION_MASK);
    ADD_CONST(MDL_CORRELATION_WEIGHT);

    ADD_CONST(MDL_CTF_INPUTPARAMS);
    ADD_CONST(MDL_CTF_MODEL);
    ADD_CONST(MDL_CTF_MODEL2);
    ADD_CONST(MDL_CTF_SAMPLING_RATE);
    ADD_CONST(MDL_CTF_VOLTAGE);
    ADD_CONST(MDL_CTF_DEFOCUSU);
    ADD_CONST(MDL_CTF_DEFOCUSV);
    ADD_CONST(MDL_CTF_DEFOCUSA);
    ADD_CONST(MDL_CTF_DEFOCUS_ANGLE);
    ADD_CONST(MDL_CTF_DEFOCUS_CHANGE);
    ADD_CONST(MDL_CTF_DEFOCUS_R2);
    ADD_CONST(MDL_CTF_DEFOCUS_RESIDUAL);
    ADD_CONST(MDL_CTF_DEFOCUS_COEFS);
    ADD_CONST(MDL_CTF_DOWNSAMPLE_PERFORMED);
    ADD_CONST(MDL_CTF_CS);
    ADD_CONST(MDL_CTF_CA);
    ADD_CONST(MDL_CTF_GROUP);
    ADD_CONST(MDL_CTF_ENERGY_LOSS);
    ADD_CONST(MDL_CTF_LENS_STABILITY);
    ADD_CONST(MDL_CTF_CONVERGENCE_CONE);
    ADD_CONST(MDL_CTF_LONGITUDINAL_DISPLACEMENT);
    ADD_CONST(MDL_CTF_TRANSVERSAL_DISPLACEMENT);
    ADD_CONST(MDL_CTF_Q0);
    ADD_CONST(MDL_CTF_K);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN_K);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN_SIGMAU);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN_SIGMAV);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN_CU);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN_CV);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN_ANGLE);
    ADD_CONST(MDL_CTF_BG_SQRT_K);
    ADD_CONST(MDL_CTF_BG_SQRT_U);
    ADD_CONST(MDL_CTF_BG_SQRT_V);
    ADD_CONST(MDL_CTF_BG_SQRT_ANGLE);
    ADD_CONST(MDL_CTF_BG_BASELINE);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN2_K);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN2_SIGMAU);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN2_SIGMAV);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN2_CU);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN2_CV);
    ADD_CONST(MDL_CTF_BG_GAUSSIAN2_ANGLE);
    ADD_CONST(MDL_CTF_CRIT_PSDCORRELATION90);
    ADD_CONST(MDL_CTF_CRIT_FIRSTZERORATIO);
    ADD_CONST(MDL_CTF_CRIT_FIRSTZEROAVG);
    ADD_CONST(MDL_CTF_CRIT_FIRSTZERODISAGREEMENT);
    ADD_CONST(MDL_CTF_CRIT_NORMALITY);
    ADD_CONST(MDL_CTF_CRIT_DAMPING);
    ADD_CONST(MDL_CTF_CRIT_PSDRADIALINTEGRAL);
    ADD_CONST(MDL_CTF_CRIT_FITTINGSCORE);
    ADD_CONST(MDL_CTF_CRIT_FITTINGCORR13);
    ADD_CONST(MDL_CTF_CRIT_ICENESS);
    ADD_CONST(MDL_CTF_CRIT_PSDVARIANCE);
    ADD_CONST(MDL_CTF_CRIT_PSDPCA1VARIANCE);
    ADD_CONST(MDL_CTF_CRIT_PSDPCARUNSTEST);
    ADD_CONST(MDL_CTF_CRIT_NONASTIGMATICVALIDITY);
    ADD_CONST(MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_DIFF_RATIO);
    ADD_CONST(MDL_CTF_CRIT_MAXFREQ);
    ADD_CONST(MDL_CTF_PHASE_SHIFT);
    ADD_CONST(MDL_CTF_VPP_RADIUS);
    ADD_CONST(MDL_COORD_CONSENSUS_SCORE);
    ADD_CONST(MDL_CRYSTAL_LATTICE_A);
    ADD_CONST(MDL_CRYSTAL_LATTICE_B);
    ADD_CONST(MDL_CRYSTAL_DISAPPEAR_THRE);
    ADD_CONST(MDL_CRYSTAL_SHFILE);
    ADD_CONST(MDL_CRYSTAL_ORTHO_PRJ);
    ADD_CONST(MDL_CRYSTAL_PROJ);
    ADD_CONST(MDL_CRYSTAL_CELLX);
    ADD_CONST(MDL_CRYSTAL_CELLY);
    ADD_CONST(MDL_CRYSTAL_SHIFTX);
    ADD_CONST(MDL_CRYSTAL_SHIFTY);
    ADD_CONST(MDL_CRYSTAL_SHIFTZ);
    ADD_CONST(MDL_CUMULATIVE_SSNR);

    ADD_CONST(MDL_DATATYPE);
    ADD_CONST(MDL_DATE);
    ADD_CONST(MDL_DEFGROUP);
    ADD_CONST(MDL_DIMRED);
    ADD_CONST(MDL_DIMENSIONS_3D);
    ADD_CONST(MDL_DIMENSIONS_2D);
    ADD_CONST(MDL_DM3_IDTAG);
    ADD_CONST(MDL_DM3_NODEID);
    ADD_CONST(MDL_DM3_NUMBER_TYPE);
    ADD_CONST(MDL_DM3_PARENTID);
    ADD_CONST(MDL_DM3_TAGCLASS);
    ADD_CONST(MDL_DM3_TAGNAME);
    ADD_CONST(MDL_DM3_SIZE);
    ADD_CONST(MDL_DM3_VALUE);

    ADD_CONST(MDL_ENABLED);

    ADD_CONST(MDL_FLIP);
    ADD_CONST(MDL_FOM);
    ADD_CONST(MDL_FRAME_ID);

    ADD_CONST(MDL_GATHER_ID);

    ADD_CONST(MDL_GRAPH_DISTANCE2MAX);
    ADD_CONST(MDL_GRAPH_DISTANCE2MAX_PREVIOUS);
    ADD_CONST(MDL_GRAPH_CC);
    ADD_CONST(MDL_GRAPH_CC_PREVIOUS);

    ADD_CONST(MDL_IDX);
    ADD_CONST(MDL_IMAGE);
    ADD_CONST(MDL_IMAGE_COVARIANCE);
    ADD_CONST(MDL_IMAGE_IDX);
    ADD_CONST(MDL_IMAGE_ORIGINAL);
    ADD_CONST(MDL_IMAGE_REF);
    ADD_CONST(MDL_IMAGE_RESIDUAL);
    ADD_CONST(MDL_IMAGE_TILTED);
    ADD_CONST(MDL_IMGMD);
    ADD_CONST(MDL_IMAGE1);
    ADD_CONST(MDL_IMAGE2);
    ADD_CONST(MDL_IMAGE3);
    ADD_CONST(MDL_IMAGE4);
    ADD_CONST(MDL_IMAGE5);

    ADD_CONST(MDL_IMED);

    ADD_CONST(MDL_INTSCALE);
    ADD_CONST(MDL_ITEM_ID);
    ADD_CONST(MDL_ITER);
    ADD_CONST(MDL_KSTEST);
    ADD_CONST(MDL_LL);
    ADD_CONST(MDL_MACRO_CMD);
    ADD_CONST(MDL_MACRO_CMD_ARGS);
    ADD_CONST(MDL_MAGNIFICATION);
    ADD_CONST(MDL_MASK);
    ADD_CONST(MDL_MAXCC);
    ADD_CONST(MDL_MAXCC_PERCENTILE);
    ADD_CONST(MDL_MAX);
    ADD_CONST(MDL_MAXCC_PREVIOUS);
    ADD_CONST(MDL_MICROGRAPH);
    ADD_CONST(MDL_MICROGRAPH_ID);
    ADD_CONST(MDL_MICROGRAPH_MOVIE);
    ADD_CONST(MDL_MICROGRAPH_MOVIE_ID);
    ADD_CONST(MDL_MICROGRAPH_PARTICLES);
    ADD_CONST(MDL_MICROGRAPH_ORIGINAL);
    ADD_CONST(MDL_MICROGRAPH_TILTED);
    ADD_CONST(MDL_MICROGRAPH_TILTED_ORIGINAL);
    ADD_CONST(MDL_MIN);
    ADD_CONST(MDL_MIRRORFRAC);
    ADD_CONST(MDL_MISSINGREGION_NR);
    ADD_CONST(MDL_MISSINGREGION_TYPE);
    ADD_CONST(MDL_MISSINGREGION_THY0);
    ADD_CONST(MDL_MISSINGREGION_THYF);
    ADD_CONST(MDL_MISSINGREGION_THX0);
    ADD_CONST(MDL_MISSINGREGION_THXF);
    ADD_CONST(MDL_MLF_CTF);
    ADD_CONST(MDL_MLF_WIENER);
    ADD_CONST(MDL_MLF_SIGNAL);
    ADD_CONST(MDL_MLF_NOISE);
    ADD_CONST(MDL_MODELFRAC);

    ADD_CONST(MDL_NEIGHBORS);
    ADD_CONST(MDL_NEIGHBOR);
    ADD_CONST(MDL_NEIGHBORHOOD_RADIUS);
    ADD_CONST(MDL_NMA);
    ADD_CONST(MDL_NMA_MODEFILE);
    ADD_CONST(MDL_NMA_COLLECTIVITY);
    ADD_CONST(MDL_NMA_MINRANGE);
    ADD_CONST(MDL_NMA_MAXRANGE);
    ADD_CONST(MDL_NMA_SCORE);
    ADD_CONST(MDL_NMA_ATOMSHIFT);
    ADD_CONST(MDL_NMA_EIGENVAL);
    ADD_CONST(MDL_NOISE_ANGLES);
    ADD_CONST(MDL_NOISE_PARTICLE_COORD);
    ADD_CONST(MDL_NOISE_COORD);
    ADD_CONST(MDL_NOISE_PIXEL_LEVEL);

    ADD_CONST(MDL_OPTICALFLOW_MEANX);
    ADD_CONST(MDL_OPTICALFLOW_MEANY);
    ADD_CONST(MDL_OPTICALFLOW_STDX);
    ADD_CONST(MDL_OPTICALFLOW_STDY);

    ADD_CONST(MDL_ORDER);
    ADD_CONST(MDL_ORIGIN_X);
    ADD_CONST(MDL_ORIGIN_Y);
    ADD_CONST(MDL_ORIGIN_Z);

    ADD_CONST(MDL_PARTICLE_ID);
    ADD_CONST(MDL_PICKING_AUTOPICKPERCENT);
    ADD_CONST(MDL_PICKING_PARTICLE_SIZE);
    ADD_CONST(MDL_PICKING_STATE);
    ADD_CONST(MDL_PICKING_MICROGRAPH_STATE);
    ADD_CONST(MDL_PICKING_TEMPLATES);
    ADD_CONST(MDL_PICKING_MANUALPARTICLES_SIZE);
    ADD_CONST(MDL_PICKING_AUTOPARTICLES_SIZE);
    ADD_CONST(MDL_PMAX);
    ADD_CONST(MDL_AVGPMAX);
    ADD_CONST(MDL_PROGRAM);
    ADD_CONST(MDL_PRJ_DIMENSIONS);
    ADD_CONST(MDL_PRJ_ANGFILE);
    ADD_CONST(MDL_PRJ_PSI_NOISE);
    ADD_CONST(MDL_PRJ_PSI_RANDSTR);
    ADD_CONST(MDL_PRJ_PSI_RANGE);
    ADD_CONST(MDL_PRJ_ROT_NOISE);
    ADD_CONST(MDL_PRJ_ROT_RANDSTR);
    ADD_CONST(MDL_PRJ_ROT_RANGE);
    ADD_CONST(MDL_PRJ_TILT_NOISE);
    ADD_CONST(MDL_PRJ_TILT_RANDSTR);
    ADD_CONST(MDL_PRJ_TILT_RANGE);
    ADD_CONST(MDL_PRJ_VOL);
    ADD_CONST(MDL_PSD);
    ADD_CONST(MDL_PSD_ENHANCED);

    ADD_CONST(MDL_RANDOMSEED);
    ADD_CONST(MDL_REF3D);
    ADD_CONST(MDL_REF);
    ADD_CONST(MDL_REF2);
    ADD_CONST(MDL_REFMD);
    ADD_CONST(MDL_RESIDUE);
    ADD_CONST(MDL_RESOLUTION_ANISOTROPY);
    ADD_CONST(MDL_RESOLUTION_DPR);
    ADD_CONST(MDL_RESOLUTION_ERRORL2);
    ADD_CONST(MDL_RESOLUTION_FRC);
    ADD_CONST(MDL_RESOLUTION_FRCRANDOMNOISE);
    ADD_CONST(MDL_RESOLUTION_FREQ);
    ADD_CONST(MDL_RESOLUTION_FREQ2);
    ADD_CONST(MDL_RESOLUTION_FREQREAL);
    ADD_CONST(MDL_RESOLUTION_FSO);
    ADD_CONST(MDL_RESOLUTION_LOG_STRUCTURE_FACTOR);
    ADD_CONST(MDL_RESOLUTION_SSNR);
    ADD_CONST(MDL_RESOLUTION_STRUCTURE_FACTOR);
    ADD_CONST(MDL_RESOLUTION_RFACTOR);
    ADD_CONST(MDL_RESOLUTION_LOCAL_RESIDUE);

    ADD_CONST(MDL_SAMPLINGRATE);
    ADD_CONST(MDL_SAMPLINGRATE_ORIGINAL);
    ADD_CONST(MDL_SAMPLINGRATE_X);
    ADD_CONST(MDL_SAMPLINGRATE_Y);
    ADD_CONST(MDL_SAMPLINGRATE_Z);
    ADD_CONST(MDL_SCALE);
    ADD_CONST(MDL_SELFILE);
    ADD_CONST(MDL_SERIE);
    ADD_CONST(MDL_SHIFT_X);
    ADD_CONST(MDL_SHIFT_Y);
    ADD_CONST(MDL_SHIFT_Z);
    ADD_CONST(MDL_SHIFT_X2);
    ADD_CONST(MDL_SHIFT_Y2);
    ADD_CONST(MDL_SHIFT_X_DIFF);
    ADD_CONST(MDL_SHIFT_Y_DIFF);
    ADD_CONST(MDL_SHIFT_DIFF0);
    ADD_CONST(MDL_SHIFT_DIFF);
    ADD_CONST(MDL_SHIFT_DIFF2);
    ADD_CONST(MDL_SIGMANOISE);
    ADD_CONST(MDL_SIGMAOFFSET);
    ADD_CONST(MDL_SIGNALCHANGE);
    ADD_CONST(MDL_STAR_COMMENT);
    ADD_CONST(MDL_STDDEV);
    ADD_CONST(MDL_SCORE_BY_ALIGNABILITY);  
    ADD_CONST(MDL_SCORE_BY_MIRROR);
    ADD_CONST(MDL_SCORE_BY_ALIGNABILITY_PRECISION_EXP);
    ADD_CONST(MDL_SCORE_BY_ALIGNABILITY_PRECISION_REF);
    ADD_CONST(MDL_SCORE_BY_ALIGNABILITY_PRECISION);
    ADD_CONST(MDL_SCORE_BY_ALIGNABILITY_ACCURACY);
    ADD_CONST(MDL_SCORE_BY_ALIGNABILITY_ACCURACY_EXP);
    ADD_CONST(MDL_SCORE_BY_ALIGNABILITY_ACCURACY_REF);
    ADD_CONST(MDL_SCORE_BY_ALIGNABILITY_NOISE);
    ADD_CONST(MDL_SCORE_BY_EMPTINESS);
    ADD_CONST(MDL_SCORE_BY_ENTROPY);
    ADD_CONST(MDL_SCORE_BY_GRANULO);
    ADD_CONST(MDL_SCORE_BY_GINI);
    ADD_CONST(MDL_SCORE_BY_HISTDIST);
    ADD_CONST(MDL_SCORE_BY_LBP);
    ADD_CONST(MDL_SCORE_BY_PCA_RESIDUAL_PROJ);
    ADD_CONST(MDL_SCORE_BY_PCA_RESIDUAL_EXP);
    ADD_CONST(MDL_SCORE_BY_PCA_RESIDUAL);
    ADD_CONST(MDL_SCORE_BY_SCREENING);
    ADD_CONST(MDL_SCORE_BY_VARIANCE);
    ADD_CONST(MDL_SCORE_BY_VAR);
    ADD_CONST(MDL_SCORE_BY_ZERNIKE);
    ADD_CONST(MDL_SCORE_BY_ZSCORE);       
    ADD_CONST(MDL_SPH_COEFFICIENTS);
    ADD_CONST(MDL_SPH_DEFORMATION);
    ADD_CONST(MDL_SPH_TSNE_COEFF1D);
    ADD_CONST(MDL_SPH_TSNE_COEFF2D);

    ADD_CONST(MDL_SUM);
    ADD_CONST(MDL_SUMWEIGHT);
    ADD_CONST(MDL_SYMNO);

    ADD_CONST(MDL_TIME);
    ADD_CONST(MDL_TRANSFORM_MATRIX);
    ADD_CONST(MDL_TOMOGRAM_VOLUME);
    ADD_CONST(MDL_TOMOGRAMMD);

    ADD_CONST(MDL_USER);

    ADD_CONST(MDL_VOLUME_SCORE_SUM);
    ADD_CONST(MDL_VOLUME_SCORE_SUM_TH);
    ADD_CONST(MDL_VOLUME_SCORE_MEAN);
    ADD_CONST(MDL_VOLUME_SCORE_MIN);
    ADD_CONST(MDL_VOLUME_SCORE1);
    ADD_CONST(MDL_VOLUME_SCORE2);
    ADD_CONST(MDL_VOLUME_SCORE3);
    ADD_CONST(MDL_VOLUME_SCORE4);

    ADD_CONST(MDL_WEIGHT);
    ADD_CONST(MDL_WEIGHT_P);
    ADD_CONST(MDL_WEIGHT_CONTINUOUS2);
    ADD_CONST(MDL_WEIGHT_JUMPER0);
    ADD_CONST(MDL_WEIGHT_JUMPER);
    ADD_CONST(MDL_WEIGHT_JUMPER2);
    ADD_CONST(MDL_WEIGHT_PRECISION_ALIGNABILITY);
    ADD_CONST(MDL_WEIGHT_ACCURACY_ALIGNABILITY);
    ADD_CONST(MDL_WEIGHT_ALIGNABILITY);
    ADD_CONST(MDL_WEIGHT_PRECISION_MIRROR);
    ADD_CONST(MDL_WEIGHT_SIGNIFICANT);
    ADD_CONST(MDL_WEIGHT_SSNR);
    ADD_CONST(MDL_WROBUST);

    ADD_CONST(MDL_XCOOR);
    ADD_CONST(MDL_XCOOR_TILT);
    ADD_CONST(MDL_XSIZE);
    ADD_CONST(MDL_X);

    ADD_CONST(MDL_YCOOR);
    ADD_CONST(MDL_YCOOR_TILT);
    ADD_CONST(MDL_Y);
    ADD_CONST(MDL_YSIZE);
    ADD_CONST(MDL_ZCOOR);
    ADD_CONST(MDL_Z);
    ADD_CONST(MDL_ZSCORE);
    ADD_CONST(MDL_ZSCORE_DEEPLEARNING1);
    ADD_CONST(MDL_GOOD_REGION_SCORE);
    ADD_CONST(MDL_ZSCORE_HISTOGRAM);
    ADD_CONST(MDL_ZSCORE_RESMEAN);
    ADD_CONST(MDL_ZSCORE_RESVAR);
    ADD_CONST(MDL_ZSCORE_RESCOV);
    ADD_CONST(MDL_ZSCORE_SHAPE1);
    ADD_CONST(MDL_ZSCORE_SHAPE2);
    ADD_CONST(MDL_ZSCORE_SNR1);
    ADD_CONST(MDL_ZSCORE_SNR2);
    ADD_CONST(MDL_ZSIZE);
    ADD_CONST(MDL_ZSIZE);
    ADD_CONST(MDL_LAST_LABEL);

    /* Label types */
    ADD_CONST(LABEL_NOTYPE);
    ADD_CONST(LABEL_INT);
    ADD_CONST(LABEL_BOOL);
    ADD_CONST(LABEL_DOUBLE);
    ADD_CONST(LABEL_VECTOR_DOUBLE);
    ADD_CONST(LABEL_STRING);
    ADD_CONST(LABEL_SIZET);
    ADD_CONST(LABEL_VECTOR_SIZET);
    ADD_CONST(TAGLABEL_NOTAG);
    ADD_CONST(TAGLABEL_TEXTFILE);
    ADD_CONST(TAGLABEL_METADATA);
    ADD_CONST(TAGLABEL_CTFPARAM);
    ADD_CONST(TAGLABEL_IMAGE);
    ADD_CONST(TAGLABEL_VOLUME);
    ADD_CONST(TAGLABEL_STACK);
    ADD_CONST(TAGLABEL_MICROGRAPH);
    ADD_CONST(TAGLABEL_PSD);
    ADD_CONST(_NONE);
    ADD_CONST(HEADER);
    ADD_CONST2("HEADER_ALL", _HEADER_ALL);
    ADD_CONST(DATA);
    ADD_CONST2("DATA_ALL", _DATA_ALL);
    ADD_CONST(xmipp_transformation::WRAP);
    ADD_CONST(ALL_IMAGES);
    ADD_CONST(FILENAMENUMBERLENGTH);
    ADD_CONST2("XMIPP_BLACK", BLACK);
    ADD_CONST2("XMIPP_RED", RED);
    ADD_CONST2("XMIPP_GREEN", GREEN);
    ADD_CONST2("XMIPP_YELLOW", YELLOW);
    ADD_CONST2("XMIPP_BLUE", BLUE);
    ADD_CONST2("XMIPP_MAGENTA", MAGENTA);
    ADD_CONST2("XMIPP_CYAN", CYAN);
    ADD_CONST2("XMIPP_WHITE", WHITE);
    ADD_CONST2("XMIPP_RND_UNIFORM", RND_UNIFORM);
    ADD_CONST2("XMIPP_RND_GAUSSIAN", RND_GAUSSIAN);

    ADD_CONST2("DT_DEFAULT", DT_Default);
    ADD_CONST2("DT_UNKNOWN", DT_Unknown);
    ADD_CONST2("DT_UCHAR", DT_UChar);
    ADD_CONST2("DT_SCHAR", DT_SChar);
    ADD_CONST2("DT_USHORT", DT_UShort);
    ADD_CONST2("DT_SHORT", DT_Short);
    ADD_CONST2("DT_UINT", DT_UInt);
    ADD_CONST2("DT_INT", DT_Int);
    ADD_CONST2("DT_LONG", DT_Long);
    ADD_CONST2("DT_FLOAT", DT_Float);
    ADD_CONST2("DT_DOUBLE", DT_Double);
    ADD_CONST2("DT_COMPLEXSHORT", DT_CShort);
    ADD_CONST2("DT_COMPLEXINT", DT_CInt);
    ADD_CONST2("DT_COMPLEXFLOAT", DT_CFloat);
    ADD_CONST2("DT_COMPLEXDOUBLE", DT_CDouble);
    ADD_CONST2("DT_BOOL", DT_Bool);
    ADD_CONST2("DT_LASTENTRY", DT_LastEntry);

    ADD_CONST(VIEW_Z_NEG);
    ADD_CONST(VIEW_Z_POS);
    ADD_CONST(VIEW_Y_NEG);
    ADD_CONST(VIEW_Y_POS);
    ADD_CONST(VIEW_X_NEG);
    ADD_CONST(VIEW_X_POS);

    ADD_CONST2("XMIPP_EQUAL_ACCURACY",XMIPP_EQUAL_ACCURACY);

    ADD_CONST(CW_CAST);
    ADD_CONST(CW_CONVERT);
    ADD_CONST(CW_ADJUST);
    ADD_CONST(CW_LAST_LABEL);

    /** RELION labels */
    ADD_CONST(RLN_AREA_ID); ///< ID for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
    ADD_CONST(RLN_AREA_NAME); ///< Name for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
    ADD_CONST(RLN_COMMENT); // The RLN_COMMENT is handled specially as well

    ADD_CONST(RLN_CTF_BFACTOR); ///< B-factor
    ADD_CONST(RLN_CTF_SCALEFACTOR); ///< linear scale-factor
    ADD_CONST(RLN_CTF_SAMPLING_RATE); ///< Sampling rate
    ADD_CONST(RLN_CTF_VOLTAGE); ///< Microscope voltage (kV)
    ADD_CONST(RLN_CTF_DEFOCUSU); ///< Defocus U (Angstroms)
    ADD_CONST(RLN_CTF_DEFOCUSV); ///< Defocus V (Angstroms)
    ADD_CONST(RLN_CTF_DEFOCUS_ANGLE); ///< Defocus angle (degrees)
    ADD_CONST(RLN_CTF_CS); ///< Spherical aberration
    ADD_CONST(RLN_CTF_CA); ///< Chromatic aberration
    ADD_CONST(RLN_CTF_DETECTOR_PIXEL_SIZE); ///< Pixel size for detector as used in CTF-determination
    ADD_CONST(RLN_CTF_ENERGY_LOSS); ///< Energy loss
    ADD_CONST(RLN_CTF_FOM); ///< ctffind3 FOM (CC) for quality of CTF-fit
    ADD_CONST(RLN_CTF_IMAGE); ///< name of an image describing the CTF model
    ADD_CONST(RLN_CTF_LENS_STABILITY); ///< Lens stability
    ADD_CONST(RLN_CTF_MAGNIFICATION); ///< Magnification used for CTF-determination
    ADD_CONST(RLN_CTF_CONVERGENCE_CONE); ///< Convergence cone
    ADD_CONST(RLN_CTF_LONGITUDINAL_DISPLACEMENT); ///< Longitudinal displacement
    ADD_CONST(RLN_CTF_TRANSVERSAL_DISPLACEMENT); ///< Transversal displacemente
    ADD_CONST(RLN_CTF_Q0); ///< Amplitude contrast
    ADD_CONST(RLN_CTF_K); ///< CTF gain
    ADD_CONST(RLN_CTF_VALUE); ///< CTF value
    ADD_CONST(RLN_CTF_PHASESHIFT); ///< Phase-shift from a phase-plate (in degrees)

    ADD_CONST(RLN_IMAGE_NAME);
    ADD_CONST(RLN_IMAGE_RECONSTRUCT_NAME);
    ADD_CONST(RLN_IMAGE_ID);
    ADD_CONST(RLN_IMAGE_ENABLED);
    ADD_CONST(RLN_IMAGE_DATATYPE);
    ADD_CONST(RLN_IMAGE_DIMENSIONALITY);
    ADD_CONST(RLN_IMAGE_BEAMTILT_X);
    ADD_CONST(RLN_IMAGE_BEAMTILT_Y);
    ADD_CONST(RLN_IMAGE_BEAMTILT_GROUP);
    ADD_CONST(RLN_IMAGE_COORD_X);
    ADD_CONST(RLN_IMAGE_COORD_Y);
    ADD_CONST(RLN_IMAGE_COORD_Z);
    ADD_CONST(RLN_IMAGE_FRAME_NR);
    ADD_CONST(RLN_IMAGE_MAGNIFICATION_CORRECTION);
    ADD_CONST(RLN_IMAGE_NORM_CORRECTION);
    ADD_CONST(RLN_IMAGE_SAMPLINGRATE);
    ADD_CONST(RLN_IMAGE_SAMPLINGRATE_X);
    ADD_CONST(RLN_IMAGE_SAMPLINGRATE_Y);
    ADD_CONST(RLN_IMAGE_SAMPLINGRATE_Z);
    ADD_CONST(RLN_IMAGE_SIZE);
    ADD_CONST(RLN_IMAGE_SIZEX);
    ADD_CONST(RLN_IMAGE_SIZEY);
    ADD_CONST(RLN_IMAGE_SIZEZ);
    ADD_CONST(RLN_IMAGE_STATS_MIN);
    ADD_CONST(RLN_IMAGE_STATS_MAX);
    ADD_CONST(RLN_IMAGE_STATS_AVG);
    ADD_CONST(RLN_IMAGE_STATS_STDDEV);
    ADD_CONST(RLN_IMAGE_STATS_SKEW);
    ADD_CONST(RLN_IMAGE_STATS_KURT);
    ADD_CONST(RLN_IMAGE_WEIGHT);

    ADD_CONST(RLN_MATRIX_1_1);
    ADD_CONST(RLN_MATRIX_1_2);
    ADD_CONST(RLN_MATRIX_1_3);
    ADD_CONST(RLN_MATRIX_2_1);
    ADD_CONST(RLN_MATRIX_2_2);
    ADD_CONST(RLN_MATRIX_2_3);
    ADD_CONST(RLN_MATRIX_3_1);
    ADD_CONST(RLN_MATRIX_3_2);
    ADD_CONST(RLN_MATRIX_3_3);

    ADD_CONST(RLN_MICROGRAPH_ID);
    ADD_CONST(RLN_MICROGRAPH_MOVIE_NAME);
    ADD_CONST(RLN_MICROGRAPH_NAME);
    ADD_CONST(RLN_MICROGRAPH_TILT_ANGLE);
    ADD_CONST(RLN_MICROGRAPH_TILT_AXIS_DIRECTION);
    ADD_CONST(RLN_MICROGRAPH_TILT_AXIS_OUTOFPLANE);

    ADD_CONST(RLN_MLMODEL_ACCURACY_ROT);
    ADD_CONST(RLN_MLMODEL_ACCURACY_TRANS);
    ADD_CONST(RLN_MLMODEL_AVE_PMAX);
    ADD_CONST(RLN_MLMODEL_CURRENT_RESOLUTION);
    ADD_CONST(RLN_MLMODEL_CURRENT_SIZE);
    ADD_CONST(RLN_MLMODEL_DATA_VS_PRIOR_REF);
    ADD_CONST(RLN_MLMODEL_DIMENSIONALITY);
    ADD_CONST(RLN_MLMODEL_DIMENSIONALITY_DATA);
    ADD_CONST(RLN_MLMODEL_DIFF2_HALVES_REF);
    ADD_CONST(RLN_MLMODEL_ESTIM_RESOL_REF);
    ADD_CONST(RLN_MLMODEL_FOURIER_COVERAGE_REF);
    ADD_CONST(RLN_MLMODEL_FOURIER_COVERAGE_TOTAL_REF);
    ADD_CONST(RLN_MLMODEL_FSC_HALVES_REF);
    ADD_CONST(RLN_MLMODEL_GROUP_NAME);
    ADD_CONST(RLN_MLMODEL_GROUP_NO);
    ADD_CONST(RLN_MLMODEL_GROUP_NR_PARTICLES);
    ADD_CONST(RLN_MLMODEL_GROUP_SCALE_CORRECTION);
    ADD_CONST(RLN_MLMODEL_INTERPOLATOR);
    ADD_CONST(RLN_MLMODEL_LL);
    ADD_CONST(RLN_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION);
    ADD_CONST(RLN_MLMODEL_NORM_CORRECTION_AVG);
    ADD_CONST(RLN_MLMODEL_NR_BODIES);
    ADD_CONST(RLN_MLMODEL_NR_CLASSES);
    ADD_CONST(RLN_MLMODEL_NR_GROUPS);
    ADD_CONST(RLN_MLMODEL_ORIGINAL_SIZE);
    ADD_CONST(RLN_MLMODEL_ORIENTABILITY_CONTRIBUTION);
    ADD_CONST(RLN_MLMODEL_PADDING_FACTOR);
    ADD_CONST(RLN_MLMODEL_PDF_CLASS);
    ADD_CONST(RLN_MLMODEL_PRIOR_OFFX_CLASS);
    ADD_CONST(RLN_MLMODEL_PRIOR_OFFY_CLASS);
    ADD_CONST(RLN_MLMODEL_PDF_ORIENT);
    ADD_CONST(RLN_MLMODEL_PIXEL_SIZE);
    ADD_CONST(RLN_MLMODEL_POWER_REF);
    ADD_CONST(RLN_MLMODEL_PRIOR_MODE);
    ADD_CONST(RLN_MLMODEL_SIGMA_OFFSET);
    ADD_CONST(RLN_MLMODEL_SIGMA_ROT);
    ADD_CONST(RLN_MLMODEL_SIGMA_TILT);
    ADD_CONST(RLN_MLMODEL_SIGMA_PSI);
    ADD_CONST(RLN_MLMODEL_REF_IMAGE);
    ADD_CONST(RLN_MLMODEL_SIGMA2_NOISE);
    ADD_CONST(RLN_MLMODEL_SIGMA2_REF);
    ADD_CONST(RLN_MLMODEL_SSNR_REF);
    ADD_CONST(RLN_MLMODEL_TAU2_FUDGE_FACTOR);
    ADD_CONST(RLN_MLMODEL_TAU2_REF);
    ADD_CONST(RLN_OPTIMISER_ACCURACY_ROT);
    ADD_CONST(RLN_OPTIMISER_ACCURACY_TRANS);
    ADD_CONST(RLN_OPTIMISER_ADAPTIVE_FRACTION);
    ADD_CONST(RLN_OPTIMISER_ADAPTIVE_OVERSAMPLING);
    ADD_CONST(RLN_OPTIMISER_AUTO_LOCAL_HP_ORDER);
    ADD_CONST(RLN_OPTIMISER_AVAILABLE_MEMORY);
    ADD_CONST(RLN_OPTIMISER_BEST_RESOL_THUS_FAR);
    ADD_CONST(RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS);
    ADD_CONST(RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS);
    ADD_CONST(RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES);
    ADD_CONST(RLN_OPTIMISER_COARSE_SIZE);
    ADD_CONST(RLN_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED);
    ADD_CONST(RLN_OPTIMISER_DATA_STARFILE);
    ADD_CONST(RLN_OPTIMISER_DO_AUTO_REFINE);
    ADD_CONST(RLN_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES);
    ADD_CONST(RLN_OPTIMISER_DO_CORRECT_CTF);
    ADD_CONST(RLN_OPTIMISER_DO_CORRECT_MAGNIFICATION);
    ADD_CONST(RLN_OPTIMISER_DO_CORRECT_NORM);
    ADD_CONST(RLN_OPTIMISER_DO_CORRECT_SCALE);
    ADD_CONST(RLN_OPTIMISER_DO_REALIGN_MOVIES);
    ADD_CONST(RLN_OPTIMISER_DO_MAP);
    ADD_CONST(RLN_OPTIMISER_DO_SOLVENT_FLATTEN);
    ADD_CONST(RLN_OPTIMISER_DO_SKIP_ALIGN);
    ADD_CONST(RLN_OPTIMISER_DO_SKIP_ROTATE);
    ADD_CONST(RLN_OPTIMISER_DO_SPLIT_RANDOM_HALVES);
    ADD_CONST(RLN_OPTIMISER_DO_ZERO_MASK);
    ADD_CONST(RLN_OPTIMISER_FIX_SIGMA_NOISE);
    ADD_CONST(RLN_OPTIMISER_FIX_SIGMA_OFFSET);
    ADD_CONST(RLN_OPTIMISER_FIX_TAU);
    ADD_CONST(RLN_OPTIMISER_HAS_CONVERGED);
    ADD_CONST(RLN_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT);
    ADD_CONST(RLN_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO);
    ADD_CONST(RLN_OPTIMISER_HIGHRES_LIMIT_EXP);
    ADD_CONST(RLN_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK);
    ADD_CONST(RLN_OPTIMISER_INCR_SIZE);
    ADD_CONST(RLN_OPTIMISER_ITERATION_NO);
    ADD_CONST(RLN_OPTIMISER_LOCAL_SYMMETRY_FILENAME);
    ADD_CONST(RLN_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES);
    ADD_CONST(RLN_OPTIMISER_MAGNIFICATION_RANGE);
    ADD_CONST(RLN_OPTIMISER_MAGNIFICATION_STEP);
    ADD_CONST(RLN_OPTIMISER_MAX_COARSE_SIZE);
    ADD_CONST(RLN_OPTIMISER_MAX_NR_POOL);
    ADD_CONST(RLN_OPTIMISER_MODEL_STARFILE);
    ADD_CONST(RLN_OPTIMISER_MODEL_STARFILE2);
    ADD_CONST(RLN_OPTIMISER_NR_ITERATIONS);
    ADD_CONST(RLN_OPTIMISER_NR_ITER_WO_RESOL_GAIN);
    ADD_CONST(RLN_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES);
    ADD_CONST(RLN_OPTIMISER_OUTPUT_ROOTNAME);
    ADD_CONST(RLN_OPTIMISER_PARTICLE_DIAMETER);
    ADD_CONST(RLN_OPTIMISER_RADIUS_MASK_3D_MAP);
    ADD_CONST(RLN_OPTIMISER_RADIUS_MASK_EXP_PARTICLES);
    ADD_CONST(RLN_OPTIMISER_RANDOM_SEED);
    ADD_CONST(RLN_OPTIMISER_REFS_ARE_CTF_CORRECTED);
    ADD_CONST(RLN_OPTIMISER_SAMPLING_STARFILE);
    ADD_CONST(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES);
    ADD_CONST(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS);
    ADD_CONST(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS);
    ADD_CONST(RLN_OPTIMISER_SOLVENT_MASK_NAME);
    ADD_CONST(RLN_OPTIMISER_SOLVENT_MASK2_NAME);
    ADD_CONST(RLN_OPTIMISER_TAU_SPECTRUM_NAME);
    ADD_CONST(RLN_OPTIMISER_USE_TOO_COARSE_SAMPLING);
    ADD_CONST(RLN_OPTIMISER_WIDTH_MASK_EDGE);
    ADD_CONST(RLN_ORIENT_FLIP);
    ADD_CONST(RLN_ORIENT_ID);
    ADD_CONST(RLN_ORIENT_ORIGIN_X);
    ADD_CONST(RLN_ORIENT_ORIGIN_X_PRIOR);
    ADD_CONST(RLN_ORIENT_ORIGIN_Y);
    ADD_CONST(RLN_ORIENT_ORIGIN_Y_PRIOR);
    ADD_CONST(RLN_ORIENT_ORIGIN_Z);
    ADD_CONST(RLN_ORIENT_ORIGIN_Z_PRIOR);
    ADD_CONST(RLN_ORIENT_ROT);
    ADD_CONST(RLN_ORIENT_ROT_PRIOR);
    ADD_CONST(RLN_ORIENT_TILT);
    ADD_CONST(RLN_ORIENT_TILT_PRIOR);
    ADD_CONST(RLN_ORIENT_PSI);
    ADD_CONST(RLN_ORIENT_PSI_PRIOR);
    ADD_CONST(RLN_ORIENT_PSI_PRIOR_FLIP_RATIO);

    ADD_CONST(RLN_PARTICLE_AUTOPICK_FOM);
    ADD_CONST(RLN_PARTICLE_CLASS);
    ADD_CONST(RLN_PARTICLE_DLL);
    ADD_CONST(RLN_PARTICLE_ID);
    ADD_CONST(RLN_PARTICLE_FOM);
    ADD_CONST(RLN_PARTICLE_KL_DIVERGENCE);
    ADD_CONST(RLN_PARTICLE_MOVIE_RUNNING_AVG);
    ADD_CONST(RLN_PARTICLE_RANDOM_SUBSET);
    ADD_CONST(RLN_PARTICLE_NAME);
    ADD_CONST(RLN_PARTICLE_NR_FRAMES);
    ADD_CONST(RLN_PARTICLE_NR_FRAMES_AVG);
    ADD_CONST(RLN_PARTICLE_ORI_NAME);
    ADD_CONST(RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES);
    ADD_CONST(RLN_PARTICLE_PMAX);

    // New helical labes in Relion 2.x
    ADD_CONST(RLN_PARTICLE_HELICAL_TUBE_ID);
    ADD_CONST(RLN_PARTICLE_HELICAL_TUBE_PITCH);
    ADD_CONST(RLN_PARTICLE_HELICAL_TRACK_LENGTH);
    ADD_CONST(RLN_MLMODEL_HELICAL_NR_ASU);
    ADD_CONST(RLN_MLMODEL_HELICAL_TWIST);
    ADD_CONST(RLN_MLMODEL_HELICAL_TWIST_MIN);
    ADD_CONST(RLN_MLMODEL_HELICAL_TWIST_MAX);
    ADD_CONST(RLN_MLMODEL_HELICAL_TWIST_INITIAL_STEP);
    ADD_CONST(RLN_MLMODEL_HELICAL_RISE);
    ADD_CONST(RLN_MLMODEL_HELICAL_RISE_MIN);
    ADD_CONST(RLN_MLMODEL_HELICAL_RISE_MAX);
    ADD_CONST(RLN_MLMODEL_HELICAL_RISE_INITIAL_STEP);
    ADD_CONST(RLN_OPTIMISER_DO_HELICAL_REFINE);
    ADD_CONST(RLN_OPTIMISER_HELICAL_TWIST_INITIAL);
    ADD_CONST(RLN_OPTIMISER_HELICAL_RISE_INITIAL);
    ADD_CONST(RLN_OPTIMISER_HELICAL_Z_PERCENTAGE);
    ADD_CONST(RLN_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER);
    ADD_CONST(RLN_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER);
    ADD_CONST(RLN_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT);
    ADD_CONST(RLN_OPTIMISER_HELICAL_SIGMA_DISTANCE);
    ADD_CONST(RLN_OPTIMISER_IGNORE_HELICAL_SYMMETRY);
    ADD_CONST(RLN_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED);
    ADD_CONST(RLN_PARTICLE_HELICAL_TUBE_ID);
    ADD_CONST(RLN_PARTICLE_HELICAL_TUBE_PITCH);
    ADD_CONST(RLN_PARTICLE_HELICAL_TRACK_LENGTH);
    ADD_CONST(RLN_SAMPLING_HELICAL_OFFSET_STEP);

    // New SGD labels in Relion 2.1
    ADD_CONST(RLN_MLMODEL_SGD_GRADIENT_IMAGE);
    ADD_CONST(RLN_OPTIMISER_DO_SGD);
    ADD_CONST(RLN_OPTIMISER_SGD_MU);
    ADD_CONST(RLN_OPTIMISER_SGD_SIGMA2FUDGE_INI);
    ADD_CONST(RLN_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE);
    ADD_CONST(RLN_OPTIMISER_SGD_SUBSET_START);
    ADD_CONST(RLN_OPTIMISER_SGD_SUBSET_SIZE);
    ADD_CONST(RLN_OPTIMISER_SGD_WRITE_EVERY_SUBSET);
    ADD_CONST(RLN_OPTIMISER_SGD_MAX_SUBSETS);
    ADD_CONST(RLN_OPTIMISER_SGD_STEPSIZE);
    ADD_CONST(RLN_OPTIMISER_HIGHRES_LIMIT_SGD);

    ADD_CONST(RLN_POSTPROCESS_BFACTOR);
    ADD_CONST(RLN_POSTPROCESS_FINAL_RESOLUTION);
    ADD_CONST(RLN_POSTPROCESS_FSC_TRUE);
    ADD_CONST(RLN_POSTPROCESS_FSC_MASKED);
    ADD_CONST(RLN_POSTPROCESS_FSC_UNMASKED);
    ADD_CONST(RLN_POSTPROCESS_FSC_RANDOM_MASKED);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_FIT_CORRELATION);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_FIT_SLOPE);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_VALUE_IN);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_VALUE_INVMTF);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT);
    ADD_CONST(RLN_POSTPROCESS_GUINIER_RESOL_SQUARED);
    ADD_CONST(RLN_POSTPROCESS_MTF_VALUE); ///< Detector MTF value
    ADD_CONST(RLN_SAMPLING_IS_3D);
    ADD_CONST(RLN_SAMPLING_IS_3D_TRANS);
    ADD_CONST(RLN_SAMPLING_HEALPIX_ORDER);
    ADD_CONST(RLN_SAMPLING_LIMIT_TILT);
    ADD_CONST(RLN_SAMPLING_OFFSET_RANGE);
    ADD_CONST(RLN_SAMPLING_OFFSET_STEP);
    ADD_CONST(RLN_SAMPLING_PERTURB);
    ADD_CONST(RLN_SAMPLING_PERTURBATION_FACTOR);
    ADD_CONST(RLN_SAMPLING_PRIOR_MODE);
    ADD_CONST(RLN_SAMPLING_PSI_STEP);
    ADD_CONST(RLN_SAMPLING_SIGMA_ROT);
    ADD_CONST(RLN_SAMPLING_SIGMA_TILT);
    ADD_CONST(RLN_SAMPLING_SIGMA_PSI);
    ADD_CONST(RLN_SAMPLING_SYMMETRY);

    ADD_CONST(RLN_SELECTED);
    ADD_CONST(RLN_SELECT_PARTICLES_ZSCORE);
    ADD_CONST(RLN_SORTED_IDX);
    ADD_CONST(RLN_PERFRAME_CUMULATIVE_WEIGHT);
    ADD_CONST(RLN_PERFRAME_RELATIVE_WEIGHT);

    ADD_CONST(RLN_RESOLUTION);
    ADD_CONST(RLN_RESOLUTION_ANGSTROM);
    ADD_CONST(RLN_RESOLUTION_INVPIXEL);
    ADD_CONST(RLN_SPECTRAL_IDX);
   
}
