/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <queue>
#include <list>
#include "data/fourier_filter.h"
#include "core/histogram.h"
#include "core/metadata_extension.h"
#include "core/multidim_array.h"
#include "core/transformations.h"
#include "core/xmipp_program.h"
#include "filters.h"
#include "mask.h"
#include "morphology.h"
#include "wavelet.h"

/* Subtract background ---------------------------------------------------- */
void substractBackgroundPlane(MultidimArray<double> &I)
{

    I.checkDimension(2);

    Matrix2D<double> A(3, 3);
    Matrix1D<double> x(3);
    Matrix1D<double> b(3);

    // Solve the plane 'x'
    A.initZeros();
    b.initZeros();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    {
        A(0, 0) += j * j;
        A(0, 1) += j * i;
        A(0, 2) += j;
        A(1, 1) += i * i;
        A(1, 2) += i;
        A(2, 2) += 1;
        b(0) += j * A2D_ELEM(I, i, j);
        b(1) += i * A2D_ELEM(I, i, j);
        b(2) += A2D_ELEM(I, i, j);
    }
    A(1, 0) = A(0, 1);
    A(2, 0) = A(0, 2);
    A(2, 1) = A(1, 2);
    solve(A, b, x);

    // Now subtract the plane
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    A2D_ELEM(I, i, j) -= x(0) * i + x(1) * j + x(2);
}

/* Subtract background ---------------------------------------------------- */
void substractBackgroundRollingBall(MultidimArray<double> &I, int radius)
{

    I.checkDimension(2);

    // Build the ball
    int arcTrimPer;
    int shrinkFactor;
    if (radius <= 10)
    {
        shrinkFactor = 1;
        arcTrimPer = 24; // trim 24% in x and y
    }
    else if (radius <= 30)
    {
        shrinkFactor = 2;
        arcTrimPer = 24; // trim 24% in x and y
    }
    else if (radius <= 100)
    {
        shrinkFactor = 4;
        arcTrimPer = 32; // trim 32% in x and y
    }
    else
    {
        shrinkFactor = 8;
        arcTrimPer = 40; // trim 40% in x and y
    }

    double smallballradius = radius / shrinkFactor;
    if (smallballradius < 1)
        smallballradius = 1;
    double r2 = smallballradius * smallballradius;
    int xtrim = (int) (arcTrimPer * smallballradius) / 100; // only use a patch of the rolling ball
    int halfWidth = ROUND(smallballradius - xtrim);
    int ballWidth = 2 * halfWidth + 1;
    MultidimArray<double> ball(ballWidth, ballWidth);
    ball.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ball)
    {
        double temp = r2 - i * i - j * j;
        ball(i, j) = temp > 0. ? sqrt(temp) : 0.;
    }

    // Shrink the image: each point in the shrinked image is the
    // minimum of its neighbourhood
    int sXdim = (XSIZE(I) + shrinkFactor - 1) / shrinkFactor;
    int sYdim = (YSIZE(I) + shrinkFactor - 1) / shrinkFactor;
    MultidimArray<double> shrinkI(sYdim, sXdim);
    shrinkI.setXmippOrigin();
    for (int ySmall = 0; ySmall < sYdim; ySmall++)
    {
        for (int xSmall = 0; xSmall < sXdim; xSmall++)
        {
            double minVal = 1e38;
            for (int j = 0, y = shrinkFactor * ySmall;
                 j < shrinkFactor && y < (int)YSIZE(I); j++, y++)
                for (int k = 0, x = shrinkFactor * xSmall;
                     k < shrinkFactor && x < (int)XSIZE(I); k++, x++)
                {
                    double thispixel = DIRECT_A2D_ELEM(I,y,x);
                    if (thispixel < minVal)
                        minVal = thispixel;
                }
            DIRECT_A2D_ELEM(shrinkI,ySmall,xSmall) = minVal;
        }
    }

    // Now roll the ball
    radius = ballWidth / 2;
    MultidimArray<double> Irolled;
    Irolled.resizeNoCopy(shrinkI);
    Irolled.initConstant(-500);
    for (int yb = -radius; yb < (int)YSIZE(shrinkI) + radius; yb++)
    {
        // Limits of the ball
        int y0 = yb - radius;
        if (y0 < 0)
            y0 = 0;
        int y0b = y0 - yb + radius; //y coordinate in the ball corresponding to y0
        int yF = yb + radius;
        if (yF >= (int)YSIZE(shrinkI))
            yF = (int)YSIZE(shrinkI) - 1;

        for (int xb = -radius; xb < (int)XSIZE(shrinkI) + radius; xb++)
        {
            // Limits of the ball
            int x0 = xb - radius;
            if (x0 < 0)
                x0 = 0;
            int x0b = x0 - xb + radius;
            int xF = xb + radius;
            if (xF >= (int)XSIZE(shrinkI))
                xF = (int)XSIZE(shrinkI) - 1;

            double z = 1e38;
            for (int yp = y0, ybp = y0b; yp <= yF; yp++, ybp++)
                for (int xp = x0, xbp = x0b; xp <= xF; xp++, xbp++)
                {
                    double zReduced = DIRECT_A2D_ELEM(shrinkI,yp,xp)
                                      - DIRECT_A2D_ELEM(ball,ybp,xbp);
                    if (z > zReduced)
                        z = zReduced;
                }
            for (int yp = y0, ybp = y0b; yp <= yF; yp++, ybp++)
                for (int xp = x0, xbp = x0b; xp <= xF; xp++, xbp++)
                {
                    double zMin = z + DIRECT_A2D_ELEM(ball,ybp,xbp);
                    if (DIRECT_A2D_ELEM(Irolled,yp,xp) < zMin)
                        DIRECT_A2D_ELEM(Irolled,yp,xp) = zMin;
                }
        }
    }

    // Now rescale the background
    MultidimArray<double> bgEnlarged;
    scaleToSize(xmipp_transformation::LINEAR, bgEnlarged, Irolled, XSIZE(I), YSIZE(I));
    bgEnlarged.copyShape(I);
    I -= bgEnlarged;
}

/* Detect background ------------------------------------------------------ */
void detectBackground(const MultidimArray<double> &vol,
                      MultidimArray<double> &mask, double alpha, double &final_mean)
{

    // 2.1.-Background detection------------------------------------------
    MultidimArray<double> bg; // We create the volumen with
    bg.resizeNoCopy(vol); // -1:not visited 0:mol 1:background
    bg.initConstant(-1); // -2:in the list

    // Ponemos las seis caras de esta variable como visitadas e inicializamos
    // la cola de pixeles por visitar
    std::queue<int> list_for_compute; // Lista del modo [x1,y1,z1,...,xi,yi,zi]
    // que contiene los pixeles por procesar
    std::vector<double> bg_values; // Vector con los valores del background
    size_t xdim = XSIZE(bg);
    size_t ydim = YSIZE(bg);
    size_t zdim = ZSIZE(bg);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(bg)
    {
        if (j == 0 || j == xdim - 1 || i == 0 || i == ydim - 1 || k == 0
            || k == zdim - 1)
        { // Visited coord
            DIRECT_A3D_ELEM(bg,k,i,j) = 1;
            // We introduce the values of the background
            bg_values.push_back(DIRECT_A3D_ELEM(vol,k,i,j));

        }
        if ((j == 1 || j == xdim - 2) && (i != 0) && (i != ydim - 1) && (k != 0)
            && (k != zdim - 1))
        { // Coord for compute for x
            list_for_compute.push(j);
            list_for_compute.push(i);
            list_for_compute.push(k);
        }
        if ((i == 1 || i == ydim - 2) && (j > 1) && (j < xdim - 2) && (k != 0)
            && (k != zdim - 1))
        { // Coord for compute for y
            list_for_compute.push(j);
            list_for_compute.push(i);
            list_for_compute.push(k);
        }
        if ((k == 1 || k == zdim - 2) && (j > 1) && (j < xdim - 2) && (i > 1)
            && (i < ydim - 2))
        { // Coord for compute for y
            list_for_compute.push(j);
            list_for_compute.push(i);
            list_for_compute.push(k);
        }
    } // end of FOR_ALL_ELEMENTS

    // We work until the list_for_compute is empty
    int n = 250; //each 250 pixels renew stats
    int cont = 250; //We start here for compute stat for first time
    double A=0; // A and B are numbers such the interval of confidence is [A,B]
    double B=0; //
    float z = icdf_gauss(1 - alpha / 2);
    while (!list_for_compute.empty())
    {

        //We compute stat when is needed
        if (cont == n)
        {
            // Compute statistics
            double avg=0;
            double stddev=0;
            computeAvgStddev(bg_values, avg, stddev);
            final_mean = avg;
            // Compute confidence interval
            A = avg - (z * stddev);
            B = avg + (z * stddev);
            cont = 0;
        } // end of if
        // Now we start to take coords from the list_for_compute
        int x_coord = list_for_compute.front();
        list_for_compute.pop();
        int y_coord = list_for_compute.front();
        list_for_compute.pop();
        int z_coord = list_for_compute.front();
        list_for_compute.pop();
        // Is visited
        DIRECT_A3D_ELEM(bg,z_coord,y_coord,x_coord) = -2;
        //We take the value
        double value = DIRECT_A3D_ELEM(vol,z_coord,y_coord,x_coord);
        // We see if is background or not
        if (A <= value && value <= B)
        {
            // We now is background
            DIRECT_A3D_ELEM(bg,z_coord,y_coord,x_coord) = 1;
            // We introduce the values of the background
            bg_values.push_back(value);
            // We update the cont variable
            cont++;

            // We add his neighbours in the list_for_compute
            for (int xx = x_coord - 1; xx <= x_coord + 1; xx++)
                for (int yy = y_coord - 1; yy <= y_coord + 1; yy++)
                    for (int zz = z_coord - 1; zz <= z_coord + 1; zz++)
                        //We see if it has been visited
                        if (DIRECT_A3D_ELEM(bg,zz,yy,xx) == -1) // not visited
                        {
                            // So we include it in the list
                            list_for_compute.push(xx);
                            list_for_compute.push(yy);
                            list_for_compute.push(zz);
                            // Is in the list
                            DIRECT_A3D_ELEM(bg,zz,yy,xx) = -2;
                        }
        } // end of if
        else
        {
            // Isn't background
            DIRECT_A3D_ELEM(bg,z_coord,y_coord,x_coord) = 0;
        }
    } // end of while
    // Now we change not visited for mol
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(bg)
    if (DIRECT_A3D_ELEM(bg,k,i,j) == -1)
        DIRECT_A3D_ELEM(bg,k,i,j) = 0;

    // End of 2.1-----------------------------------------------------------
    // 2.2.-Matematical Morphology
    MultidimArray<double> bg_mm; // We create the output volumen
    bg_mm.initZeros(vol);
    closing3D(bg, bg_mm, 26, 0, 1);
    // Output
    //typeCast(bg_mm,mask);
    mask = bg_mm;
}

/* Contranst enhancement --------------------------------------------------- */
void contrastEnhancement(Image<double> *I)
{
    (*I)().rangeAdjust(0, 255);
}

/* Region growing for images ----------------------------------------------- */
typedef struct
{
    int ii;
    int jj;
}
Coordinate2D;

void regionGrowing2D(const MultidimArray<double> &I_in,
                     MultidimArray<double> &I_out, int i, int j, float stop_colour,
                     float filling_colour, bool less, int neighbourhood)
{
    I_in.checkDimension(2);

    std::list<Coordinate2D> iNeighbours; /* A list for neighbour pixels */
    int iCurrenti;
    int iCurrentj; /* Coordinates of the current pixel considered */

    /* First task is copying the input image into the output one */
    I_out = I_in;

    /**** Then the region growing is done ****/
    /* Insert at the beginning of the list the seed coordinates */
    Coordinate2D coord;
    coord.ii = i;
    coord.jj = j;
    iNeighbours.push_front(coord);

    /* Fill the seed coordinates */
    A2D_ELEM(I_out, i, j) = filling_colour;

    while (!iNeighbours.empty())
    {
        /* Take the current pixel to explore */
        coord = iNeighbours.front();
        iNeighbours.pop_front();
        iCurrenti = coord.ii;
        iCurrentj = coord.jj;

#define CHECK_POINT(i,j) \
    { \
  int I=i; \
        int J=j; \
  if (INSIDEXY(I_out,J,I))  { \
   double pixel=A2D_ELEM(I_out,I,J);\
   if (pixel!=filling_colour) \
    if ((less && pixel < stop_colour) || \
     (!less && pixel > stop_colour)) { \
     coord.ii=I; \
     coord.jj=J; \
     A2D_ELEM (I_out,coord.ii,coord.jj)=filling_colour; \
     iNeighbours.push_front(coord); \
    } \
  } \
    }

        /* Make the exploration of the pixel's neighbours */
        CHECK_POINT(iCurrenti, iCurrentj - 1);
        CHECK_POINT(iCurrenti, iCurrentj + 1);
        CHECK_POINT(iCurrenti - 1, iCurrentj);
        CHECK_POINT(iCurrenti + 1, iCurrentj);
        if (neighbourhood == 8)
        {
            CHECK_POINT(iCurrenti - 1, iCurrentj - 1);
            CHECK_POINT(iCurrenti - 1, iCurrentj + 1);
            CHECK_POINT(iCurrenti + 1, iCurrentj - 1);
            CHECK_POINT(iCurrenti + 1, iCurrentj + 1);
        }
    }
}

/* Region growing for volumes ----------------------------------------------- */
typedef struct
{
    int ii;
    int jj;
    int kk;
}
Coordinate3D;

void regionGrowing3D(const MultidimArray<double> &V_in,
                     MultidimArray<double> &V_out, int k, int i, int j, float stop_colour,
                     float filling_colour, bool less)
{
    V_in.checkDimension(3);

    std::list<Coordinate3D> iNeighbours; /* A list for neighbour voxels */
    int iCurrentk;
    int iCurrenti;
    int iCurrentj; /* Coordinates of the current voxel considered */

    /* First task is copying the input volume into the output one */
    V_out = V_in;

    /**** Then the region growing is done in output volume ****/
    /* Insert at the beginning of the list the seed coordinates */
    Coordinate3D coord;
    coord.ii = i;
    coord.jj = j;
    coord.kk = k;
    iNeighbours.push_front(coord);

    /* Fill the seed coordinates */
    A3D_ELEM(V_out, k, i, j) = filling_colour;

    while (!iNeighbours.empty())
    {
        /* Take the current pixel to explore */
        coord = iNeighbours.front();
        iNeighbours.pop_front();
        iCurrenti = coord.ii;
        iCurrentj = coord.jj;
        iCurrentk = coord.kk;

        /* a macro for doing exploration of a voxel. If the voxel has a value
         lower than stop_colour, its filled with filling colour and added to the
         list for exploring its neighbours */
#define CHECK_POINT_3D(k,i,j) \
 {\
     int I=i; \
        int J=j; \
        int K=k; \
  if (INSIDEXYZ(V_out,J,I,K))  { \
   double voxel=A3D_ELEM(V_out,K,I,J); \
   if (voxel!=filling_colour) \
    if ((less && voxel < stop_colour)|| \
     (!less &&voxel > stop_colour)) { \
     coord.ii=I; \
     coord.jj=J; \
     coord.kk=K; \
     A3D_ELEM (V_out,coord.kk,coord.ii,coord.jj)=filling_colour; \
     iNeighbours.push_front(coord); \
    } \
  }\
 }

        /* Make the exploration of the pixel neighbours */
        CHECK_POINT_3D(iCurrentk, iCurrenti, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti - 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk, iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti + 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk, iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti, iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti, iCurrentj);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj + 1);
    }
}


void regionGrowing3DEqualValue(const MultidimArray<double> &V_in, MultidimArray<int> &V_out, int filling_value=0)
{
    V_in.checkDimension(3);

    std::list<Coordinate3D> iNeighbours; /* A list for neighbour voxels */
    int iCurrentk;
    int iCurrenti;
    int iCurrentj; /* Coordinates of the current voxel considered */

    /* First task is copying the input volume into the output one */
    V_out.resizeNoCopy(V_in);
    V_out.initConstant(1);

    /**** Then the region growing is done in output volume ****/
    /* Insert at the beginning of the list the seed coordinates */
    Coordinate3D coord;
    coord.ii = STARTINGY(V_in);
    coord.jj = STARTINGX(V_in);
    coord.kk = STARTINGZ(V_in);
    iNeighbours.push_front(coord);

    /* Fill the seed coordinates */
    double bgColor = A3D_ELEM(V_in, STARTINGZ(V_in), STARTINGY(V_in), STARTINGX(V_in));
    A3D_ELEM(V_out, STARTINGZ(V_out), STARTINGY(V_out), STARTINGX(V_out)) = 0;

    while (!iNeighbours.empty())
    {
        /* Take the current pixel to explore */
        coord = iNeighbours.front();
        iNeighbours.pop_front();
        iCurrenti = coord.ii;
        iCurrentj = coord.jj;
        iCurrentk = coord.kk;

        /* a macro for doing exploration of a voxel. If the voxel has a value
         lower than stop_colour, its filled with filling colour and added to the
         list for exploring its neighbours */
#undef CHECK_POINT_3D
#define CHECK_POINT_3D(k,i,j) \
 {\
        int I=i; \
        int J=j; \
        int K=k; \
  if (INSIDEXYZ(V_out,J,I,K))  { \
   if (A3D_ELEM(V_in,K,I,J)==bgColor && A3D_ELEM(V_out,K,I,J)==1) { \
     coord.ii=I; \
     coord.jj=J; \
     coord.kk=K; \
     A3D_ELEM (V_out,coord.kk,coord.ii,coord.jj)=filling_value; \
     iNeighbours.push_front(coord); \
    } \
  }\
 }\

        /* Make the exploration of the pixel neighbours */
        CHECK_POINT_3D(iCurrentk, iCurrenti, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti - 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk, iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti + 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk, iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk, iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti, iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti, iCurrentj);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj + 1);
    }
}


void distanceTransform(const MultidimArray<int> &in, MultidimArray<int> &out,
                       bool wrap)
{
    std::list<int> toExplore; /* A list of points to explore */

    in.checkDimension(2);

    out.resize(in);
    out.initConstant(XSIZE(in) + YSIZE(in));

#define CHECK_POINT_DIST(i,j,proposedDistance) \
    { \
        int ip=i; \
        int jp=j; \
        if (wrap) \
        { \
            ip=intWRAP(ip,STARTINGY(out),FINISHINGY(out));\
            jp=intWRAP(jp,STARTINGX(out),FINISHINGX(out));\
        } \
        if (ip>=STARTINGY(out) && ip<=FINISHINGY(out) && \
            jp>=STARTINGX(out) && jp<=FINISHINGX(out)) \
            if (out(ip,jp)>proposedDistance) \
            { \
                out(ip,jp)=proposedDistance; \
                toExplore.push_back(ip); \
                toExplore.push_back(jp); \
                toExplore.push_back(proposedDistance); \
            } \
    }

    // Look for all elements in the binary mask and set the corresponding
    // distance to 0
    FOR_ALL_ELEMENTS_IN_ARRAY2D(in)
    if (in(i, j))
    {
        out(i, j) = 0;
        CHECK_POINT_DIST(i-1, j, 1);
        CHECK_POINT_DIST(i+1, j, 1);
        CHECK_POINT_DIST(i, j-1, 1);
        CHECK_POINT_DIST(i, j+1, 1);
    }

    while (!toExplore.empty())
    {
        int i = toExplore.front();
        toExplore.pop_front();
        int j = toExplore.front();
        toExplore.pop_front();
        int proposedDistance = toExplore.front();
        toExplore.pop_front();

        if (proposedDistance == out(i, j))
        {
            // If this is the current distance (i.e., it has not
            // been supersceded by a later value
            CHECK_POINT_DIST(i-1, j, proposedDistance+1);
            CHECK_POINT_DIST(i+1, j, proposedDistance+1);
            CHECK_POINT_DIST(i, j-1, proposedDistance+1);
            CHECK_POINT_DIST(i, j+1, proposedDistance+1);
        }
    }
}

/* Label image ------------------------------------------------------------ */
int labelImage2D(const MultidimArray<double> &I, MultidimArray<double> &label,
                 int neighbourhood)
{
    I.checkDimension(2);

    label = I;
    int colour = 32000;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(label)
    {
        if (label(i, j) != 1)
            continue;
        regionGrowing2D(label, label, i, j, 0, colour, false, neighbourhood);
        colour++;
    }
    FOR_ALL_ELEMENTS_IN_ARRAY2D(label)
    if (label(i, j) != 0)
        label(i, j) = label(i, j) - 31999;
    return colour - 32000;
}

/* Label volume ------------------------------------------------------------ */
int labelImage3D(const MultidimArray<double> &V, MultidimArray<double> &label)
{
    V.checkDimension(3);

    label = V;
    int colour = 32000;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(label)
    {
        if (label(k, i, j) != 1)
            continue;
        regionGrowing3D(label, label, k, i, j, 0, colour, false);
        colour++;
    }
    FOR_ALL_ELEMENTS_IN_ARRAY3D(label)
    if (A3D_ELEM(label,k, i, j) != 0)
        A3D_ELEM(label,k, i, j) = A3D_ELEM(label,k, i, j) - 31999;
    return colour - 32000;
}

/* Remove small components ------------------------------------------------- */
void removeSmallComponents(MultidimArray<double> &I, int size,
                           int neighbourhood)
{
    MultidimArray<double> label;
    int imax;
    if (ZSIZE(I)==1)
    	imax=labelImage2D(I, label, neighbourhood);
    else
    	imax=labelImage3D(I, label);
    MultidimArray<int> nlabel(imax + 1);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(label)
    {
    	int l=(int)DIRECT_MULTIDIM_ELEM(label,n);
    	DIRECT_A1D_ELEM(nlabel,l)++;
    }
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(label)
    {
    	int l=(int)DIRECT_MULTIDIM_ELEM(label,n);
    	if (DIRECT_A1D_ELEM(nlabel,l)<size)
    		DIRECT_MULTIDIM_ELEM(I,n)=0;
    }
}

/* Keep biggest component -------------------------------------------------- */
void keepBiggestComponent(MultidimArray<double> &I, double percentage,
                          int neighbourhood)
{
    MultidimArray<double> label;
    int imax;
    if (ZSIZE(I)==1)
    	imax=labelImage2D(I, label, neighbourhood);
    else
    	imax=labelImage3D(I, label);
    MultidimArray<int> nlabel(imax + 1);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(label)
    {
    	int l=(int)DIRECT_MULTIDIM_ELEM(label,n);
    	if (l>0)
    		A1D_ELEM(nlabel,l)++;
    }

    MultidimArray <int> best;
    nlabel.indexSort(best);
    best -= 1;
    int nbest = XSIZE(best) - 1;
    double total = nlabel.sum();
    double explained = nlabel(best(nbest));
    while (explained < percentage * total)
    {
        nbest--;
        explained += nlabel(best(nbest));
    }

    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(label)
    {
    	int l=(int)DIRECT_MULTIDIM_ELEM(label,n);
        bool among_the_best = false;
        for (int k = nbest; k < imax + 1; k++)
            if (l == A1D_ELEM(best,k))
            {
                among_the_best = true;
                break;
            }
        if (!among_the_best)
    		DIRECT_MULTIDIM_ELEM(I,n)=0;
    }
}

/* Fill object ------------------------------------------------------------- */
void fillBinaryObject(MultidimArray<double> &I, int neighbourhood)
{
    I.checkDimension(2);

    MultidimArray<double> label;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    I(i, j) = 1 - I(i, j);
    labelImage2D(I, label, neighbourhood);
    double l0 = label(STARTINGY(I), STARTINGX(I));
    FOR_ALL_ELEMENTS_IN_ARRAY2D(label)
    if (label(i, j) == l0)
        I(i, j) = 0;
    else
        I(i, j) = 1;
}

/* Variance filter ----------------------------------------------------------*/
void varianceFilter(MultidimArray<double> &I, int kernelSize, bool relative)
{
    int kernelSize_2 = kernelSize/2;
    MultidimArray<double> kernel;
    kernel.resize(kernelSize,kernelSize);
    kernel.setXmippOrigin();

    // std::cout << " Creating the variance matrix " << std::endl;
    MultidimArray<double> mVar(YSIZE(I),XSIZE(I));
    mVar.setXmippOrigin();
    double stdKernel;
    double avgKernel;
    double min_val;
    double max_val;
    int x0;
    int y0;
    int xF;
    int yF;

    // I.computeStats(avgImg, stdImg, min_im, max_im);    
    
    for (int i=kernelSize_2; i<=(int)YSIZE(I)-kernelSize_2; i+=kernelSize_2)
        for (int j=kernelSize_2; j<=(int)XSIZE(I)-kernelSize_2; j+=kernelSize_2)
            {
                x0 = j-kernelSize_2;
                y0 = i-kernelSize_2;
                xF = j+kernelSize_2-1;
                yF = i+kernelSize_2-1;

                if (x0 < 0)
                    x0 = 0;
                if (xF > XSIZE(I))
                    xF = XSIZE(I);
                if (y0 < 0)
                    y0 = 0;
                if (yF > YSIZE(I))
                    yF = YSIZE(I);

                I.window(kernel, y0, x0, yF, xF);
                kernel.computeStats(avgKernel, stdKernel, min_val, max_val);
                
                DIRECT_A2D_ELEM(mVar, i, j) = stdKernel;
            }

    // filtering to fill the matrices (convolving with a Gaussian)
    FourierFilter filter;
    filter.FilterShape = REALGAUSSIAN;
    filter.FilterBand = LOWPASS;
    filter.w1 = kernelSize_2;
    filter.applyMaskSpace(mVar);

    if (relative) // normalize to the mean of the variance
    {    
        double avgVar;
        double stdVar;
        double minVar;
        double maxVar;
        mVar.computeStats(avgVar, stdVar, minVar, maxVar);
        mVar = mVar/avgVar;
    }

    I = mVar;

    // filling the borders with the nearest variance value
    // the corners only are filled once (with Ycoord)
    for (int i=0; i<YSIZE(I); ++i)
        for (int j=0; j<kernelSize; j++)
        {
            if (i<kernelSize)
                DIRECT_A2D_ELEM(I, i, j) = DIRECT_A2D_ELEM(mVar, kernelSize, kernelSize);
            else if (i>YSIZE(I)-kernelSize)
                DIRECT_A2D_ELEM(I, i, j) = DIRECT_A2D_ELEM(mVar, YSIZE(I)-kernelSize, kernelSize);
            else
                DIRECT_A2D_ELEM(I, i, j) = DIRECT_A2D_ELEM(mVar, i, kernelSize);
        }
    for (int i=0; i<YSIZE(I); ++i)
        for (int j=XSIZE(I)-kernelSize; j<XSIZE(I); j++)
        {
            if (i<kernelSize)
                DIRECT_A2D_ELEM(I, i, j) = DIRECT_A2D_ELEM(mVar, kernelSize, XSIZE(I)-kernelSize);
            else if (i>YSIZE(I)-kernelSize)
                DIRECT_A2D_ELEM(I, i, j) = DIRECT_A2D_ELEM(mVar, YSIZE(I)-kernelSize, XSIZE(I)-kernelSize);
            else
                DIRECT_A2D_ELEM(I, i, j) = DIRECT_A2D_ELEM(mVar, i, XSIZE(I)-kernelSize);
        }
    for (int j=kernelSize; j<XSIZE(I)-kernelSize; ++j)
        for (int i=0; i<kernelSize; i++)
            DIRECT_A2D_ELEM(I, i, j) = DIRECT_A2D_ELEM(mVar, kernelSize, j);
    for (int j=kernelSize; j<XSIZE(I)-kernelSize; ++j)
        for (int i=YSIZE(I)-kernelSize; i<YSIZE(I); i++)
            DIRECT_A2D_ELEM(I, i, j) = DIRECT_A2D_ELEM(mVar, YSIZE(I)-kernelSize, j);
    
}

/* Noise filter (returns a binary mask where both variance and mean are high)*/
void noisyZonesFilter(MultidimArray<double> &I, int kernelSize)
{
    int kernelSize_2 = kernelSize/2;
    MultidimArray<double> kernel;
    kernel.resize(kernelSize,kernelSize);
    kernel.setXmippOrigin();

    MultidimArray<double> mAvg=I;
    MultidimArray<double> mVar=I;
    double stdKernel;
    double varKernel;
    double avgKernel;
    double min_val;
    double max_val;
    int x0;
    int y0;
    int xF;
    int yF;
    
    for (int i=kernelSize_2; i<(int)YSIZE(I); i+=kernelSize_2)
        for (int j=kernelSize_2; j<(int)XSIZE(I); j+=kernelSize_2)
            {
                x0 = j-kernelSize_2;
                y0 = i-kernelSize_2;
                xF = j+kernelSize_2-1;
                yF = i+kernelSize_2-1;

                if (x0 < 0)
                    x0 = 0;
                if (xF > XSIZE(I))
                    xF = XSIZE(I);
                if (y0 < 0)
                    y0 = 0;
                if (yF > YSIZE(I))
                    yF = YSIZE(I);

                I.window(kernel, y0, x0, yF, xF);
                kernel.computeStats(avgKernel, stdKernel, min_val, max_val);
                varKernel = stdKernel*stdKernel;

                DIRECT_A2D_ELEM(mAvg, i, j) = avgKernel*avgKernel;
                DIRECT_A2D_ELEM(mVar, i, j) = varKernel;
            }

    // filtering to fill the matrices (convolving with a Gaussian)
    FourierFilter filter;
    filter.FilterShape = REALGAUSSIAN;
    filter.FilterBand = LOWPASS;
    filter.w1 = kernelSize_2;
    filter.applyMaskSpace(mAvg);
    filter.applyMaskSpace(mVar);

    // // Draw to debug
    // Image<double> imAvg(mAvg), imVar(mVar);
    // imAvg.write("AvgFilter.mrc");
    // imVar.write("VarFilter.mrc");

    // Working in a auxilary windows to avoid borders bad defined
    auto ysize = static_cast<int>YSIZE(I);
    auto xsize = static_cast<int>XSIZE(I);
    MultidimArray<double> mAvgAux(ysize-kernelSize,xsize-kernelSize);
    MultidimArray<double> mVarAux(ysize-kernelSize,xsize-kernelSize);
    mAvgAux.setXmippOrigin();
    mVarAux.setXmippOrigin();
    mAvg.window(mAvgAux,STARTINGY(mAvg)+kernelSize_2, STARTINGX(mAvg)+kernelSize_2,
                       FINISHINGY(mAvg)-kernelSize_2, FINISHINGX(mAvg)-kernelSize_2);
    mVar.window(mVarAux,STARTINGY(mVar)+kernelSize_2, STARTINGX(mVar)+kernelSize_2,
                       FINISHINGY(mVar)-kernelSize_2, FINISHINGX(mVar)-kernelSize_2);

    // Refiltering to get a smoother distribution
    // filter.w1 = XSIZE(I)/40;
    // filter.applyMaskSpace(mAvgAux);
    // filter.applyMaskSpace(mVarAux);

    // Binarization
    MultidimArray<double> mAvgAuxBin = mAvgAux;
    MultidimArray<double> mVarAuxBin = mVarAux;
    EntropySegmentation(mVarAuxBin);
    // EntropySegmentation(mAvgAuxBin);
    float th = EntropySegmentation(mAvgAuxBin);
    mAvgAuxBin.binarize(th*0.92);
    mAvgAuxBin = 1-mAvgAuxBin;
    // std::cout << "binarize threshold = " << th << std::endl;

    // Returning to the previous windows size
    MultidimArray<double> mAvgBin = mAvg;
    MultidimArray<double> mVarBin = mVar;
    mAvgAuxBin.window(mAvgBin, STARTINGY(mVar), STARTINGX(mVar),
                              FINISHINGY(mVar), FINISHINGX(mVar));
    mVarAuxBin.window(mVarBin, STARTINGY(mVar), STARTINGX(mVar),
                              FINISHINGY(mVar), FINISHINGX(mVar));

    // // Draw to debug
    // Image<double> imAvgBin(mAvgBin), imVarBin(mVarBin);
    // imAvgBin.write("noisyZoneFilter_AVGmask.mrc");
    // imVarBin.write("noisyZoneFilter_VARmask.mrc");

    // Combining both masks
    I = 1-(mVarBin*mAvgBin);
}

/* Gini Coefficient -- (applies a variance filter to the input Image) ------ */
double giniCoeff(MultidimArray<double> &I, int varKernelSize)
{
    MultidimArray<double> im = I;

    // std::cout << " - Starting fft filtering " << std::endl;

    FourierFilter filter;
    filter.FilterShape = REALGAUSSIAN;
    filter.FilterBand = LOWPASS;
    filter.w1 = 4;
    filter.applyMaskSpace(im);
   
    // Image<double> imG(im);
    // imG.write("I_Gauss.mrc");

    // std::cout << " - Calling varianceFilter() " << std::endl;
    varianceFilter(I, varKernelSize, true);
    im = I;

    // std::cout << " - Starting 2nd fft filtering " << std::endl;
    filter.w1 = varKernelSize/8;
    filter.applyMaskSpace(I);

    // Image<double> imGV(im);
    // imGV.write("I_Gauss_Var.mrc");
   
    im -= im.computeMin();
    im /= im.computeMax();

    // Image<double> imGVN(im);
    // imGVN.write("I_Gauss_Var_Norm.mrc");
   
    // std::cout << " - Starting histogram analysis " << std::endl;
    Histogram1D hist;
    hist.clear();
    compute_hist(im, hist, 256);

    // std::cout << " - Computing Gini coeff " << std::endl;
    MultidimArray<double> sortedList=hist;
    hist.sort(sortedList);
    double height=0;
    double area=0;
    for (int i=0; i<XSIZE(sortedList); i++)
    {
        height += DIRECT_MULTIDIM_ELEM(sortedList,i);
        area += height - DIRECT_MULTIDIM_ELEM(sortedList,i)/2.0;
    }
        
    double fair_area = height*XSIZE(hist)/2.0;

    double giniValue = (fair_area-area)/fair_area;

    return giniValue;
}

/* Otsu Segmentation ------------------------------------------------------- */
double OtsuSegmentation(MultidimArray<double> &V)
{
    // V.checkDimension(3);

    // Compute the probability density function
    Histogram1D hist;
    hist.clear();
    compute_hist(V, hist, 200);
    hist /= hist.sum();

    // Compute the cumulative 0th and 1st order moments
    double x;
    MultidimArray<double> mom0;
    MultidimArray<double> mom1;
    mom0.initZeros(XSIZE(hist));
    mom1.initZeros(XSIZE(hist));
    mom0(0) = hist(0);
    hist.index2val(0, x);
    mom1(0) = hist(0) * x;
    for (size_t i = 1; i < XSIZE(mom0); i++)
    {
        mom0(i) = mom0(i - 1) + hist(i);
        hist.index2val(i, x);
        mom1(i) = mom1(i - 1) + hist(i) * x;
    }

    // Maximize sigma2B
    double bestSigma2B = -1;
    int ibestSigma2B = -1;
    for (size_t i = 0; i < XSIZE(hist) - 1; i++)
    {
        double w1 = mom0(i);
        double w2 = 1 - mom0(i);
        double mu1 = mom1(i);
        double mu2 = mom1(XSIZE(mom1) - 1) - mom1(i);
        double sigma2B = w1 * w2 * (mu1 - mu2) * (mu1 - mu2);
        if (sigma2B > bestSigma2B)
        {
            bestSigma2B = sigma2B;
            ibestSigma2B = i;
        }
    }

    hist.index2val(ibestSigma2B, x);
    V.binarize(x);

    return x;
}

/* Entropy Segmentation ---------------------------------------------------- */
double EntropySegmentation(MultidimArray<double> &V)
{
    // V.checkDimension(3);

    // Compute the probability density function
    Histogram1D hist;
    hist.clear();
    compute_hist(V, hist, 200);
    hist /= hist.sum();

    // Compute the cumulative 0th and 1st order moments
    double x;
    MultidimArray<double> mom0;
    mom0.initZeros(XSIZE(hist));
    mom0(0) = hist(0);
    for (size_t i = 1; i < XSIZE(mom0); i++)
        mom0(i) = mom0(i - 1) + hist(i);

    // Entropy for black and white parts of the histogram
    const double epsilon = 1e-15;
    MultidimArray<double> h1;
    MultidimArray<double> h2;
    h1.initZeros(XSIZE(hist));
    h2.initZeros(XSIZE(hist));
    for (size_t i = 0; i < XSIZE(hist); i++)
    {
        // Entropy h1
        double w1 = mom0(i);
        if (w1 > epsilon)
            for (size_t ii = 0; ii <= i; ii++)
                if (hist(ii) > epsilon)
                {
                    double aux = hist(ii) / w1;
                    h1(i) -= aux * log10(aux);
                }

        // Entropy h2
        double w2 = 1 - mom0(i);
        if (w2 > epsilon)
            for (size_t ii = i + 1; ii < XSIZE(hist); ii++)
                if (hist(ii) > epsilon)
                {
                    double aux = hist(ii) / w2;
                    h2(i) -= aux * log10(aux);
                }
    }

    // Find histogram index with maximum entropy
    double Hmax = h1(0) + h2(0);
    size_t iHmax = 0;
    for (size_t i = 1; i < XSIZE(hist) - 1; i++)
    {
        double H = h1(i) + h2(i);
        if (H > Hmax)
        {
            Hmax = H;
            iHmax = i;
        }
    }

    hist.index2val(iHmax, x);
    V.binarize(x);

    return x;
}

/* Otsu+Entropy Segmentation ----------------------------------------------- */
double EntropyOtsuSegmentation(MultidimArray<double> &V, double percentil,
                               bool binarizeVolume)
{
    //V.checkDimension(3);

    // Compute the probability density function
    Histogram1D hist;
    hist.clear();
    compute_hist(V, hist, 300);
    hist /= hist.sum();

    // Compute the cumulative 0th and 1st order moments
    double x;
    MultidimArray<double> mom0;
    MultidimArray<double> mom1;
    mom0.initZeros(XSIZE(hist));
    mom1.initZeros(XSIZE(hist));
    mom0(0) = hist(0);
    hist.index2val(0, x);
    mom1(0) = hist(0) * x;
    for (size_t i = 1; i < XSIZE(mom0); i++)
    {
        mom0(i) = mom0(i - 1) + hist(i);
        hist.index2val(i, x);
        mom1(i) = mom1(i - 1) + hist(i) * x;
    }

    // Entropy for black and white parts of the histogram
    const double epsilon = 1e-15;
    MultidimArray<double> h1;
    MultidimArray<double> h2;
    h1.initZeros(XSIZE(hist));
    h2.initZeros(XSIZE(hist));
    for (size_t i = 0; i < XSIZE(hist); i++)
    {
        // Entropy h1
        double w1 = mom0(i);
        if (w1 > epsilon)
            for (size_t ii = 0; ii <= i; ii++)
                if (hist(ii) > epsilon)
                {
                    double aux = hist(ii) / w1;
                    h1(i) -= aux * log10(aux);
                }

        // Entropy h2
        double w2 = 1 - mom0(i);
        if (w2 > epsilon)
            for (size_t ii = i + 1; ii < XSIZE(hist); ii++)
                if (hist(ii) > epsilon)
                {
                    double aux = hist(ii) / w2;
                    h2(i) -= aux * log10(aux);
                }
    }

    // Compute sigma2B and H
    MultidimArray<double> sigma2B;
    MultidimArray<double> H;
    MultidimArray<double> HSigma2B;
    sigma2B.initZeros(XSIZE(hist) - 1);
    H.initZeros(XSIZE(hist) - 1);
    HSigma2B.initZeros(XSIZE(hist) - 1);
    for (size_t i = 0; i < XSIZE(hist) - 1; i++)
    {
        double w1 = mom0(i);
        double w2 = 1 - mom0(i);
        double mu1 = mom1(i);
        double mu2 = mom1(XSIZE(mom1) - 1) - mom1(i);
        sigma2B(i) = w1 * w2 * (mu1 - mu2) * (mu1 - mu2);
        H(i) = h1(i) + h2(i);
        HSigma2B(i) = -log10(sigma2B(i)) / H(i);
        // The logic behind this expression is
        // Otsu:    max sigma2B -> max log10(sigma2B) -> min -log10(sigma2B)
        // Entropy: max H       -> max H              -> min 1/H
    }

    // Sort HSigma2B and take a given percentage of it
    MultidimArray<double> HSigma2Bsorted;
    HSigma2B.sort(HSigma2Bsorted);
    int iTh = ROUND(XSIZE(HSigma2B)*percentil);
    double threshold = HSigma2Bsorted(iTh);

    // Find the first value within HSigma2B falling below this threshold
    iTh = 0;
    while (A1D_ELEM(HSigma2B,iTh) >= threshold)
        iTh++;
    iTh--;
    if (iTh <= 0)
        x = threshold;
    else
        hist.index2val(iTh, x);

    if (binarizeVolume)
        V.binarize(x);
    return x;
}

/* Fast correntropy -------------------------------------------------------- */
double fastCorrentropy(const MultidimArray<double> &x,
                       const MultidimArray<double> &y, double sigma,
                       const GaussianInterpolator &G, const MultidimArray<int> &mask)
{
    double retvalxy = 0;
    double isigma = 1.0 / sigma;
    int maskSum = 0;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(x)
    {
        if (DIRECT_MULTIDIM_ELEM(mask,n))
        {
            retvalxy += G.getValue(
                            isigma
                            * (DIRECT_MULTIDIM_ELEM(x,n)
                               - DIRECT_MULTIDIM_ELEM(y,n)));
            ++maskSum;
        }
    }
    if (0 == maskSum) {
        REPORT_ERROR(ERR_NUMERICAL, "maskSum is zero (0), which would lead to division by zero");
    }
    return (retvalxy / maskSum);
}

/* Imed distance ---------------------------------------------------------- */
double imedDistance(const MultidimArray<double>& I1, const MultidimArray<double>& I2)
{
	// [x,y]=meshgrid([-3:1:3],[-3:1:3])
	// format long	// w=1/sqrt(2*pi)*exp(-0.5*(x.*x+y.*y))
	double *refW;
	double w[49]={
			0.000049233388666,   0.000599785460091,   0.002688051941039,   0.004431848411938,   0.002688051941039,   0.000599785460091,   0.000049233388666,
			0.000599785460091,   0.007306882745281,   0.032747176537767,   0.053990966513188,   0.032747176537767,   0.007306882745281,   0.000599785460091,
			0.002688051941039,   0.032747176537767,   0.146762663173740,   0.241970724519143,   0.146762663173740,   0.032747176537767,   0.002688051941039,
			0.004431848411938,   0.053990966513188,   0.241970724519143,   0.398942280401433,   0.241970724519143,   0.053990966513188,   0.004431848411938,
			0.002688051941039,   0.032747176537767,   0.146762663173740,   0.241970724519143,   0.146762663173740,   0.032747176537767,   0.002688051941039,
			0.000599785460091,   0.007306882745281,   0.032747176537767,   0.053990966513188,   0.032747176537767,   0.007306882745281,   0.000599785460091,
			0.000049233388666,   0.000599785460091,   0.002688051941039,   0.004431848411938,   0.002688051941039,   0.000599785460091,   0.000049233388666};

	MultidimArray<double> diffImage = I1 - I2;
    double 	imed = 0;
    int imiddle=YSIZE(I1)/2;
    int jmiddle=XSIZE(I1)/2;
    int R2max=imiddle*imiddle;
    int ysize=(int)YSIZE(I1);
    int xsize=(int)XSIZE(I1);
    for (int i=3; i<ysize-3; ++i)
    {
    	int i2=(i-imiddle)*(i-imiddle);
        for (int j=3; j<xsize-3; ++j)
        {
        	int j2=(j-jmiddle)*(j-jmiddle);
        	if (i2+j2>R2max) // Measure only within the maximum circle
        		continue;
        	double diffi=DIRECT_A2D_ELEM(diffImage,i,j);
        	int	index=0;
        	for (int ii=-3; ii<=3; ++ii)
        	{
        		refW = &w[index];
        		index = index + 7;
        	    double *diffj=&DIRECT_A2D_ELEM(diffImage,i+ii,j-3);
        	    double	imedAux=(*refW++)*(*diffj++);
            	imedAux+=(*refW++)*(*diffj++);
            	imedAux+=(*refW++)*(*diffj++);
            	imedAux+=(*refW++)*(*diffj++);
            	imedAux+=(*refW++)*(*diffj++);
            	imedAux+=(*refW++)*(*diffj++);
            	imedAux+=(*refW++)*(*diffj++);

            	imed+=imedAux*diffi;
        	}
        }
    }
	return sqrt(imed);
}

/* Imed distance ---------------------------------------------------------- */
double imedNormalizedDistance(const MultidimArray<double>& I1, const MultidimArray<double>& I2)
{
	// [x,y]=meshgrid([-3:1:3],[-3:1:3])
	// format long
	// w=1/sqrt(2*pi)*exp(-0.5*(x.*x+y.*y))
	const double w[7][7]={
			{0.000049233388666,   0.000599785460091,   0.002688051941039,   0.004431848411938,   0.002688051941039,   0.000599785460091,   0.000049233388666},
			{0.000599785460091,   0.007306882745281,   0.032747176537767,   0.053990966513188,   0.032747176537767,   0.007306882745281,   0.000599785460091},
			{0.002688051941039,   0.032747176537767,   0.146762663173740,   0.241970724519143,   0.146762663173740,   0.032747176537767,   0.002688051941039},
			{0.004431848411938,   0.053990966513188,   0.241970724519143,   0.398942280401433,   0.241970724519143,   0.053990966513188,   0.004431848411938},
			{0.002688051941039,   0.032747176537767,   0.146762663173740,   0.241970724519143,   0.146762663173740,   0.032747176537767,   0.002688051941039},
			{0.000599785460091,   0.007306882745281,   0.032747176537767,   0.053990966513188,   0.032747176537767,   0.007306882745281,   0.000599785460091},
			{0.000049233388666,   0.000599785460091,   0.002688051941039,   0.004431848411938,   0.002688051941039,   0.000599785460091,   0.000049233388666}
	};
    int imiddle=YSIZE(I1)/2;
    int jmiddle=XSIZE(I1)/2;
    int R2max=imiddle*imiddle;
    // Measure the average of both images
    double avg1=0;
    double avg2=0;
    double N=0;
    int ydim=YSIZE(I1);
    int xdim=XSIZE(I1);
    for (int i=3; i<ydim-3; ++i)
    {
    	int i2=(i-imiddle)*(i-imiddle);
        for (int j=3; j<xdim-3; ++j)
        {
        	int j2=(j-jmiddle)*(j-jmiddle);
        	if (i2+j2>R2max) // Measure only within the maximum circle
        		continue;
        	avg1+=DIRECT_A2D_ELEM(I1,i,j);
        	avg2+=DIRECT_A2D_ELEM(I2,i,j);
        	++N;
        }
    }
    avg1/=N;
    avg2/=N;
    double diffAvg=avg1-avg2;

    double imed = 0;
    double imed1=0;
    double imed2=0;
    for (int i=3; i<ydim-3; ++i)
    {
    	int i2=(i-imiddle)*(i-imiddle);
        for (int j=3; j<xdim-3; ++j)
        {
        	int j2=(j-jmiddle)*(j-jmiddle);
        	if (i2+j2>R2max) // Measure only within the maximum circle
        		continue;
        	double subtractedI1ij=DIRECT_A2D_ELEM(I1,i,j)-avg1;
        	double subtractedI2ij=DIRECT_A2D_ELEM(I2,i,j)-avg2;
        	double diffi=DIRECT_A2D_ELEM(I1,i,j)-DIRECT_A2D_ELEM(I2,i,j)-diffAvg;
        	for (int ii=-3; ii<=3; ++ii)
        	{
        		int i_ii=i+ii;
            	for (int jj=-3; jj<=3; ++jj)
            	{
            		int j_jj=j+jj;
                	//double diffj=DIRECT_A2D_ELEM(I1,i_ii,j_jj)-DIRECT_A2D_ELEM(I2,i_ii,j_jj)-diffAvg;
                	double diffj=DIRECT_A2D_ELEM(I1,i,j)-DIRECT_A2D_ELEM(I2,i_ii,j_jj)-diffAvg;
                	double wiijj=w[ii+3][jj+3];
            		imed+=wiijj*diffi*diffj;
            		imed1+=wiijj*subtractedI1ij*(DIRECT_A2D_ELEM(I1,i_ii,j_jj)-avg1);
            		imed2+=wiijj*subtractedI2ij*(DIRECT_A2D_ELEM(I2,i_ii,j_jj)-avg2);
            	}
        	}
        }
    }
	return imed/sqrt((imed1+imed2)/2);
}



/* Correlation measured inside a mask defined by the standard deviation value of the pixels of first image ---- */
double correlationMasked(const MultidimArray<double>& I1, const MultidimArray<double>& I2)
{
	double mean1;
    double std1;
    double th1;
	I1.computeAvgStdev(mean1,std1);
	th1 = std1;

	// Estimate the mean and stddev within the mask
	double N1=0;
	double sumMI1=0;
    double sumMI2=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I1)
	{
		double p1=DIRECT_MULTIDIM_ELEM(I1,n);
		double p2=DIRECT_MULTIDIM_ELEM(I2,n);
		if(p1>=th1){
			sumMI1+=p1;
			sumMI2+=p2;
			N1+=1.0;
		}
	}

	double sumMI1I2=0.0;
    double sumMI1I1=0.0;
    double sumMI2I2=0.0;
	double iN1=1;
    double avgM1=0;
    double avgM2=0;
    double corrM1M2;
	if (N1>0){
		iN1=1.0/N1;
		avgM1=sumMI1*iN1;
		avgM2=sumMI2*iN1;
	}
	else
		return 0;

	double p1a;
    double p2a;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I1)
	{
		double p1=DIRECT_MULTIDIM_ELEM(I1,n);
		double p2=DIRECT_MULTIDIM_ELEM(I2,n);
		if (p1>th1){
			p1a=p1-avgM1;
			p2a=p2-avgM2;
			sumMI1I1 +=p1a*p1a;
			sumMI2I2 +=p2a*p2a;
			sumMI1I2+=p1a*p2a;
		}
	}
	corrM1M2=sumMI1I2/sqrt(sumMI1I1*sumMI2I2);

	return corrM1M2;
}



/* Weighted correlation measured with a weight defined by the difference between both images ---- */
double correlationWeighted(MultidimArray<double>& I1, MultidimArray<double>& I2)
{

	MultidimArray<double> Idiff;
    MultidimArray<double> I2Aligned;
	Matrix2D<double> M;

	I1.setXmippOrigin();
	I2.setXmippOrigin();
	I2Aligned=I2;
	alignImages(I1, I2Aligned, M, false);
	Idiff=I1;
	Idiff-=I2Aligned;

	double mean;
    double std;
	Idiff.computeAvgStdev(mean,std);
	Idiff.selfABS();
	double threshold=std;

	// Estimate the mean and stddev within the mask
	double N=0;
	double sumWI1=0;
    double sumWI2=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Idiff)
	{
		if (DIRECT_MULTIDIM_ELEM(Idiff,n)>threshold)
		{
			double p1=DIRECT_MULTIDIM_ELEM(I1,n);
			double p2Alg=DIRECT_MULTIDIM_ELEM(I2Aligned,n);
			sumWI1+=p1;
			sumWI2+=p2Alg;
			N+=1.0;
		}
	}

	double sumWI1I2=0.0;
    double sumWI1I1=0.0;
    double sumWI2I2=0.0;
	double iN;
    double avgW1;
    double avgW2;
    double corrW1W2;
	if(N>0){
		iN=1.0/N;
		avgW1=sumWI1*iN;
		avgW2=sumWI2*iN;
	}
	double p1a;
    double p2a;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I1)
	{
		if (DIRECT_MULTIDIM_ELEM(Idiff,n)>threshold)
		{
			double p1=DIRECT_MULTIDIM_ELEM(I1,n);
			double p2Alg=DIRECT_MULTIDIM_ELEM(I2Aligned,n);
			p1a=p1-avgW1;
			p2a=p2Alg-avgW2;

			double w=DIRECT_MULTIDIM_ELEM(Idiff,n);
			double wp1a=w*p1a;
			double wp2a=w*p2a;

			sumWI1I1 +=wp1a*p1a;
			sumWI2I2 +=wp2a*p2a;
			sumWI1I2 +=wp1a*p2a;
		}
	}
	corrW1W2=sumWI1I2/sqrt(sumWI1I1*sumWI2I2);

	return corrW1W2;
}




double svdCorrelation(const MultidimArray<double>& I1, const MultidimArray<double>& I2, const MultidimArray< int >* mask)
{
	// Copy input images into matrix2Ds
	Matrix2D<double> mI1;
    Matrix2D<double> mI2;
	I1.copy(mI1);
	I2.copy(mI2);

	// SVD Decomposition of I1
	Matrix2D<double> U1;
    Matrix2D<double> V1;
	Matrix1D<double> W1;
	svdcmp(mI1,U1,W1,V1);

	// Express I2 in the basis of I1
	Matrix2D<double> aux;
    Matrix2D<double> W2;
	matrixOperation_AB(mI2,V1,aux);
	matrixOperation_AtB(U1,aux,W2);

	// Make W2 to be diagonal
	FOR_ALL_ELEMENTS_IN_MATRIX2D(W2)
	if (i!=j)
		MAT_ELEM(W2,i,j)=0.0;

	// Recover I2p
	matrixOperation_ABt(W2,V1,aux);
	matrixOperation_AB(U1,aux,mI2);

	// Measure correlation
	MultidimArray<double> I2p;
	I2p=mI2;
//	Image<double> save;
//	save()=I1;
//	save.write("PPPI1.xmp");
//	save()=I2;
//	save.write("PPPI2.xmp");
//	save()=I2p;
//	save.write("PPPI2p.xmp");
//	std::cout << "Correlation= " << correlationIndex(I1,I2p,mask) << std::endl;
//	std::cout << "Press any key\n";
//	char c; std::cin >> c;
	// double ratio=I2p.computeStddev()/I2.computeStddev();
	return correlationIndex(I1,I2p,mask); //*ratio;
}

/* Covariance matrix ------------------------------------------------------ */
void covarianceMatrix(const MultidimArray<double> &I, Matrix2D<double> &C)
{
	Matrix2D<double> mI;
	I.copy(mI);
	subtractColumnMeans(mI);
	matrixOperation_AtA(mI,C);
	C*=1.0/(YSIZE(I)-1.0);
}

/* Best shift -------------------------------------------------------------- */
template float bestShift(MultidimArray<float>&, float&, float&, const MultidimArray<int>*, int);
template<typename T>
T bestShift(MultidimArray<T> &Mcorr,
               T &shiftX, T &shiftY, const MultidimArray<int> *mask, int maxShift)
{
    int imax = INT_MIN;
    int jmax;
    int i_actual;
    int j_actual;
    double xmax;
    double ymax;
    double avecorr;
    double stdcorr;
    double dummy;
    bool neighbourhood = true;

    /*
     Warning: for masks with a small number of non-zero pixels, this routine is NOT reliable...
     Anyway, maybe using a mask is not a good idea at al...
     */

    // Adjust statistics within shiftmask to average 0 and stddev 1
    if (mask != nullptr)
    {
        if ((*mask).sum() < 2)
        {
            shiftX = shiftY = 0.;
            return -1;
        }
        else
        {
            computeStats_within_binary_mask(*mask, Mcorr, dummy, dummy, avecorr,
                                            stdcorr);
            double istdcorr = 1.0 / stdcorr;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mcorr)
            if (DIRECT_MULTIDIM_ELEM(*mask, n))
                DIRECT_MULTIDIM_ELEM(Mcorr, n) =
                    (DIRECT_MULTIDIM_ELEM(Mcorr, n) - avecorr)
                    * istdcorr;
            else
                DIRECT_MULTIDIM_ELEM(Mcorr, n) = 0.;
        }
    }
    else
        Mcorr.statisticsAdjust((T)0, (T)1);

    // Look for maximum shift
    if (maxShift==-1)
    	Mcorr.maxIndex(imax, jmax);
    else
    {
    	int maxShift2=maxShift*maxShift;
    	double bestCorr=std::numeric_limits<T>::lowest();
    	for (int i=-maxShift; i<=maxShift; i++)
    		for (int j=-maxShift; j<=maxShift; j++)
    		{
    			if (i*i+j*j>maxShift2) // continue if the Euclidean distance is too far
    				continue;
    			else if (A2D_ELEM(Mcorr, i, j)>bestCorr)
    			{
    				imax=i;
    				jmax=j;
    				bestCorr=A2D_ELEM(Mcorr, imax, jmax);
    			}
    		}
    }
    double max = A2D_ELEM(Mcorr, imax, jmax);

    // Estimate n_max around the maximum
    int n_max = -1;
    while (neighbourhood)
    {
        n_max++;
        for (int i = -n_max; i <= n_max && neighbourhood; i++)
        {
            i_actual = i + imax;
            if (i_actual < STARTINGY(Mcorr) || i_actual > FINISHINGY(Mcorr))
            {
                neighbourhood = false;
                break;
            }
            for (int j = -n_max; j <= n_max && neighbourhood; j++)
            {
                j_actual = j + jmax;
                if (j_actual < STARTINGX(Mcorr) || j_actual > FINISHINGX(Mcorr))
                {
                    neighbourhood = false;
                    break;
                }
                else if (max / 1.414 > A2D_ELEM(Mcorr, i_actual, j_actual))
                {
                    neighbourhood = false;
                    break;
                }
            }
        }
    }

    // We have the neighbourhood => looking for the gravity centre
    xmax = ymax = 0.;
    double sumcorr = 0.;
    if (imax-n_max<STARTINGY(Mcorr))
        n_max=std::min(imax-STARTINGY(Mcorr),n_max);
    if (imax+n_max>FINISHINGY(Mcorr))
        n_max=std::min(FINISHINGY(Mcorr)-imax,n_max);
    if (jmax-n_max<STARTINGY(Mcorr))
        n_max=std::min(jmax-STARTINGX(Mcorr),n_max);
    if (jmax+n_max>FINISHINGY(Mcorr))
        n_max=std::min(FINISHINGX(Mcorr)-jmax,n_max);
    for (int i = -n_max; i <= n_max; i++)
    {
        i_actual = i + imax;
        for (int j = -n_max; j <= n_max; j++)
        {
            j_actual = j + jmax;
            double val = A2D_ELEM(Mcorr, i_actual, j_actual);
            ymax += i_actual * val;
            xmax += j_actual * val;
            sumcorr += val;
        }
    }
    if (sumcorr != 0)
    {
        shiftX = xmax / sumcorr;
        shiftY = ymax / sumcorr;
    }
    return max;
}

double bestShift(const MultidimArray<double> &I1, const MultidimArray< std::complex<double> > &FFTI1,
		       const MultidimArray<double> &I2,
               double &shiftX, double &shiftY, CorrelationAux &aux,
               const MultidimArray<int> *mask, int maxShift)
{
    I1.checkDimension(2);
    I2.checkDimension(2);

    MultidimArray<double> Mcorr;

    correlation_matrix(FFTI1, I2, Mcorr, aux);
    return bestShift(Mcorr, shiftX, shiftY, mask, maxShift);
}

double bestShift(const MultidimArray<double> &I1, const MultidimArray<double> &I2,
               double &shiftX, double &shiftY, CorrelationAux &aux,
               const MultidimArray<int> *mask, int maxShift)
{
    aux.transformer1.FourierTransform((MultidimArray<double> &)I1, aux.FFT1, false);
    return bestShift(I1,aux.FFT1,I2,shiftX,shiftY,aux,mask,maxShift);
}


double bestShift(const MultidimArray< std::complex<double> > &FFTI1,
				const MultidimArray< std::complex<double> > &FFTI2,
				MultidimArray<double> &Mcorr,
               double &shiftX, double &shiftY, CorrelationAux &aux,
               const MultidimArray<int> *mask, int maxShift)
{
	correlation_matrix(FFTI1, FFTI2, Mcorr, aux);
	return bestShift(Mcorr, shiftX, shiftY, mask, maxShift);
}

/* Best shift -------------------------------------------------------------- */
void bestShift(const MultidimArray<double> &I1, const MultidimArray<double> &I2,
               double &shiftX, double &shiftY, double &shiftZ, CorrelationAux &aux,
               const MultidimArray<int> *mask)
{
    I1.checkDimension(3);
    I2.checkDimension(3);

    int imax;
    int jmax;
    int kmax;
    int i_actual;
    int j_actual;
    int k_actual;
    double max;
    double xmax;
    double ymax;
    double zmax;
    double sumcorr;
    double avecorr;
    double stdcorr;
    double dummy;
    bool neighbourhood = true;
    MultidimArray<double> Mcorr;

    correlation_matrix(I1, I2, Mcorr, aux);

    /*
     Warning: for masks with a small number of non-zero pixels, this routine is NOT reliable...
     Anyway, maybe using a mask is not a good idea at al...
     */

    // Adjust statistics within shiftmask to average 0 and stddev 1
    if (mask != nullptr)
    {
        if ((*mask).sum() < 2)
        {
            shiftX = shiftY = shiftZ = 0.;
            return;
        }
        else
        {
            computeStats_within_binary_mask(*mask, Mcorr, dummy, dummy, avecorr,
                                            stdcorr);
            double istdcorr = 1.0 / stdcorr;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mcorr)
            if (DIRECT_MULTIDIM_ELEM(*mask, n))
                DIRECT_MULTIDIM_ELEM(Mcorr, n) =
                    (DIRECT_MULTIDIM_ELEM(Mcorr, n) - avecorr)
                    * istdcorr;
            else
                DIRECT_MULTIDIM_ELEM(Mcorr, n) = 0.;
        }
    }
    else
        Mcorr.statisticsAdjust(0., 1.);
    Mcorr.maxIndex(kmax, imax, jmax);
    max = A3D_ELEM(Mcorr, kmax, imax, jmax);

    // Estimate n_max around the maximum
    int n_max = -1;
    while (neighbourhood)
    {
        n_max++;
        for (int k = -n_max; k <= n_max && neighbourhood; k++)
        {
            k_actual = k + kmax;
            if (k_actual < STARTINGZ(Mcorr) || k_actual > FINISHINGZ(Mcorr))
            {
                neighbourhood = false;
                break;
            }
            for (int i = -n_max; i <= n_max && neighbourhood; i++)
            {
                i_actual = i + imax;
                if (i_actual < STARTINGY(Mcorr) || i_actual > FINISHINGY(Mcorr))
                {
                    neighbourhood = false;
                    break;
                }
                for (int j = -n_max; j <= n_max && neighbourhood; j++)
                {
                    j_actual = j + jmax;
                    if (j_actual < STARTINGX(Mcorr)
                        || j_actual > FINISHINGX(Mcorr))
                    {
                        neighbourhood = false;
                        break;
                    }
                    else if (max
                             / 1.414 > A3D_ELEM(Mcorr, k_actual, i_actual, j_actual))
                    {
                        neighbourhood = false;
                        break;
                    }
                }
            }
        }
    }

    // We have the neighbourhood => looking for the gravity centre
    zmax = xmax = ymax = sumcorr = 0.;
    if (kmax-n_max<STARTINGZ(Mcorr))
        n_max=std::min(kmax-STARTINGZ(Mcorr),n_max);
    if (kmax+n_max>FINISHINGZ(Mcorr))
        n_max=std::min(FINISHINGZ(Mcorr)-kmax,n_max);
    if (imax-n_max<STARTINGY(Mcorr))
        n_max=std::min(imax-STARTINGY(Mcorr),n_max);
    if (imax+n_max>FINISHINGY(Mcorr))
        n_max=std::min(FINISHINGY(Mcorr)-imax,n_max);
    if (jmax-n_max<STARTINGY(Mcorr))
        n_max=std::min(jmax-STARTINGX(Mcorr),n_max);
    if (jmax+n_max>FINISHINGY(Mcorr))
        n_max=std::min(FINISHINGX(Mcorr)-jmax,n_max);
    for (int k = -n_max; k <= n_max; k++)
    {
        k_actual = k + kmax;
        for (int i = -n_max; i <= n_max; i++)
        {
            i_actual = i + imax;
            for (int j = -n_max; j <= n_max; j++)
            {
                j_actual = j + jmax;
                double val = A3D_ELEM(Mcorr, k_actual, i_actual, j_actual);
                zmax += k_actual * val;
                ymax += i_actual * val;
                xmax += j_actual * val;
                sumcorr += val;
            }
        }
    }
    if (sumcorr != 0)
    {
        shiftX = xmax / sumcorr;
        shiftY = ymax / sumcorr;
        shiftZ = zmax / sumcorr;
    }
}

/* Best non-wrapping shift ------------------------------------------------- */
//#define DEBUG
void bestNonwrappingShift(const MultidimArray<double> &I1,
                          const MultidimArray<double> &I2, double &shiftX, double &shiftY,
                          CorrelationAux &aux)
{
    aux.transformer1.FourierTransform((MultidimArray<double> &)I1, aux.FFT1, false);
    bestNonwrappingShift(I1,aux.FFT1,I2,shiftX,shiftY,aux);
}

void bestNonwrappingShift(const MultidimArray<double> &I1, const MultidimArray< std::complex<double> >&FFTI1,
                          const MultidimArray<double> &I2, double &shiftX, double &shiftY,
                          CorrelationAux &aux)
{
    I1.checkDimension(2);
    I2.checkDimension(2);

    bestShift(I1, FFTI1, I2, shiftX, shiftY, aux);
    double bestCorr;
    double corr;
    MultidimArray<double> Iaux;

    translate(1, Iaux, I1, vectorR2(-shiftX, -shiftY), xmipp_transformation::DONT_WRAP);
    //I1.translate(vectorR2(-shiftX,-shiftY),Iaux, DONT_WRAP);
    bestCorr = corr = fastCorrelation(I2, Iaux);
    double finalX = shiftX;
    double finalY = shiftY;
#ifdef DEBUG

    std::cout << "shiftX=" << shiftX << " shiftY=" << shiftY
    << " corr=" << corr << std::endl;
    ImageXmipp save;
    save()=I1;
    save.write("PPPI1.xmp");
    save()=I2;
    save.write("PPPI2.xmp");
    save()=Iaux;
    save.write("PPPpp.xmp");
#endif

    Iaux.initZeros();
    double testX = (shiftX > 0) ? (shiftX - XSIZE(I1)) : (shiftX + XSIZE(I1));
    double testY = shiftY;
    translate(1, Iaux, I1, vectorR2(-testX, -testY), xmipp_transformation::DONT_WRAP);
    //I1.translate(vectorR2(-testX,-testY),Iaux,DONT_WRAP);
    corr = fastCorrelation(I2, Iaux);
    if (corr > bestCorr)
        finalX = testX;
#ifdef DEBUG

    std::cout << "shiftX=" << testX << " shiftY=" << testY
    << " corr=" << corr << std::endl;
    save()=Iaux;
    save.write("PPPmp.xmp");
#endif

    Iaux.initZeros();
    testX = shiftX;
    testY = (shiftY > 0) ? (shiftY - YSIZE(I1)) : (shiftY + YSIZE(I1));
    translate(1, Iaux, I1, vectorR2(-testX, -testY), xmipp_transformation::DONT_WRAP);
    //I1.translate(vectorR2(-testX,-testY),Iaux,DONT_WRAP);
    corr = fastCorrelation(I2, Iaux);
    if (corr > bestCorr)
        finalY = testY;
#ifdef DEBUG

    std::cout << "shiftX=" << testX << " shiftY=" << testY
    << " corr=" << corr << std::endl;
    save()=Iaux;
    save.write("PPPpm.xmp");
#endif

    Iaux.initZeros();
    testX = (shiftX > 0) ? (shiftX - XSIZE(I1)) : (shiftX + XSIZE(I1));
    testY = (shiftY > 0) ? (shiftY - YSIZE(I1)) : (shiftY + YSIZE(I1));
    translate(1, Iaux, I1, vectorR2(-testX, -testY), xmipp_transformation::DONT_WRAP);
    //I1.translate(vectorR2(-testX,-testY),Iaux,DONT_WRAP);
    corr = fastCorrelation(I2, Iaux);
    if (corr > bestCorr)
    {
        finalX = testX;
        finalY = testY;
    }
#ifdef DEBUG
    std::cout << "shiftX=" << testX << " shiftY=" << testY
    << " corr=" << corr << std::endl;
    save()=Iaux;
    save.write("PPPmm.xmp");
#endif

    shiftX = finalX;
    shiftY = finalY;
}
#undef DEBUG

/* Best shift -------------------------------------------------------------- */
double bestShiftRealSpace(const MultidimArray<double> &I1, MultidimArray<double> &I2,
               double &shiftX, double &shiftY,
               const MultidimArray<int> *mask, int maxShift, double shiftStep)
{
    I1.checkDimension(2);
    I2.checkDimension(2);

    double bestCorr=-1e38;

    MultidimArray<double> alignedI2;
    MultidimArray<double> bestI2;
    int maxShift2=maxShift*maxShift;
    Matrix1D<double> shift(2);
    for (double y=-maxShift; y<=maxShift; y+=shiftStep)
    	for (double x=-maxShift; x<=maxShift; x+=shiftStep)
    	{
    		if (y*y+x*x>maxShift2)
    			continue;
    		YY(shift)=y;
    		XX(shift)=x;
    		translate(xmipp_transformation::LINEAR,alignedI2,I2,shift,xmipp_transformation::DONT_WRAP,0.0);

    		double corr=correlationIndex(I1,alignedI2,mask);
    		if (corr>bestCorr)
    		{
    			bestI2=alignedI2;
    			bestCorr=corr;
    			shiftY=y;
    			shiftX=x;
    		}
    	}
    I2=bestI2;

    return bestCorr;
}

/* Align two images -------------------------------------------------------- */
AlignmentAux::AlignmentAux()
{
    plans = nullptr;
}
AlignmentAux::~AlignmentAux()
{
    delete plans;
}

void computeAlignmentTransforms(const MultidimArray<double>& I, AlignmentTransforms &ITransforms,
		AlignmentAux &aux, CorrelationAux &aux2)
{
	aux2.transformer1.FourierTransform((MultidimArray<double> &)I, ITransforms.FFTI, false);
    polarFourierTransform<true>(I, ITransforms.polarFourierI, false, XSIZE(I) / 5, XSIZE(I) / 2, aux.plans, 1);
}

constexpr double SHIFT_THRESHOLD = 0.95;             // Shift threshold in pixels.
constexpr float  ROTATE_THRESHOLD  =	1.0	;		// Rotate threshold in degrees.
constexpr double  INITIAL_SHIFT_THRESHOLD =	SHIFT_THRESHOLD + 1.0;		// Shift threshold in pixels.
constexpr float  INITIAL_ROTATE_THRESHOLD = ROTATE_THRESHOLD + 1.0	;	// Rotate threshold in degrees.

double alignImages(const MultidimArray<double>& Iref, const AlignmentTransforms& IrefTransforms, MultidimArray<double>& I,
                   Matrix2D<double>&M, bool wrap, AlignmentAux &aux, CorrelationAux &aux2,
                   RotationalCorrelationAux &aux3)
{
    I.checkDimension(2);

    aux.ARS.initIdentity(3);
    aux.ASR.initIdentity(3);
    aux.IauxSR = I;
    aux.IauxRS = I;

    aux.rotationalCorr.resize(2 * IrefTransforms.polarFourierI.getSampleNoOuterRing() - 1);
    aux3.local_transformer.setReal(aux.rotationalCorr);

	double shiftXSR=INITIAL_SHIFT_THRESHOLD;
    double shiftYSR=INITIAL_SHIFT_THRESHOLD;
    double bestRotSR=INITIAL_ROTATE_THRESHOLD;
	double shiftXRS=INITIAL_SHIFT_THRESHOLD;
    double shiftYRS=INITIAL_SHIFT_THRESHOLD;
    double bestRotRS=INITIAL_ROTATE_THRESHOLD;

    // Align the image with the reference
    for (int i = 0; i < 3; i++)
    {
        // Shift then rotate
		if (((shiftXSR > SHIFT_THRESHOLD) || (shiftXSR < (-SHIFT_THRESHOLD))) ||
			((shiftYSR > SHIFT_THRESHOLD) || (shiftYSR < (-SHIFT_THRESHOLD))))
		{
			bestNonwrappingShift(Iref, IrefTransforms.FFTI, aux.IauxSR, shiftXSR, shiftYSR, aux2);
			MAT_ELEM(aux.ASR,0,2) += shiftXSR;
			MAT_ELEM(aux.ASR,1,2) += shiftYSR;
			applyGeometry(xmipp_transformation::LINEAR, aux.IauxSR, I, aux.ASR, xmipp_transformation::IS_NOT_INV, wrap);
		}

		if (bestRotSR > ROTATE_THRESHOLD)
		{
			polarFourierTransform<true>(aux.IauxSR, aux.polarFourierI, true,
											XSIZE(Iref) / 5, XSIZE(Iref) / 2, aux.plans, 1);

			bestRotSR = best_rotation(IrefTransforms.polarFourierI, aux.polarFourierI, aux3);
			rotation2DMatrix(bestRotSR, aux.R);
			aux.ASR = aux.R * aux.ASR;
			applyGeometry(xmipp_transformation::LINEAR, aux.IauxSR, I, aux.ASR, xmipp_transformation::IS_NOT_INV, wrap);
		}

        // Rotate then shift
		if (bestRotRS > ROTATE_THRESHOLD)
		{
			polarFourierTransform<true>(aux.IauxRS, aux.polarFourierI, true,
											XSIZE(Iref) / 5, XSIZE(Iref) / 2, aux.plans, 1);
			bestRotRS = best_rotation(IrefTransforms.polarFourierI, aux.polarFourierI, aux3);
			rotation2DMatrix(bestRotRS, aux.R);
			aux.ARS = aux.R * aux.ARS;
			applyGeometry(xmipp_transformation::LINEAR, aux.IauxRS, I, aux.ARS, xmipp_transformation::IS_NOT_INV, wrap);
		}

		if (((shiftXRS > SHIFT_THRESHOLD) || (shiftXRS < (-SHIFT_THRESHOLD))) ||
			((shiftYRS > SHIFT_THRESHOLD) || (shiftYRS < (-SHIFT_THRESHOLD))))
		{
			bestNonwrappingShift(Iref, IrefTransforms.FFTI, aux.IauxRS, shiftXRS, shiftYRS, aux2);
			MAT_ELEM(aux.ARS,0,2) += shiftXRS;
			MAT_ELEM(aux.ARS,1,2) += shiftYRS;
			applyGeometry(xmipp_transformation::LINEAR, aux.IauxRS, I, aux.ARS, xmipp_transformation::IS_NOT_INV, wrap);
		}
    }

    double corrRS = correlationIndex(aux.IauxRS, Iref);
    double corrSR = correlationIndex(aux.IauxSR, Iref);
    double corr;
    if (corrRS > corrSR)
    {
        I = aux.IauxRS;
        M = aux.ARS;
        corr = corrRS;
    }
    else
    {
        I = aux.IauxSR;
        M = aux.ASR;
        corr = corrSR;
    }
    return corr;
}

double alignImages(const MultidimArray<double>& Iref, MultidimArray<double>& I,
                   Matrix2D<double>&M, bool wrap, AlignmentAux &aux, CorrelationAux &aux2,
                   RotationalCorrelationAux &aux3)
{
    Iref.checkDimension(2);
    AlignmentTransforms IrefTransforms;
    computeAlignmentTransforms(Iref,IrefTransforms, aux, aux2);
    return alignImages(Iref, IrefTransforms, I, M, wrap, aux, aux2, aux3);
}

double alignImages(const MultidimArray<double>& Iref, MultidimArray<double>& I,
                   Matrix2D<double>&M, bool wrap)
{
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    return alignImages(Iref, I, M, wrap, aux, aux2, aux3);
}

double alignImagesConsideringMirrors(const MultidimArray<double>& Iref, const AlignmentTransforms& IrefTransforms,
                                     MultidimArray<double>& I, Matrix2D<double> &M, AlignmentAux& aux,
                                     CorrelationAux& aux2, RotationalCorrelationAux &aux3, bool wrap,
                                     const MultidimArray<int>* mask)
{
    MultidimArray<double> Imirror;
    Matrix2D<double> Mmirror;
    Imirror = I;
    Imirror.selfReverseX();
    Imirror.setXmippOrigin();

    double corr=alignImages(Iref, IrefTransforms, I, M, wrap, aux, aux2, aux3);
    double corrMirror=alignImages(Iref, IrefTransforms, Imirror, Mmirror, wrap, aux, aux2, aux3);
    if (mask!=nullptr)
    {
    	corr = correlationIndex(Iref, I, mask);
    	corrMirror = correlationIndex(Iref, Imirror, mask);
    }
    double bestCorr = corr;
    if (corrMirror > bestCorr)
    {
        bestCorr = corrMirror;
        I = Imirror;
        M = Mmirror;
        MAT_ELEM(M,0,0) *= -1;
        MAT_ELEM(M,1,0) *= -1;
    }
    return bestCorr;
}

double alignImagesConsideringMirrors(const MultidimArray<double>& Iref, MultidimArray<double>& I,
                   Matrix2D<double>&M, bool wrap)
{
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    return alignImagesConsideringMirrors(Iref, I, M, aux, aux2, aux3, wrap, nullptr);
}

double alignImagesConsideringMirrors(const MultidimArray<double>& Iref,
                                     MultidimArray<double>& I, Matrix2D<double> &M, AlignmentAux& aux,
                                     CorrelationAux& aux2, RotationalCorrelationAux &aux3, bool wrap,
                                     const MultidimArray<int>* mask)
{
    AlignmentTransforms IrefTransforms;
    computeAlignmentTransforms(Iref,IrefTransforms, aux, aux2);
    return alignImagesConsideringMirrors(Iref, IrefTransforms, I, M, aux, aux2, aux3, wrap, mask);
}

void alignSetOfImages(MetaData &MD, MultidimArray<double>& Iavg, int Niter,
                      bool considerMirror)
{
    Image<double> I;
    MultidimArray<double> InewAvg;
    FileName fnImg;
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    Matrix2D<double> M;
    size_t Nimgs;
    size_t Xdim;
    size_t Ydim;
    size_t Zdim;
    getImageSize(MD, Xdim, Ydim, Zdim, Nimgs);
    for (int n = 0; n < Niter; ++n)
    {
        bool lastIteration = (n == (Niter - 1));
        InewAvg.initZeros(Ydim, Xdim);
        InewAvg.setXmippOrigin();
        for (size_t objId : MD.ids())
        {
            MD.getValue(MDL_IMAGE, fnImg, objId);
            I.read(fnImg);
            I().setXmippOrigin();
            double corr;
            if (considerMirror)
                corr = alignImagesConsideringMirrors(Iavg, I(), M, aux, aux2,
                                                     aux3, xmipp_transformation::WRAP);
            else
                corr = alignImages(Iavg, I(), M, xmipp_transformation::WRAP, aux, aux2, aux3);
            InewAvg += I();
            if (n == 0)
                Iavg = InewAvg;
            if (lastIteration)
            {
                double scale;
                double shiftx;
                double shifty;
                double psi;
                bool flip;
                transformationMatrix2Parameters2D(M, flip, scale, shiftx,
                                                  shifty, psi);
                MD.setValue(MDL_FLIP, flip, objId);
                MD.setValue(MDL_SHIFT_X, shiftx, objId);
                MD.setValue(MDL_SHIFT_Y, shifty, objId);
                MD.setValue(MDL_ANGLE_PSI, psi, objId);
                MD.setValue(MDL_MAXCC, corr, objId);
            }
        }
        InewAvg /= Nimgs;
        Iavg = InewAvg;
        centerImage(Iavg, aux2, aux3, 4);
    }
}

double fastBestRotation(const MultidimArray<double>& IrefCyl,
                        const MultidimArray<double>& Icyl, CorrelationAux &aux,
                        VolumeAlignmentAux &aux2, double deltaAng)
{
    correlation_matrix(IrefCyl, Icyl, aux2.corr, aux, false);
    STARTINGZ(aux2.corr) = STARTINGY(aux2.corr) = STARTINGX(aux2.corr) = 0;
    double bestCorr = A3D_ELEM(aux2.corr,0,STARTINGY(aux2.corr),0);
    double bestAngle = 0;
    for (int i = STARTINGY(aux2.corr) + 1; i <= FINISHINGY(aux2.corr); i++)
    {
        double corr = A3D_ELEM(aux2.corr,0,i,0);
        if (corr > bestCorr)
        {
            bestCorr = corr;
            bestAngle = i;
        }
    }
    bestAngle *= deltaAng * 180.0 / PI;
    return -bestAngle;
}

double fastBestRotationAroundZ(const MultidimArray<double>& IrefCyl,
                               const MultidimArray<double>& I, CorrelationAux &aux,
                               VolumeAlignmentAux &aux2)
{
    double deltaAng = atan(2.0 / XSIZE(I));
    Matrix1D<double> v(3);
    XX(v) = 0;
    YY(v) = 0;
    ZZ(v) = 1;
    volume_convertCartesianToCylindrical(I, aux2.Icyl, 3, XSIZE(I) / 2, 1, 0,
                                         2 * PI, deltaAng, v);
    return fastBestRotation(IrefCyl,aux2.Icyl,aux,aux2,deltaAng);
}

double fastBestRotationAroundY(const MultidimArray<double>& IrefCyl,
                               const MultidimArray<double>& I, CorrelationAux &aux,
                               VolumeAlignmentAux &aux2)
{
    double deltaAng = atan(2.0 / XSIZE(I));
    Matrix1D<double> v(3);
    XX(v) = 0;
    YY(v) = 1;
    ZZ(v) = 0;
    volume_convertCartesianToCylindrical(I, aux2.Icyl, 3, XSIZE(I) / 2, 1, 0,
                                         2 * PI, deltaAng, v);
    return -fastBestRotation(IrefCyl,aux2.Icyl,aux,aux2,deltaAng);
}

double fastBestRotationAroundX(const MultidimArray<double>& IrefCyl,
                               const MultidimArray<double>& I, CorrelationAux &aux,
                               VolumeAlignmentAux &aux2)
{
    double deltaAng = atan(2.0 / XSIZE(I));
    Matrix1D<double> v(3);
    XX(v) = 1;
    YY(v) = 0;
    ZZ(v) = 0;
    volume_convertCartesianToCylindrical(I, aux2.Icyl, 3, XSIZE(I) / 2, 1, 0,
                                         2 * PI, deltaAng, v);
    return fastBestRotation(IrefCyl,aux2.Icyl,aux,aux2,deltaAng);
}

void fastBestRotation(const MultidimArray<double>& IrefCylZ,
                      const MultidimArray<double>& IrefCylY,
                      const MultidimArray<double>& IrefCylX,
                      const MultidimArray<double>& I,
                      const MultidimArray<double>& Icurrent,
                      MultidimArray<double>& Ifinal,
                      char axis,
                      Matrix2D<double> &R,
                      CorrelationAux &aux, VolumeAlignmentAux &aux2)
{
    double bestAngle;
    if (axis=='Z')
        bestAngle = fastBestRotationAroundZ(IrefCylZ, Icurrent, aux, aux2);
    else if (axis=='Y')
        bestAngle = fastBestRotationAroundY(IrefCylY, Icurrent, aux, aux2);
    else
        bestAngle = fastBestRotationAroundX(IrefCylX, Icurrent, aux, aux2);
    Matrix2D<double> Raux;
    rotation3DMatrix(bestAngle, axis, Raux);
    R=Raux*R;
    applyGeometry(xmipp_transformation::LINEAR, Ifinal, I, R, xmipp_transformation::IS_NOT_INV, xmipp_transformation::WRAP);
}

void fastBestRotation(const MultidimArray<double>& IrefCylZ,
                      const MultidimArray<double>& IrefCylY,
                      const MultidimArray<double>& IrefCylX,
                      MultidimArray<double>& I,
                      const String &eulerAngles,
                      Matrix2D<double> &R,
                      CorrelationAux &aux, VolumeAlignmentAux &aux2)
{
    R.initIdentity(4);
    fastBestRotation(IrefCylZ,IrefCylY,IrefCylX,I,I,aux2.I1,eulerAngles[0],R, aux,aux2);
    fastBestRotation(IrefCylZ,IrefCylY,IrefCylX,I,aux2.I1,aux2.I12,eulerAngles[1],R,aux,aux2);
    fastBestRotation(IrefCylZ,IrefCylY,IrefCylX,I,aux2.I12,aux2.I123,eulerAngles[2],R,aux,aux2);
    I=aux2.I123;
}

double bestRotationAroundZ(const MultidimArray<double>& Iref,
                           const MultidimArray<double>& I, CorrelationAux &aux,
                           VolumeAlignmentAux &aux2)
{
    Iref.checkDimension(3);
    I.checkDimension(3);
    if (!I.sameShape(Iref))
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     "Both volumes should be of the same shape");

    double deltaAng = atan(2.0 / XSIZE(I));
    Matrix1D<double> v(3);
    XX(v) = 0;
    YY(v) = 0;
    ZZ(v) = 1;
    volume_convertCartesianToCylindrical(Iref, aux2.IrefCyl, 3, XSIZE(I) / 2, 1,
                                         0, 2 * PI, deltaAng, v);
    return fastBestRotationAroundZ(aux2.IrefCyl, I, aux, aux2);
}

/* Estimate 2D Gaussian ---------------------------------------------------- */
/* See Brandle, Chen, Bischof, Lapp. Robust parametric and semi-parametric
 spot fitting for spot array images. 2000 */
double unnormalizedGaussian2D(const Matrix1D<double> &r,
                              const Matrix1D<double> &mu, const Matrix2D<double> &sigmainv)
{
    double x = XX(r) - XX(mu);
    double y = YY(r) - YY(mu);
    return exp(
               -0.5
               * (sigmainv(0, 0) * x * x + 2 * sigmainv(0, 1) * x * y
                  + sigmainv(1, 1) * y * y));
}

void estimateGaussian2D(const MultidimArray<double> &I, double &a, double &b,
                        Matrix1D<double> &mu, Matrix2D<double> &sigma, bool estimateMu,
                        int iterations)
{
    I.checkDimension(2);

    MultidimArray<double> z(I);

    // Estimate b
    Histogram1D hist;
    compute_hist(z, hist, 100);
    b = hist.percentil(5);

    // Iteratively estimate all parameters
    for (int n = 0; n < iterations; n++)
    {
        // Reestimate z
        FOR_ALL_ELEMENTS_IN_ARRAY2D(z)
        z(i, j) = XMIPP_MAX(I(i,j)-b,0);

        // Sum of z
        double T = z.sum();

        // Estimate center
        mu.initZeros(2);
        if (estimateMu)
        {
            FOR_ALL_ELEMENTS_IN_ARRAY2D(z)
            {
                double val = z(i, j);
                XX(mu) += val * j;
                YY(mu) += val * i;
            }
            mu /= T;
        }

        // Estimate sigma
        sigma.initZeros(2, 2);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(z)
        {
            double val = z(i, j);
            double j_mu = j - XX(mu);
            double i_mu = i - YY(mu);
            sigma(0, 0) += val * j_mu * j_mu;
            sigma(0, 1) += val * i_mu * j_mu;
            sigma(1, 1) += val * i_mu * i_mu;
        }
        sigma(1, 0) = sigma(0, 1);
        sigma /= T;

        // Estimate amplitude
        Matrix2D<double> sigmainv = sigma.inv();
        Matrix1D<double> r(2);
        double G2 = 0;
        a = 0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(z)
        {
            XX(r) = j;
            YY(r) = i;
            double G = unnormalizedGaussian2D(r, mu, sigmainv);
            a += z(i, j) * G;
            G2 += G * G;
        }
        a /= G2;

        // Reestimate b
        FOR_ALL_ELEMENTS_IN_ARRAY2D(z)
        {
            XX(r) = j;
            YY(r) = i;
            double G = unnormalizedGaussian2D(r, mu, sigmainv);
            z(i, j) = I(i, j) - a * G;
        }
        compute_hist(z, hist, 100);
        b = hist.percentil(5);
    }
}

/* Fourier-Bessel decomposition. ------------------------------------------- */
void fourierBesselDecomposition(const MultidimArray<double> &img_in,
                                double r2, int k1, int k2)
{
    img_in.checkDimension(2);

    for (int k = k1; k <= k2; k++)
    {
        int k_1 = k - 1;

        // Compute h and a,b coefficients
        double h = 0;
        double my5 = 0;
        if (k_1 != 0)
        {
            double my = 1 + PI * r2 / 2 / k_1;
            double my4 = my * k_1;
            double ntot = 4 * my4;
            h = 2 * PI / ntot;
        }
        else
        {
            double my = 1 + PI * r2 / 2;
            double my4 = my;
            double ntot = 4 * my4;
            h = 2 * PI / ntot;
        }

        MultidimArray<double> sine(CEIL(my5));
        FOR_ALL_ELEMENTS_IN_ARRAY1D(sine)
        sine(i) = sin((i + 1) * h);

    }
}

/* Harmonic decomposition. ------------------------------------------------- */
//void harmonicDecomposition(const MultidimArray<double> &img_in,
//                           MultidimArray<double> &v_out)
//{}


/* Shah energy ------------------------------------------------------------- */
/* This function computes the current functional energy */
double Shah_energy(const MultidimArray<double> &img,
                   const MultidimArray<double> &surface_strength,
                   const MultidimArray<double> &edge_strength, double K,
                   const Matrix1D<double> &W)
{
    img.checkDimension(2);

    int Ydim1 = YSIZE(img) - 1;
    int Xdim1 = XSIZE(img) - 1;

    double Kinv = 1.0 / K;

    /* Calculate surface energy */
    double E1 = 0.0;
    double E2 = 0.0;
    double E3 = 0.0;
    double E4 = 0.0;
    double w0 = VEC_ELEM(W,0);
    double w1 = VEC_ELEM(W,1);
    double w2 = VEC_ELEM(W,2);
    double w3 = VEC_ELEM(W,3);
    for (int i = 1; i < Ydim1; i++)
    {
        int ip1 = i + 1;
        int im1 = i - 1;
        for (int j = 1; j < Xdim1; j++)
        {
            int jp1 = j + 1;
            int jm1 = j - 1;
            /* Calculate data matching terms */
            double D = dAij(img, i, j);
            double F = dAij(surface_strength, i, j);
            double S = dAij(edge_strength, i, j);
            double diff = D - F;
            E1 += w0 * diff * diff;
            E3 += w2 * K * S * S;

            /* Calculate first derivative terms */
            double Fx = 0.5
                        * (dAij(surface_strength, i, jp1)
                           - dAij(surface_strength, i, jm1));
            double Fy = 0.5
                        * (dAij(surface_strength, ip1, j)
                           - dAij(surface_strength, im1, j));
            double Sx =
                0.5
                * (dAij(edge_strength, i, jp1)
                   - dAij(edge_strength, i, jm1));
            double Sy =
                0.5
                * (dAij(edge_strength, ip1, j)
                   - dAij(edge_strength, im1, j));
            double S_1 = 1 - S;
            E2 += w1 * S_1 * S_1 * (Fx * Fx + Fy * Fy);
            E4 += w3 * Kinv * (Sx * Sx + Sy * Sy);
        }
    }

    return E1 + E2 + E3 + E4; // Total energy
}

/* Update Surface Shah ----------------------------------------------------- */
/* This routine performs one update to the edge estimate based
 on a finite differences solution to the following equation:
 0 = dE/df = dF/df - d(dF/dfx)/dx - d(dF/dfy)/dy
 + dd(dF/dfxx)/dxx + dd(dF/dfxy)/dxy + dd(dF/dfyy)/dyy */
double Update_surface_Shah(MultidimArray<double> &img,
                           MultidimArray<double> &surface_strength,
                           MultidimArray<double> &edge_strength, const Matrix1D<double> &W)
{
    img.checkDimension(2);

    double Diff = 0.0;
    double Norm = 0.0;
    size_t Ydim1 = YSIZE(img) - 1;
    size_t Xdim1 = XSIZE(img) - 1;

    /* Update surface estimate */
    double w0 = VEC_ELEM(W,0);
    double w1 = VEC_ELEM(W,1);
    for (size_t i = 1; i < Ydim1; i++)
    {
        int ip1 = i + 1;
        int im1 = i - 1;
        for (size_t j = 1; j < Xdim1; j++)
        {
            int jp1 = j + 1;
            int jm1 = j - 1;
            /* Calculate edge partial derivative terms */
            double S = dAij(edge_strength, i, j);
            double Sx =
                0.5
                * (dAij(edge_strength, i, jp1)
                   - dAij(edge_strength, i, jm1));
            double Sy =
                0.5
                * (dAij(edge_strength, ip1, j)
                   - dAij(edge_strength, im1, j));

            double nS = 1 - S;
            double nS2 = nS * nS;

            /* Calculate surface partial derivative terms (excluding central pixel) */
            double F;
            double D;
            F = D = dAij(img, i, j);
            double SS_i_jp1 = dAij(surface_strength, i, jp1);
            double SS_i_jm1 = dAij(surface_strength, i, jm1);
            double SS_ip1_j = dAij(surface_strength, ip1, j);
            double SS_im1_j = dAij(surface_strength, im1, j);
            double Fx = 0.5 * (SS_i_jp1 - SS_i_jm1);
            double Fy = 0.5 * (SS_ip1_j - SS_im1_j);
            double Fxx = SS_i_jp1 + SS_i_jm1;
            double Fyy = SS_ip1_j + SS_im1_j;

            /* Calculate surface partial derivative weights */
            double wFx = 4 * w1 * nS * Sx;
            double wFy = 4 * w1 * nS * Sy;
            double wFxx = -2 * w1 * nS2;
            double wFyy = -2 * w1 * nS2;

            /* Calculate new surface value */
            double Constant = -2 * w0 * D;
            double Central = -2 * w0 + 2 * wFxx + 2 * wFyy;
            double Neighbors = wFx * Fx + wFy * Fy + wFxx * Fxx + wFyy * Fyy;

            if (fabs(Central) > XMIPP_EQUAL_ACCURACY)
                F = (Constant + Neighbors) / Central;
            F = CLIP(F, 0, 1);

            // Compute the difference.
            double SS_i_j = dAij(surface_strength, i, j);
            Diff += fabs(SS_i_j - F);
            Norm += fabs(SS_i_j);

            // Update the new value.
            dAij(surface_strength, i, j) = F;
        }
    }
    return Diff / Norm; // Return the relative difference.
}

/* Update Edge Shah -------------------------------------------------------- */
/* This routine performs one update to the edge estimate based
 on a finite differences solution to the following equation:
 0 = dE/ds = dF/ds - d(dF/dsx)/dx - d(dF/dsy)/dy */
double Update_edge_Shah(MultidimArray<double> &img,
                        MultidimArray<double> &surface_strength,
                        MultidimArray<double> &edge_strength, double K,
                        const Matrix1D<double> &W)
{
    img.checkDimension(2);

    double Diff = 0.0;
    double Norm = 0.0;
    size_t Ydim1 = YSIZE(img) - 1;
    size_t Xdim1 = XSIZE(img) - 1;
    double Kinv = 1.0 / K;

    /* Update edge estimate */
    double w1 = VEC_ELEM(W,1);
    double w2 = VEC_ELEM(W,2);
    double w3 = VEC_ELEM(W,3);
    for (size_t i = 1; i < Ydim1; i++)
    {
        int ip1 = i + 1;
        int im1 = i - 1;
        for (size_t j = 1; j < Xdim1; j++)
        {
            int jp1 = j + 1;
            int jm1 = j - 1;

            /* Calculate first and second derivative terms */
            double Fx = 0.5
                        * (dAij(surface_strength, i, jp1)
                           - dAij(surface_strength, i, jm1));
            double Fy = 0.5
                        * (dAij(surface_strength, ip1, j)
                           - dAij(surface_strength, im1, j));
            double Constant = w1 * (Fx * Fx + Fy * Fy);

            /* Calculate weights for central pixel and neighbors */
            double Central = w2 * K + w3 * Kinv * 4;
            double Neighbors = w3 * Kinv
                               * (dAij(edge_strength, im1, j) + dAij(edge_strength, ip1, j)
                                  + dAij(edge_strength, i, jm1)
                                  + dAij(edge_strength, i, jp1));

            /* Calculate new S value */
            double Old_edge_strength = dAij(edge_strength, i, j);
            double S = (Constant + Neighbors) / (Constant + Central);
            if (S < 0)
                dAij(edge_strength, i, j) *= 0.5;
            else if (S > 1)
                dAij(edge_strength, i, j) = 0.5
                                            * (dAij(edge_strength, i, j) + 1);
            else
                dAij(edge_strength, i, j) = S;

            // Compute the difference.
            Diff += fabs(dAij(edge_strength, i, j) - Old_edge_strength);
            Norm += fabs(Old_edge_strength);
        }
    }
    return Diff / Norm; // Return the relative difference.
}

/* Smoothing Shah ---------------------------------------------------------- */
constexpr double SHAH_CONVERGENCE_THRESHOLD = 0.0001;
void smoothingShah(MultidimArray<double> &img,
                   MultidimArray<double> &surface_strength,
                   MultidimArray<double> &edge_strength, const Matrix1D<double> &W,
                   int OuterLoops, int InnerLoops, int RefinementLoops,
                   bool adjust_range)
{

    img.checkDimension(2);

    typeCast(img, surface_strength);
    if (adjust_range)
        surface_strength.rangeAdjust(0, 1);
    edge_strength.resizeNoCopy(img);

    for (int k = 1; k <= RefinementLoops; k++)
    {
        // Initialize Edge Image.
        edge_strength.initZeros();

        double diffsurface = MAXFLOAT; // Reset surface difference
        for (int i = 0;
             ((i < OuterLoops) && OuterLoops)
             || ((diffsurface > SHAH_CONVERGENCE_THRESHOLD)
                 && !OuterLoops); i++)
        {

            /* std::cout << "Iteration ..." << i+1;*/
            /* Iteratively update surface estimate */
            for (int j = 0; j < InnerLoops; j++)
                diffsurface = Update_surface_Shah(img, surface_strength,
                                                  edge_strength, W);

            /* Iteratively update edge estimate */
            for (int j = 0; j < InnerLoops; j++)
                Update_edge_Shah(img, surface_strength, edge_strength, k, W);

            /* Calculate new functional energy */
            Shah_energy(img, surface_strength, edge_strength, k, W);
        }
    }
}

/* Tomographic diffusion --------------------------------------------------- */
//#define DEBUG
double tomographicDiffusion(MultidimArray<double>& V,
                            const Matrix1D<double>& alpha, double lambda)
{
    V.checkDimension(3);

    double alphax = XX(alpha);
    double alphay = YY(alpha);
    double alphaz = ZZ(alpha);
    double diffx;
    double diffy;
    double diffz;

    // Compute regularization error
    double regError = 0;
    for (size_t z = 1; z < ZSIZE(V) - 1; z++)
        for (size_t y = 1; y < YSIZE(V) - 1; y++)
            for (size_t x = 1; x < XSIZE(V) - 1; x++)
            {
                diffx = DIRECT_A3D_ELEM(V,z,y,x+1) - DIRECT_A3D_ELEM(V,z,y,x-1);
                diffy = DIRECT_A3D_ELEM(V,z,y+1,x) - DIRECT_A3D_ELEM(V,z,y-1,x);
                diffz = DIRECT_A3D_ELEM(V,z+1,y,x) - DIRECT_A3D_ELEM(V,z-1,y,x);
                regError += sqrt(
                                alphax * diffx * diffx + alphay * diffy * diffy
                                + alphaz * diffz * diffz);
            }
    regError *= 0.5;

    // Compute the gradient of the regularization error
    MultidimArray<double> gradient;
    gradient.initZeros(V);
    for (size_t z = 2; z < ZSIZE(V) - 2; z++)
        for (size_t y = 2; y < YSIZE(V) - 2; y++)
            for (size_t x = 2; x < XSIZE(V) - 2; x++)
            {
                // First term
                double V000 = DIRECT_A3D_ELEM(V,z,y,x);
                double V_200 = DIRECT_A3D_ELEM(V,z,y,x-2);
                double V_110 = DIRECT_A3D_ELEM(V,z,y+1,x-1);
                double V_1_10 = DIRECT_A3D_ELEM(V,z,y-1,x-1);
                double V_101 = DIRECT_A3D_ELEM(V,z+1,y,x-1);
                double V_10_1 = DIRECT_A3D_ELEM(V,z-1,y,x-1);
                diffx = V000 - V_200;
                diffy = V_110 - V_1_10;
                diffz = V_101 - V_10_1;
                double t1 = diffx
                            / sqrt(
                                alphax * diffx * diffx + alphay * diffy * diffy
                                + alphaz * diffz * diffz);

                // Second term
                double V200 = DIRECT_A3D_ELEM(V,z,y,x+2);
                double V110 = DIRECT_A3D_ELEM(V,z,y+1,x+1);
                double V1_10 = DIRECT_A3D_ELEM(V,z,y-1,x+1);
                double V101 = DIRECT_A3D_ELEM(V,z+1,y,x+1);
                double V10_1 = DIRECT_A3D_ELEM(V,z-1,y,x+1);
                diffx = V200 - V000;
                diffy = V110 - V1_10;
                diffz = V101 - V10_1;
                double t2 = diffx
                            / sqrt(
                                alphax * diffx * diffx + alphay * diffy * diffy
                                + alphaz * diffz * diffz);

                // Third term
                double V0_20 = DIRECT_A3D_ELEM(V,z,y-2,x);
                double V0_11 = DIRECT_A3D_ELEM(V,z+1,y-1,x);
                double V0_1_1 = DIRECT_A3D_ELEM(V,z-1,y-1,x);
                diffx = V1_10 - V_1_10;
                diffy = V000 - V0_20;
                diffz = V0_11 - V0_1_1;
                double t3 = diffy
                            / sqrt(
                                alphax * diffx * diffx + alphay * diffy * diffy
                                + alphaz * diffz * diffz);

                // Fourth term
                double V020 = DIRECT_A3D_ELEM(V,z,y+2,x);
                double V011 = DIRECT_A3D_ELEM(V,z+1,y+1,x);
                double V01_1 = DIRECT_A3D_ELEM(V,z-1,y+1,x);
                diffx = V110 - V_110;
                diffy = V020 - V000;
                diffz = V011 - V01_1;
                double t4 = diffy
                            / sqrt(
                                alphax * diffx * diffx + alphay * diffy * diffy
                                + alphaz * diffz * diffz);

                // Fifth term
                double V00_2 = DIRECT_A3D_ELEM(V,z-2,y,x);
                diffx = V10_1 - V_10_1;
                diffy = V01_1 - V0_1_1;
                diffz = V000 - V00_2;
                double t5 = diffz
                            / sqrt(
                                alphax * diffx * diffx + alphay * diffy * diffy
                                + alphaz * diffz * diffz);

                // Sixth term
                double V002 = DIRECT_A3D_ELEM(V,z+2,y,x);
                diffx = V101 - V_101;
                diffy = V011 - V0_11;
                diffz = V002 - V000;
                double t6 = diffz
                            / sqrt(
                                alphax * diffx * diffx + alphay * diffy * diffy
                                + alphaz * diffz * diffz);

                // Compute gradient
                DIRECT_A3D_ELEM(gradient,z,y,x) = 0.5
                                                  * (alphax * (t1 - t2) + alphay * (t3 - t4)
                                                     + alphaz * (t5 - t6));
            }
#ifdef DEBUG
    VolumeXmipp save;
    save()=V;
    save.write("PPPvolume.vol");
    save()=gradient;
    save.write("PPPgradient.vol");
    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
#endif

    // Update volume
    for (size_t z = 2; z < ZSIZE(V) - 2; z++)
        for (size_t y = 2; y < YSIZE(V) - 2; y++)
            for (size_t x = 2; x < XSIZE(V) - 2; x++)
                DIRECT_A3D_ELEM(V,z,y,x) -= lambda
                                            * DIRECT_A3D_ELEM(gradient,z,y,x);

    // Finish
    return regError;
}
#undef DEBUG

/* Rotational invariant moments -------------------------------------------- */
void rotationalInvariantMoments(const MultidimArray<double> &img,
                                const MultidimArray<int> *mask, MultidimArray<double> &v_out)
{
    img.checkDimension(2);

    // Prepare some variables
    double m_11 = 0;
    double m_02 = 0;
    double m_20 = 0;
    double m_12 = 0;
    double m_21 = 0;
    double m_03 = 0;
    double m_30 = 0;
    double normalize_x = 2.0 / XSIZE(img);
    double normalize_y = 2.0 / YSIZE(img);

    // Compute low-level moments
    FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
    {
        if (mask != nullptr)
            if (!(*mask)(i, j))
                continue;
        double val = img(i, j);
        double dx = j * normalize_x;
        double dy = i * normalize_y;
        double dx2 = dx * dx;
        double dx3 = dx2 * dx;
        double dy2 = dy * dy;
        double dy3 = dy2 * dy;
        // m_00+=val;
        m_11 += val * dx * dy;
        m_02 += val * dy2;
        m_20 += val * dx2;
        m_12 += val * dx * dy2;
        m_21 += val * dx2 * dy;
        m_03 += val * dy3;
        m_30 += val * dx3;
    }
    //m_11/=m_00; m_02/=m_00; m_20/=m_00;
    //m_12/=m_00; m_21/=m_00; m_03/=m_00; m_30/=m_00;

    // Compute high-level rotational invariant moments
    v_out.resize(5);
    v_out(0) = m_20 + m_02;
    v_out(1) = (m_20 - m_02) * (m_20 - m_02) + 4 * m_11 * m_11;
    v_out(2) = (m_30 - 3 * m_12) * (m_30 - 3 * m_12)
               + (3 * m_21 - m_03) * (3 * m_21 - m_03);
    v_out(3) = (m_30 + m_12) * (m_30 + m_12) + (m_21 + m_03) * (m_21 + m_03);
    v_out(4) =
        (m_30 - 3 * m_12) * (m_30 + m_12)
        * ((m_30 + m_12) * (m_30 + m_12)
           - 3 * (m_21 + m_03) * (m_21 + m_03))
        + (3 * m_21 - m_03) * (m_21 + m_03)
        * (3 * (m_30 + m_12) * (m_30 + m_12)
           - (m_21 + m_03) * (m_21 + m_03));
    /*
     v_out( 5)=(m_20+m_02)*(
     (m_30+m_12)*(m_30+m_12)
     -3*(m_21+m_03)*(m_21+m_03))
     +4*m_11*(m_30+m_12)*(m_03+m_21);
     v_out( 6)=(3*m_21-m_03)*(m_12+m_30)*(
     (m_30+m_12)*(m_30+m_12)
     -3*(m_21+m_03)*(m_21+m_03))
     -(m_30-3*m_12)*(m_12+m_03)*(
     3*(m_30+m_12)*(m_30+m_12)
     -(m_21+m_03)*(m_21+m_03));
     */
}

/* Inertia moments --------------------------------------------------------- */
void inertiaMoments(const MultidimArray<double> &img,
                    const MultidimArray<int> *mask, Matrix1D<double> &v_out,
                    Matrix2D<double> &u)
{
    img.checkDimension(2);

    // Prepare some variables
    double m_11 = 0;
    double m_02 = 0;
    double m_20 = 0;
    double normalize_x = 2.0 / XSIZE(img);
    double normalize_y = 2.0 / YSIZE(img);

    // Compute low-level moments
    FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
    {
        if (mask != nullptr)
            if (!(*mask)(i, j))
                continue;
        double val = img(i, j);
        double dx = j * normalize_x;
        double dy = i * normalize_y;
        double dx2 = dx * dx;
        double dy2 = dy * dy;
        m_11 += val * dx * dy;
        m_02 += val * dy2;
        m_20 += val * dx2;
    }

    // Compute the eigen values of the inertia matrix
    // [m_02 m_11
    //  m_11 m_20]
    Matrix2D<double> A(2, 2);
    A(0, 0) = m_02;
    A(0, 1) = A(1, 0) = m_11;
    A(1, 1) = m_20;
    Matrix2D<double> v;
    svdcmp(A, u, v_out, v);
    v_out = v_out.sort();
}

/* Fill triangle ----------------------------------------------------------- */
void fillTriangle(MultidimArray<double> &img, int *tx, int *ty, double color)
{
    img.checkDimension(2);

    /*
     * Order in y values
     */
    int y1 = 0;
    while (!y1)
    {
        y1 = 1;
        for (int y2 = 0; y2 < 2; y2++)
        {
            if (ty[y2] > ty[y2 + 1]
                || (ty[y2] == ty[y2 + 1] && tx[y2] < tx[y2 + 1]))
            {
                int x1 = ty[y2];
                ty[y2] = ty[y2 + 1];
                ty[y2 + 1] = x1;
                x1 = tx[y2];
                tx[y2] = tx[y2 + 1];
                tx[y2 + 1] = x1;
                y1 = 0;
            }
        }
    }

    int dx1 = tx[1] - tx[0];
    int dx2 = tx[2] - tx[0];
    int dy1 = ty[1] - ty[0];
    int dy2 = ty[2] - ty[0];

    int sx1 = SGN0(dx1);
    int sx2 = SGN0(dx2);
    int sy1 = SGN0(dy1);

    int ix1 = ABS(dx1);
    int ix2 = ABS(dx2);
    int iy1 = ABS(dy1);
    int iy2 = ABS(dy2);

    int inc1 = XMIPP_MAX(ix1, iy1);
    int inc2 = XMIPP_MAX(ix2, iy2);

    int x1;
    int x2;
    int y2;
    int xl;
    int xr;
    x1 = x2 = y1 = y2 = 0;
    xl = xr = tx[0];
    int y = ty[0];

    while (y != ty[1])
    {
        int doit1 = 0;
        int doit2 = 0;

        while (!doit1)
        { /* Wait until y changes */
            x1 += ix1;
            y1 += iy1;
            if (x1 > inc1)
            {
                x1 -= inc1;
                xl += sx1;
            }
            if (y1 > inc1)
            {
                y1 -= inc1;
                y += sy1;
                doit1 = 1;
            }
        }

        while (!doit2)
        { /* Wait until y changes */
            x2 += ix2;
            y2 += iy2;
            if (x2 > inc2)
            {
                x2 -= inc2;
                xr += sx2;
            }
            if (y2 > inc2)
            {
                y2 -= inc2;
                doit2 = 1;
            }
        }

        for (int myx = xl; myx <= xr; myx++)
            img(y, myx) = color;
    }

    dx1 = tx[2] - tx[1];
    dy1 = ty[2] - ty[1];

    sx1 = SGN0(dx1);
    sy1 = SGN0(dy1);

    ix1 = ABS(dx1);
    iy1 = ABS(dy1);

    inc1 = XMIPP_MAX(ix1, iy1);
    xl = tx[1];
    x1 = 0;

    while (y != ty[2])
    {
        int doit1 = 0;
        int doit2 = 0;

        while (!doit1)
        { /* Wait until y changes */
            x1 += ix1;
            y1 += iy1;
            if (x1 > inc1)
            {
                x1 -= inc1;
                xl += sx1;
            }
            if (y1 > inc1)
            {
                y1 -= inc1;
                y += sy1;
                doit1 = 1;
            }
        }

        while (!doit2)
        { /* Wait until y changes */
            x2 += ix2;
            y2 += iy2;
            if (x2 > inc2)
            {
                x2 -= inc2;
                xr += sx2;
            }
            if (y2 > inc2)
            {
                y2 -= inc2;
                doit2 = 1;
            }
        }

        for (int myx = xl; myx <= xr; myx++)
            img(y, myx) = color;
    }
}

/* Local thresholding ------------------------------------------------------ */
void localThresholding(MultidimArray<double> &img, double C, double dimLocal,
                       MultidimArray<int> &result, MultidimArray<int> *mask)
{

    // Convolve the input image with the kernel
    MultidimArray<double> convolved;
    convolved.initZeros(img);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(convolved)
    {
        if (mask != nullptr)
            if (!(*mask)(i, j))
                continue;
        int ii0 = XMIPP_MAX(STARTINGY(convolved), FLOOR(i - dimLocal));
        int jj0 = XMIPP_MAX(STARTINGX(convolved), FLOOR(j - dimLocal));
        int iiF = XMIPP_MIN(FINISHINGY(convolved), CEIL(i + dimLocal));
        int jjF = XMIPP_MIN(FINISHINGX(convolved), CEIL(j + dimLocal));
        double N = 0;
        for (int ii = ii0; ii <= iiF; ii++)
            for (int jj = jj0; jj <= jjF; jj++)
            {
                if (mask == nullptr)
                {
                    convolved(i, j) += img(ii, jj);
                    ++N;
                }
                else if ((*mask)(i, j))
                {
                    convolved(i, j) += img(ii, jj);
                    ++N;
                }
            }
        if (N != 0)
            convolved(i, j) /= N;
    }

    // Subtract the original from the convolved image and threshold
    result.initZeros(img);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(convolved)
    {
        if (mask != nullptr)
            if (!(*mask)(i, j))
                continue;
        if (img(i, j) - convolved(i, j) > C)
            result(i, j) = 1;
    }
}

/* Center translationally -------------------------------------------------- */
void centerImageTranslationally(MultidimArray<double> &I, CorrelationAux &aux)
{
    I.checkDimension(2);

    MultidimArray<double> Ix = I;
    Ix.selfReverseX();
    Ix.setXmippOrigin();
    MultidimArray<double> Iy = I;
    Iy.selfReverseY();
    Iy.setXmippOrigin();
    MultidimArray<double> Ixy = Ix;
    Ixy.selfReverseY();
    Ixy.setXmippOrigin();

    double meanShiftX = 0;
    double meanShiftY = 0;
    double shiftX;
    double shiftY;
    bestNonwrappingShift(I, Ix, meanShiftX, meanShiftY, aux);
    bestNonwrappingShift(I, Iy, shiftX, shiftY, aux);
    meanShiftX += shiftX;
    meanShiftY += shiftY;
    bestNonwrappingShift(I, Ixy, shiftX, shiftY, aux);
    meanShiftX += shiftX;
    meanShiftY += shiftY;
    meanShiftX /= 3;
    meanShiftY /= 3;

    Matrix1D<double> shift(2);
    VECTOR_R2(shift, -meanShiftX, -meanShiftY);
    MultidimArray<double> aux2 = I;
    translate(3, I, aux2, shift);
    //I.selfTranslateBSpline(3,shift);
}

/* Center rotationally ----------------------------------------------------- */
void centerImageRotationally(MultidimArray<double> &I,
                             RotationalCorrelationAux &aux)
{
    I.checkDimension(2);

    MultidimArray<double> Ix = I;
    Ix.selfReverseX();
    Ix.setXmippOrigin();

    Polar_fftw_plans *plans = nullptr;
    Polar<std::complex<double> > polarFourierI;
    Polar<std::complex<double> > polarFourierIx;
    polarFourierTransform<true>(Ix, polarFourierIx, false, XSIZE(Ix) / 5,
                                    XSIZE(Ix) / 2, plans);
    polarFourierTransform<true>(I, polarFourierI, true, XSIZE(I) / 5,
                                    XSIZE(I) / 2, plans);

    MultidimArray<double> rotationalCorr;
    rotationalCorr.resize(2 * polarFourierI.getSampleNoOuterRing() - 1);
    aux.local_transformer.setReal(rotationalCorr);
    double bestRot = best_rotation(polarFourierIx, polarFourierI, aux);

    MultidimArray<double> auxI = I;
    rotate(3, I, auxI, -bestRot / 2, xmipp_transformation::WRAP);
    //I.selfRotateBSpline(3,-bestRot/2,WRAP);
}

/* Center both rotationally and translationally ---------------------------- */
//#define DEBUG
Matrix2D<double> centerImage(MultidimArray<double> &I, CorrelationAux &aux,
                 RotationalCorrelationAux &aux2, int Niter, bool limitShift)
{
    I.checkDimension(2);

    I.setXmippOrigin();
    double avg = I.computeAvg();
    I -= avg;

    MultidimArray<double> Ix;
    MultidimArray<double> Iy;
    MultidimArray<double> Ixy;
    MultidimArray<double> Iaux;
    Matrix2D<double> A;
    A.initIdentity(3);
    Iaux = I;

    MultidimArray<int> mask;
    mask.initZeros(I);
    BinaryCircularMask(mask, XSIZE(I) / 2);

    MultidimArray<double> lineY;
    MultidimArray<double> lineX;
    lineY.initZeros(YSIZE(I));
    STARTINGX(lineY) = STARTINGY(I);
    lineX.initZeros(XSIZE(I));
    STARTINGX(lineX) = STARTINGX(I);

    Polar_fftw_plans *plans = nullptr;
    Polar<std::complex<double> > polarFourierI;
    Polar<std::complex<double> > polarFourierIx;
    MultidimArray<double> rotationalCorr;
    Matrix2D<double> R;
    for (int i = 0; i < Niter; i++)
    {
        // Mask Iaux
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        if (!A2D_ELEM(mask,i,j))
            A2D_ELEM(Iaux,i,j) = 0;

        // Center translationally
        Ix = Iaux;
        Ix.selfReverseX();
        Ix.setXmippOrigin();
        Iy = Iaux;
        Iy.selfReverseY();
        Iy.setXmippOrigin();
        Ixy = Ix;
        Ixy.selfReverseY();
        Ixy.setXmippOrigin();

        double meanShiftX = 0;
        double meanShiftY = 0;
        double shiftX;
        double shiftY;
        double Nx = 0;
        double Ny = 0;
        bestNonwrappingShift(Iaux, Ix, shiftX, shiftY, aux);
#ifdef DEBUG

        ImageXmipp save;
        save()=Ix;
        save.write("PPPx.xmp");
        std::cout << "con Ix: " << shiftX << " " << shiftY << std::endl;
#endif

        if (fabs(shiftX) < XSIZE(I) / 3 || !limitShift)
        {
            meanShiftX += shiftX;
            Nx++;
        }
        if (fabs(shiftY) < YSIZE(I) / 3 || !limitShift)
        {
            meanShiftY += shiftY;
            Ny++;
        }
        bestNonwrappingShift(Iaux, Iy, shiftX, shiftY, aux);
#ifdef DEBUG

        save()=Iy;
        save.write("PPPy.xmp");
        std::cout << "con Iy: " << shiftX << " " << shiftY << std::endl;
#endif

        if (fabs(shiftX) < XSIZE(I) / 3 || !limitShift)
        {
            meanShiftX += shiftX;
            Nx++;
        }
        if (fabs(shiftY) < YSIZE(I) / 3 || !limitShift)
        {
            meanShiftY += shiftY;
            Ny++;
        }
        bestNonwrappingShift(Iaux, Ixy, shiftX, shiftY, aux);
#ifdef DEBUG

        save()=Ixy;
        save.write("PPPxy.xmp");
        std::cout << "con Ixy: " << shiftX << " " << shiftY << std::endl;
#endif

        if (fabs(shiftX) < XSIZE(I) / 3 || !limitShift)
        {
            meanShiftX += shiftX;
            Nx++;
        }
        if (fabs(shiftY) < YSIZE(I) / 3 || !limitShift)
        {
            meanShiftY += shiftY;
            Ny++;
        }
        if (Nx > 0)
            meanShiftX /= Nx;
        if (Ny > 0)
            meanShiftY /= Ny;

        MAT_ELEM(A,0,2) += -meanShiftX / 2;
        MAT_ELEM(A,1,2) += -meanShiftY / 2;
        Iaux.initZeros();
        applyGeometry(xmipp_transformation::LINEAR, Iaux, I, A, xmipp_transformation::IS_NOT_INV, xmipp_transformation::WRAP);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        if (!A2D_ELEM(mask,i,j))
            A2D_ELEM(Iaux,i,j) = 0;

#ifdef DEBUG

        std::cout << "Iter " << i << std::endl;
        std::cout << "shift=" << -meanShiftX << "," << -meanShiftY << std::endl;
        save()=I;
        save.write("PPP.xmp");
        save()=Iaux;
        save.write("PPPshift.xmp");
#endif

        // Center rotationally
        Ix = Iaux;
        Ix.selfReverseX();
        Ix.setXmippOrigin();

        polarFourierTransform<true>(Iaux, polarFourierI, true, XSIZE(I) / 5,
                                        XSIZE(I) / 2, plans);
        rotationalCorr.resizeNoCopy(
            2 * polarFourierI.getSampleNoOuterRing() - 1);
        aux2.local_transformer.setReal(rotationalCorr);

        polarFourierTransform<true>(Ix, polarFourierIx, false,
                                        XSIZE(Ix) / 5, XSIZE(Ix) / 2, plans);
        double bestRot = best_rotation(polarFourierIx, polarFourierI, aux2);
        bestRot = realWRAP(bestRot,0,180);
        if (bestRot > 90)
            bestRot = bestRot - 180;

        rotation2DMatrix(bestRot / 2, R);
        A = R * A;
        Iaux.initZeros();
        applyGeometry(xmipp_transformation::LINEAR, Iaux, I, A, xmipp_transformation::IS_NOT_INV, xmipp_transformation::WRAP);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        if (!A2D_ELEM(mask,i,j))
            A2D_ELEM(Iaux,i,j) = 0;

#ifdef DEBUG

        std::cout << "rot=" << -bestRot/2 << std::endl;
        save()=Iaux;
        save.write("PPProt.xmp");
#endif

        // Remove horizontal/vertical ambiguity
        lineX.initZeros();
        lineY.initZeros();
        FOR_ALL_ELEMENTS_IN_ARRAY2D(Iaux)
        {
            double val = A2D_ELEM(Iaux,i,j);
            if (j == 0)
                lineY(i) = val;
            else if (lineY(i) < val)
                lineY(i) = val;
            if (i == 0)
                lineX(j) = val;
            else if (lineX(j) < val)
                lineX(j) = val;
        }

        double thX = lineX.computeMin()
                     + 0.75 * (lineX.computeMax() - lineX.computeMin());
        double thY = lineY.computeMin()
                     + 0.75 * (lineY.computeMax() - lineY.computeMin());
        int x0 = STARTINGX(lineX);
        while (lineX(x0) < thX)
            x0++;
        int y0 = STARTINGX(lineY);
        while (lineY(y0) < thY)
            y0++;
        int xF = FINISHINGX(lineX);
        while (lineX(xF) < thX)
            xF--;
        int yF = FINISHINGX(lineY);
        while (lineY(yF) < thY)
            yF--;
        if ((xF - x0) > (yF - y0))
        {
            rotation2DMatrix(-90., R);
            A = R * A;
        }
        applyGeometry(xmipp_transformation::LINEAR, Iaux, I, A, xmipp_transformation::IS_NOT_INV, xmipp_transformation::WRAP);
#ifdef DEBUG

        lineX.write("PPPlineX.txt");
        lineY.write("PPPlineY.txt");
        std::cout << "dev X=" << xF-x0 << std::endl;
        std::cout << "dev Y=" << yF-y0 << std::endl;
        save()=Iaux;
        save.write("PPPhorver.xmp");
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
#endif

    }
    applyGeometry(xmipp_transformation::BSPLINE3, Iaux, I, A, xmipp_transformation::IS_NOT_INV, xmipp_transformation::WRAP);
    I = Iaux;
    I += avg;
    delete plans;
    return A;
}
#undef DEBUG

/** Force positive -------------------------------------------------------- */
void forcePositive(MultidimArray<double> &V)
{
    MultidimArray<char> mask(ZSIZE(V), YSIZE(V), XSIZE(V));
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V)
    {
        double x = DIRECT_MULTIDIM_ELEM(V, n);
        DIRECT_MULTIDIM_ELEM(mask, n) = (x <= 0);
    }
    boundMedianFilter(V, mask);
}

void computeEdges(const MultidimArray<double>& vol,
                  MultidimArray<double> &vol_edge)
{
    MultidimArray<double> BSpline_coefs;
    BSpline_coefs.initZeros(vol);
    produceSplineCoefficients(3, BSpline_coefs, vol);
    vol_edge.initZeros(vol);

    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
    {
        double V_dx;
        double V_dy;
        double V_dz;
        V_dx = interpolatedElementBSplineDiffX(BSpline_coefs, j, i, k, 3);
        V_dy = interpolatedElementBSplineDiffY(BSpline_coefs, j, i, k, 3);
        V_dz = interpolatedElementBSplineDiffZ(BSpline_coefs, j, i, k, 3);
        A3D_ELEM(vol_edge,k,i,j) = sqrt(
                                       (V_dx * V_dx) + (V_dy * V_dy) + (V_dz * V_dz));

    }
}

void forceDWTSparsity(MultidimArray<double> &V, double eps)
{
	int size0=XSIZE(V);
	int sizeF=(int)NEXT_POWER_OF_2(size0);
    selfScaleToSize(xmipp_transformation::BSPLINE3,V,sizeF,sizeF,sizeF);
    MultidimArray<double> vol_wavelets;
    MultidimArray<double> vol_wavelets_abs;
    set_DWT_type(DAUB12);
    DWT(V,vol_wavelets);
    vol_wavelets_abs=vol_wavelets;
    vol_wavelets_abs.selfABS();
    double *begin=MULTIDIM_ARRAY(vol_wavelets_abs);
    double *end=MULTIDIM_ARRAY(vol_wavelets_abs)+MULTIDIM_SIZE(vol_wavelets_abs);
    std::sort(begin,end);
    double threshold1=DIRECT_MULTIDIM_ELEM(vol_wavelets_abs,
                                           (long int)((1.0-eps)*MULTIDIM_SIZE(vol_wavelets_abs)));
    vol_wavelets.threshold("abs_below", threshold1, 0.0);
    IDWT(vol_wavelets,V);
    selfScaleToSize(xmipp_transformation::BSPLINE3,V,size0, size0, size0);
}

/////////////// FILTERS IMPLEMENTATIONS /////////////////

/** Define the parameters for use inside an Xmipp program */
void BadPixelFilter::defineParams(XmippProgram * program)
{
    program->addParamsLine("== Bad pixels ==");
    program->addParamsLine(
        "  [ --bad_pixels <type>]            : Applied filters on bad pixels of the image.");
    program->addParamsLine("         where <type>  ");
    program->addParamsLine(
        "            negative              : Applied at those negative values. Positive values are untouched.");
    program->addParamsLine(
        "            mask <mask_file>      : Applied at those pixels given by mask.");
    program->addParamsLine(
        "            outliers <factor>     : Applied at those pixels out of the range [mean - factor*std, mean + factor*std].");
    program->addParamsLine("         alias -b; ");
}

/** Read from program command line */
void BadPixelFilter::readParams(XmippProgram * program)
{
    type = NEGATIVE;
    // Check operation to do
    String typeStr = program->getParam("--bad_pixels");
    if (typeStr == "negative")
        ; //nothing to do type already equal to NEGATIVE
    else if (typeStr == "mask")
    {
        mask = new Image<char>;
        mask->read(program->getParam("--bad_pixels", "mask"));
        type = MASK;
    }
    else if (typeStr == "outliers")
    {
        factor = program->getDoubleParam("--bad_pixels", "outliers");
        type = OUTLIER;
    }
}

/** Apply the filter to an image or volume*/
void BadPixelFilter::apply(MultidimArray<double> &img)
{
    switch (type)
    {
    case NEGATIVE:
        forcePositive(img);
        break;
    case MASK:
        boundMedianFilter(img, mask->data);
        break;
    case OUTLIER:
        pixelDesvFilter(img, factor);
        break;

    }
}

void LogFilter::defineParams(XmippProgram * program)
{
    program->addParamsLine("== Log filter (for scanners) uses equation fa-fb*log(x+fc) ==");
    program->addParamsLine("  [--log]");
    program->addParamsLine("  [--fa <a>]");
    program->addParamsLine("  [--fb <b>]");
    program->addParamsLine("  [--fc <c>]");
}

/** Read from program command line */
void LogFilter::readParams(XmippProgram * program)
{
    a = program->getDoubleParam("--fa");
    b = program->getDoubleParam("--fb");
    c = program->getDoubleParam("--fc");
}

/** Apply the filter to an image or volume*/
void LogFilter::apply(MultidimArray<double> &img)
{
        logFilter(img, a,b,c);
}

void DenoiseTVFilter::defineParams(XmippProgram * program)
{
    program->addParamsLine("== Total Variation Denoising method for micrographs ==");
    program->addParamsLine("  [--denoiseTV]");
    program->addParamsLine("  [--maxIterTV <maxIter>]");
}

/** Read from program command line */
void DenoiseTVFilter::readParams(XmippProgram * program)
{
    maxIter = program->getIntParam("--maxIterTV");
}

/** Apply the filter to an image or volume*/
void DenoiseTVFilter::apply(MultidimArray<double> &img)
{
    denoiseTVFilter(img, maxIter);
}

/** Define the parameters for use inside an Xmipp program */
void BackgroundFilter::defineParams(XmippProgram * program)
{
    program->addParamsLine("== Background removal ==");
    program->addParamsLine(
        "  [ --background <type=plane> ]            : Filters to remove the background.");
    program->addParamsLine("         where <type>  ");
    program->addParamsLine(
        "            plane                          : Remove the plane that best fits the pixels.");
    program->addParamsLine(
        "            rollingball <radius>           : The background is computed as a rolling ball operation.");
    program->addParamsLine("         alias -g; ");
}

/** Read from program command line */
void BackgroundFilter::readParams(XmippProgram * program)
{
    type = PLANE;
    // Check operation to do
    String typeStr = program->getParam("--background");
    if (typeStr == "plane") //Nothing to do, plane by default
        ;
    else if (typeStr == "rollingball")
    {
        type = ROLLINGBALL;
        radius = program->getIntParam("--background", "rollingball");
    }
}

/** Apply the filter to an image or volume*/
void BackgroundFilter::apply(MultidimArray<double> &img)
{
    switch (type)
    {
    case PLANE:
        substractBackgroundPlane(img);
        break;
    case ROLLINGBALL:
        substractBackgroundRollingBall(img, radius);
        break;

    }
}

/** Define the parameters for use inside an Xmipp program */
void MedianFilter::defineParams(XmippProgram * program)
{
    program->addParamsLine("== Median ==");
    program->addParamsLine("  [ --median ]            : Use median filter.");
    program->addParamsLine("         alias -m; ");
}

/** Read from program command line */
void MedianFilter::readParams(XmippProgram * program)
{ //Do nothing by now
}

/** Apply the filter to an image or volume*/
void MedianFilter::apply(MultidimArray<double> &img)
{
    static MultidimArray<double> tmp;
    tmp = img;
    medianFilter3x3(tmp, img);
}

/** Define the parameters for use inside an Xmipp program */
void DiffusionFilter::defineParams(XmippProgram * program)
{
    program->addParamsLine("== Anisotropic diffusion ==");
    program->addParamsLine(
        "  [--diffusion]              : Use anisotropic diffusion filter.");
    program->addParamsLine(
        "  [--shah_iter+ <outer=10> <inner=1> <refinement=1>]  : Diffusion outer, inner and refinement iterations");
    program->addParamsLine("     requires --diffusion;");
    program->addParamsLine(
        "  [--shah_weight+ <w0=0> <w1=50> <w2=50> <w3=0.02>]:Diffusion weights");
    program->addParamsLine(
        "                             :  w0 = data matching ");
    program->addParamsLine(
        "                             :  w1 = 1st derivative smooth ");
    program->addParamsLine(
        "                             :  w2 = edge strength ");
    program->addParamsLine(
        "                             :  w3 = edge smoothness ");
    program->addParamsLine("     requires --diffusion;");
    program->addParamsLine(
        "  [--shah_only_edge+]        : Produce the edge image of the diffusion");
    program->addParamsLine("     requires --diffusion;");
}

/** Read from program command line */
void DiffusionFilter::readParams(XmippProgram * program)
{
    Shah_weight.resizeNoCopy(4);
    Shah_outer = program->getIntParam("--shah_iter", 0);
    Shah_inner = program->getIntParam("--shah_iter", 1);
    Shah_refinement = program->getIntParam("--shah_iter", 2);
    Shah_weight(0) = program->getDoubleParam("--shah_weight", 0);
    Shah_weight(1) = program->getDoubleParam("--shah_weight", 1);
    Shah_weight(2) = program->getDoubleParam("--shah_weight", 2);
    Shah_weight(3) = program->getDoubleParam("--shah_weight", 3);
    Shah_edge = program->checkParam("--shah_only_edge");
}

void DiffusionFilter::show()
{
    std::cout << " Shah difussion\n" << " Outer iterations " << Shah_outer
    << std::endl << " Inner iterations " << Shah_inner << std::endl
    << " Refinement interations " << Shah_refinement << std::endl
    << " Weight " << Shah_weight.transpose() << std::endl;
    if (Shah_edge)
        std::cout << " Generating edge image\n";

}

/** Apply the filter to an image or volume*/
void DiffusionFilter::apply(MultidimArray<double> &img)
{
    MultidimArray<double> surface_strength;
    MultidimArray<double> edge_strength;
    smoothingShah(img, surface_strength, edge_strength, Shah_weight, Shah_outer,
                  Shah_inner, Shah_refinement);
    if (Shah_edge)
        img = edge_strength;
    else
        img = surface_strength;
}

/** Define the parameters for use inside an Xmipp program */
void BasisFilter::defineParams(XmippProgram * program)
{
    program->addParamsLine("== Basis filter ==");
    program->addParamsLine(
        "  [--basis <file> <N=-1>]           : Stack file with the basis, N is the number of elements to consider");
}

/** Read from program command line */
void BasisFilter::readParams(XmippProgram * program)
{
    fnBasis = program->getParam("--basis");
    Nbasis = program->getIntParam("--basis", 1);
    basis.read(fnBasis);
    if (Nbasis > 0)
        basis().resize(Nbasis, 1, YSIZE(basis()), XSIZE(basis()));
}

void BasisFilter::show()
{
    std::cout << " Basis filter\n" << " Basis file " << fnBasis << std::endl
    << " Number of basis " << Nbasis << std::endl;
}

/** Apply the filter to an image or volume*/
void BasisFilter::apply(MultidimArray<double> &img)
{
    const MultidimArray<double> &mBasis = basis();
    if (XSIZE(img) != XSIZE(mBasis) || YSIZE(img) != YSIZE(mBasis))
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     "Images and basis are of different size");

    MultidimArray<double> result;
    result.initZeros(img);
    for (size_t nn = 0; nn < NSIZE(mBasis); ++nn)
    {
        double cnn = 0;
        double *ptrBasis = &NZYX_ELEM(mBasis,nn,0,0,0);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img)
        cnn += DIRECT_MULTIDIM_ELEM(img,n) * (*ptrBasis++);
        ptrBasis = &NZYX_ELEM(mBasis,nn,0,0,0);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img)
        DIRECT_MULTIDIM_ELEM(result,n) += cnn * (*ptrBasis++);
    }
    img = result;
}

/** Define the parameters for use inside an Xmipp program */
void RetinexFilter::defineParams(XmippProgram * program)
{
    program->addParamsLine("== Retinex ==");
    program->addParamsLine(
        "  [ --retinex <percentile=0.9> <mask_file=\"\"> <eps=1>] : Retinex filter, remove a given percentile of the");
    program->addParamsLine("         : Laplacian within the mask, epsilon controls the behaviour at low frequency");
}

/** Read from program command line */
void RetinexFilter::readParams(XmippProgram * program)
{
    percentile = program->getDoubleParam("--retinex",0);
    String fnMask = program->getParam("--retinex",1);
    if (fnMask=="")
    	mask = nullptr;
    else
    {
        mask = new Image<int>;
        mask->read(fnMask);
    }
    eps = program->getDoubleParam("--retinex",2);
}
void RetinexFilter::laplacian(const MultidimArray<double> &img,
   MultidimArray< std::complex<double> > &fimg, bool direct)
{
    if (ZSIZE(img)>1)
    {
        Matrix1D<double> w(3);
        for (size_t k=0; k<ZSIZE(fimg); k++)
        {
            FFT_IDX2DIGFREQ(k,ZSIZE(img),ZZ(w));
            double filter_k=6+eps-2*cos(TWOPI*ZZ(w));
            for (size_t i=0; i<YSIZE(fimg); i++)
            {
                FFT_IDX2DIGFREQ(i,YSIZE(img),YY(w));
                double filter_ki=filter_k-2*cos(TWOPI*YY(w));
                for (size_t j=0; j<XSIZE(fimg); j++)
                {
                    FFT_IDX2DIGFREQ(j,XSIZE(img),XX(w));
                    double filter_kij=filter_ki-2*cos(TWOPI*XX(w));
                    if (!direct && filter_kij>0)
                        filter_kij=1.0/filter_kij;
                    DIRECT_A3D_ELEM(fimg,k,i,j)*=filter_kij;
                }
            }
        }
    }
    else
    {
        Matrix1D<double> w(2);
		for (size_t i=0; i<YSIZE(fimg); i++)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(img),YY(w));
			double filter_i=4+eps-2*cos(TWOPI*YY(w));
			for (size_t j=0; j<XSIZE(fimg); j++)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(img),XX(w));
				double filter_ij=filter_i-2*cos(TWOPI*XX(w));
				if (!direct && filter_ij>0)
					filter_ij=1.0/filter_ij;
				DIRECT_A2D_ELEM(fimg,i,j)*=filter_ij;
			}
		}
    }
}

/** Apply the filter to an image or volume*/
void RetinexFilter::apply(MultidimArray<double> &img)
{
    FourierTransformer fft;
    MultidimArray< std::complex<double> > fimg;

    // Apply forward Laplacian
    fft.FourierTransform(img,fimg,false);
    laplacian(img,fimg,true);
    fft.inverseFourierTransform();

    // Compute the threshold outside the mask
    double *sortedValues=nullptr;
    size_t Nsorted=0;
    if (mask==nullptr)
    {
        Nsorted=MULTIDIM_SIZE(img);
        sortedValues=new double[Nsorted];
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img)
            sortedValues[n]=fabs(DIRECT_MULTIDIM_ELEM(img,n));
    }
    else
    {
        const MultidimArray<int> &mmask=(*mask)();
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img)
            if (DIRECT_MULTIDIM_ELEM(mmask,n)==0)
               Nsorted+=1;

        sortedValues=new double[Nsorted];
        size_t i=0;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img)
            if (DIRECT_MULTIDIM_ELEM(mmask,n)==0)
                 sortedValues[i++]=fabs(DIRECT_MULTIDIM_ELEM(img,n));
    }
    std::sort(sortedValues,sortedValues+Nsorted);
    double threshold=sortedValues[(size_t)(percentile*Nsorted)];
    delete[] sortedValues;

    // Apply threshold
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img)
        if (fabs(DIRECT_MULTIDIM_ELEM(img,n))<threshold)
            DIRECT_MULTIDIM_ELEM(img,n)=0.0;

    fft.FourierTransform(img,fimg,false);
    laplacian(img,fimg,false);
    fft.inverseFourierTransform();
}

// Denoising using Spectral Projected Gradient (SPG) optimization
/** Compute energy.
 *  * @ingroup Filters
 *  energy function
 *  0.5*|vst(X)-Y|^2 + mu*TV(X) = data_term + mu*TV
 *
 * X=current image
 * Y=sensed image
 * mu=regularization weight (on TV norm)
 */
double denoiseTVenergy(double mu,
              const MultidimArray<double>& X,
              const MultidimArray<double>& Y,
              double s,
              double q,
              double lambda,
              double sigmag,
              double g)
{
    // parameter beta which role is to prevent division with zeros
    const double beta = 0.00001;
    const double beta2 = beta*beta;

    double m = 0;
    double tv = 0;
    const double K1 = ((3.0/8.0)*lambda*lambda + sigmag*sigmag - lambda*g) / (s*s);
    const double K2 = lambda * (q/(s*s));
    const double K3 = 2.0 / lambda;
    double Xij;
    double Yij;
    double dXx;
    double dXy;
    double d;
    double msqrt;
    int i1;
    int j1;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(X)
    {
        Xij = DIRECT_A2D_ELEM(X,i,j);
        Yij = DIRECT_A2D_ELEM(Y,i,j);
        i1 = i + 1;
        if (i1==YSIZE(X))
           i1 = 0;
        j1 = j + 1;
        if (j1==XSIZE(X))
           j1 = 0;

        dXx = DIRECT_A2D_ELEM(X,i,j1) - Xij;
        dXy = DIRECT_A2D_ELEM(X,i1,j) - Xij;
        tv += sqrt(dXx*dXx + dXy*dXy + beta2);

        d = K2*Xij + K1;
        msqrt = K3 * sqrt(std::max(0.0, d)) - Yij;
        m += msqrt*msqrt;
    }
    return 0.5*m + mu*tv;
}

/** Compute gradient.
 *  * @ingroup Filters
 *  Jacobian of Energy_TV
 *  0.5*|vst(X)-Y|^2 + mu*TV(X)
 *
 * X=current image
 * Y=sensed image
 * mu=regularization weight (on TV norm)
 */
void denoiseTVgradient(double mu,
                const MultidimArray<double>& X,
                const MultidimArray<double>& Y,
                double s,
                double q,
                double lambda,
                double sigmag,
                double g,
                MultidimArray<double>& gradientI
                )
{
    // parameter beta which role is to prevent division with zeros
    const double beta = 0.00001;
    const double beta2 = beta*beta;

    MultidimArray<double> d = X;
    double K1 = ((3.0/8.0)*lambda*lambda + sigmag*sigmag - lambda*g) / (s*s);
    double K2 = lambda * (q/s*s);
    double K3 = (2.0 / (lambda*lambda)) * (q / (s*s)) * lambda;
    double Xij;
    double Yij;
    double dXx;
    double dXy;
    double dij;
    double d_left;
    double d_up;
    double X_left;
    double X_right;
    double X_up;
    double X_down;
    double dTV;
    double dE;
    int i1;
    int j1;
    int i_1;
    int j_1;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(X)
    {
        Xij = DIRECT_A2D_ELEM(X,i,j);
        i1 = i + 1;
        if (i1==YSIZE(X))
           i1 = 0;
        j1 = j + 1;
        if (j1==XSIZE(X))
           j1 = 0;

        dXx = DIRECT_A2D_ELEM(X,i,j1) - Xij;
        dXy = DIRECT_A2D_ELEM(X,i1,j) - Xij;
        DIRECT_A2D_ELEM(d,i,j) = 1.0 / sqrt(dXx*dXx + dXy*dXy + beta2);
    }

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(d)
    {
        dij = DIRECT_A2D_ELEM(d,i,j);
        Xij = DIRECT_A2D_ELEM(X,i,j);
        Yij = DIRECT_A2D_ELEM(Y,i,j);
        i1 = i + 1;
        if (i1==YSIZE(d))
           i1 = 0;
        j1 = j + 1;
        if (j1==XSIZE(d))
           j1 = 0;
        i_1 = i - 1;
        if (i_1==-1)
           i_1 = YSIZE(d) - 1;
        j_1 = j - 1;
        if (j_1==-1)
           j_1 = XSIZE(d) - 1;

        d_left = DIRECT_A2D_ELEM(d,i,j_1);
        d_up = DIRECT_A2D_ELEM(d,i_1,j);
        X_left = DIRECT_A2D_ELEM(X,i,j_1);
        X_right = DIRECT_A2D_ELEM(X,i,j1);
        X_up = DIRECT_A2D_ELEM(X,i_1,j);
        X_down = DIRECT_A2D_ELEM(X,i1,j);

        dTV = Xij * (2.0*dij + d_left + d_up) - X_left*d_left - X_up*d_up - dij*(X_right + X_down);

        if (K2*Xij + K1 > 0)
            dE = K3 - (q / (s*s)) * Yij / sqrt(Xij * (q / (s*s)) * lambda + K1);
        else
            dE = 0;

        DIRECT_A2D_ELEM(gradientI,i,j) = dE + mu*dTV;
    }
}

/** Function which makes projections on feasible set.
 *  * @ingroup Filters
 */
void denoiseTVproj(const MultidimArray<double>& X,
          const MultidimArray<double>& Y,
          double theta,
          MultidimArray<double>& dold)
{
    double div;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(X)
    {
        div = DIRECT_A2D_ELEM(X,i,j) - (DIRECT_A2D_ELEM(Y,i,j)*theta);
        if (div < 0)
           DIRECT_A2D_ELEM(dold,i,j) = 0 - DIRECT_A2D_ELEM(X,i,j);
        else if (div > 1)
           DIRECT_A2D_ELEM(dold,i,j) = 1 - DIRECT_A2D_ELEM(X,i,j);
        else
           DIRECT_A2D_ELEM(dold,i,j) = div - DIRECT_A2D_ELEM(X,i,j);
    }
}

/** Compute denoising.
 *  * @ingroup Filters
 */
void denoiseTVFilter(MultidimArray<double> &xnew, int maxIter)
{
    // parameters for denoising
    double lambda = 1.0;    // Poisson scaling factor
    double sigmag = 5.8;  // Gaussian component N(g,sigma^2)
    double g = 0.0;

    // *** Anscombe transformation of input image ***
    double q = 255.0;

    double xnewmin = xnew.computeMin();
    double xnewscale = 255.0 / (xnew.computeMax() - xnewmin);

    // apply generalized Anscombe VST tranformation
    double K1 = ((3.0/8.0)*lambda*lambda + sigmag*sigmag - lambda*g);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(xnew)
    {
        DIRECT_MULTIDIM_ELEM(xnew, n) = (DIRECT_MULTIDIM_ELEM(xnew, n) - xnewmin) * xnewscale;
        DIRECT_MULTIDIM_ELEM(xnew, n) = 2.0 / lambda * sqrt(std::max(0.0, lambda*DIRECT_MULTIDIM_ELEM(xnew, n) + K1));
    }

    double mx = xnew.computeMax();

    // xold = initial degraded image
    MultidimArray<double> xold = xnew;
    MultidimArray<double> grold = xnew;
    MultidimArray<double> grnew = xnew;
    MultidimArray<double> dold = xnew;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xnew)
    {
        DIRECT_A2D_ELEM(xold,i,j) = DIRECT_A2D_ELEM(xnew,i,j) / mx;
    }
    MultidimArray<double> origInput = xold;

    // *** Spectral Projected Gradient (SPG) optimization ***
    // parameters of optimization
    double tol = 0;
    double theta;
    double thetamin = 0.001;
    double thetamax = 1000;
    double gamma = 0.0001;
    double sigma1 = 0.1;
    double sigma2 = 0.9;
    double mu = 0.03;

    // initial objective function value of xold
    double fold = denoiseTVenergy(mu, xold, origInput, mx, q, lambda, sigmag, g);

    // gradient of objective function
    denoiseTVgradient(mu, xold, origInput, mx, q, lambda, sigmag, g, grold);

    denoiseTVproj(xold, grold, 1.0, dold);

    double delta;
    double ksi;
    double fnew;
    double s_norm;
    double p;
    double s2;
    double xij;
    for (int kk=1; kk <= maxIter; kk++)
    {
        delta = 0;

        // make a step
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xnew)
        {
            DIRECT_A2D_ELEM(xnew,i,j) = DIRECT_A2D_ELEM(xold,i,j) + DIRECT_A2D_ELEM(dold,i,j);
            delta += DIRECT_A2D_ELEM(grold,i,j) * DIRECT_A2D_ELEM(dold,i,j);
        }

        // objective function of xnew
        fnew = denoiseTVenergy(mu, xnew, origInput, mx, q, lambda, sigmag, g);

        ksi = 1.0;

        while (fnew > fold + gamma*ksi*delta)
        {
            double ksitsl = -0.5 * (ksi*ksi) * delta / (fnew - fold - ksi*delta);
            if ((ksitsl >= sigma1)  && (ksitsl <= sigma2*ksi))
                ksi = ksitsl;
            else
                ksi = ksi / 2.0;

            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xnew)
            {
                DIRECT_A2D_ELEM(xnew,i,j) = DIRECT_A2D_ELEM(xold,i,j) + ksi*DIRECT_A2D_ELEM(dold,i,j);
            }

            // because we updated xnew, update also function value
            fnew = denoiseTVenergy(mu, xnew, origInput, mx, q, lambda, sigmag, g);
        }

        denoiseTVgradient(mu, xnew, origInput, mx, q, lambda, sigmag, g, grnew);

        s_norm = -100;
        p = 0;
        s2 = 0;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xnew)
        {
            xij = (DIRECT_A2D_ELEM(xnew,i,j) - DIRECT_A2D_ELEM(xold,i,j));
            p += xij * (DIRECT_A2D_ELEM(grnew,i,j) - DIRECT_A2D_ELEM(grold,i,j));
            s2 += xij * xij;
            if (xij > s_norm) s_norm = xij;
        }

        if (s_norm < tol)
            break;

        if (p <= 0)
            theta = thetamax;
        else
            theta = std::min(thetamax, std::max(thetamin, s2/p));

        denoiseTVproj(xnew, grnew, theta, dold);

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xold)
        {
            DIRECT_A2D_ELEM(xold,i,j) = DIRECT_A2D_ELEM(xnew,i,j);
            DIRECT_A2D_ELEM(grold,i,j) = DIRECT_A2D_ELEM(grnew,i,j);
        }

        // we already computed it, let's use it for next round
        fold = fnew;
    }
}

void matlab_filter(Matrix1D<double> &B, Matrix1D<double> &A,
		    const MultidimArray<double> &X, MultidimArray<double> &Y, MultidimArray<double> &Z)
{
	    if (VEC_ELEM(A,0)!=1.0)
	    {
	    	A/=VEC_ELEM(A,0);
	    	B/=VEC_ELEM(B,0);
	    }

		size_t input_size = XSIZE(X);
		size_t filter_order = std::max(VEC_XSIZE(A),VEC_XSIZE(B));
		if (VEC_XSIZE(B)!=filter_order)
			B.resize(filter_order);
		if (VEC_XSIZE(A)!=filter_order)
			A.resize(filter_order);
		Y.resize(input_size);
		Z.initZeros(filter_order);

		for (size_t i = 0; i < input_size; ++i)
		{
			size_t order = filter_order - 1;
			while (order)
			{
				if (i >= order)
				{
					DIRECT_A1D_ELEM(Z,order - 1) = VEC_ELEM(B,order) * DIRECT_A1D_ELEM(X,i - order) -
							                       VEC_ELEM(A,order) * DIRECT_A1D_ELEM(Y,i - order) + DIRECT_A1D_ELEM(Z,order);
				}
				--order;
			}
			DIRECT_A1D_ELEM(Y,i) = VEC_ELEM(B,0) * DIRECT_A1D_ELEM(X,i) + DIRECT_A1D_ELEM(Z,0);
		}
}

template <typename T>
double mutualInformation(const MultidimArray< T >& x,
                         const MultidimArray< T >& y,
                         int nx,
                         int ny,
                         const MultidimArray< int >* mask)
{
    SPEED_UP_temps;

    long n = 0;
    Histogram1D histx;
    Histogram1D histy;
    Histogram2D histxy;
    MultidimArray< T > aux_x;
    MultidimArray< T > aux_y;
    MultidimArray< double > mx;
    MultidimArray< double > my;
    MultidimArray< double > mxy;
    int xdim;
    int ydim;
    int zdim;
    double retval = 0.0;

    xdim=XSIZE(x);
    ydim=YSIZE(x);
    zdim=ZSIZE(x);
    aux_x.resize(xdim * ydim * zdim);
    xdim=XSIZE(y);
    ydim=YSIZE(y);
    zdim=ZSIZE(y);
    aux_y.resize(xdim * ydim * zdim);

    FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(x, y)
    {
        if (mask != nullptr)
            if (!(*mask)(k, i, j))
                continue;

        aux_x(n) = A3D_ELEM(x, k, i, j);
        aux_y(n) = A3D_ELEM(y, k, i, j);
        n++;
    }

    aux_x.resize(n);
    aux_y.resize(n);

    if (n != 0)
    {
        if (nx == 0)
            //Assume Gaussian distribution
            nx = (int)((log((double) n) / LOG2) + 1);

        if (ny == 0)
            //Assume Gaussian distribution
            ny = (int)((log((double) n) / LOG2) + 1);

        compute_hist(aux_x, histx, nx);
        compute_hist(aux_y, histy, ny);
        compute_hist(aux_x, aux_y, histxy, nx, ny);

        mx = histx;
        my = histy;
        mxy = histxy;
        for (int i = 0; i < nx; i++)
        {
            double histxi = (histx(i)) / n;
            for (int j = 0; j < ny; j++)
            {
                double histyj = (histy(j)) / n;
                double histxyij = (histxy(i, j)) / n;
                if (histxyij > 0)
                    retval += histxyij * log(histxyij / (histxi * histyj)) /
                              LOG2;
            }
        }

        return retval;
    }
    else
        return 0;
}
