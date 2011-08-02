/* ============================================================================
 *  This software is part of Frog, a chemo informatics class able to build 
 *  3D coordinates for small compounds
 *   Copyright (C) 2006-2007 P. Tuffery, B.O. Villoutreix, Th. Bohme Leite, D. Gomes, M. Miteva, J. Chomilier
 *
 *   Frog2 (C) 2009-2010 by P. Tuffery, M. Miteva, F. Guyon
 *
 *   Using this software, please cite:
 *       Frog2: Efficient 3D conformation ensemble generator for small compounds.
 *       Miteva MA, Guyon F, Tuffery P.
 *       Nucleic Acids Res. 2010 Jul;38(Web Server issue):W622-7. Epub 2010 May 5.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *   ==========================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DtFloat double
typedef DtFloat DtMatrix4x4[4][4];
typedef DtFloat DtPoint3[3];
/* ===============================================================
 * ---------- Matrice de rotation autour d'un axe a b. -----------
 * ---------- PREMULTIPLIE -> Y = XM (X vecteur ligne) -----------
 * =============================================================== */
void MkArbitraryAxisRotMat4x4(DtPoint3 A,DtPoint3 B,DtFloat angle,DtMatrix4x4 M)
{
  extern double sqrt(), sin(), cos();
  
  DtFloat rx, ry, rz, lnorme, rx2, rxry, rxrz, ry2, ryrz, rz2,rxry1mcosa, rxrz1mcosa, ryrz1mcosa;

  DtFloat cosa, sina;

  cosa = cos((double) angle);
  sina = sin((double) angle);

  rx = B[0] - A[0];
  ry = B[1] - A[1];
  rz = B[2] - A[2];
  lnorme = (DtFloat) sqrt(rx*rx + ry*ry + rz*rz);
  rx /= lnorme;
  ry /= lnorme;
  rz /= lnorme;

  rx2 = rx*rx;
  ry2 = ry*ry;
  rz2 = rz*rz;
  rxry = rx*ry;
  rxrz = rx*rz;
  ryrz = ry*rz;
  rxry1mcosa = rxry * (1. - cosa);
  rxrz1mcosa = rxrz * (1. - cosa);
  ryrz1mcosa = ryrz * (1. - cosa);

  M[0][0] = rx2 + (1. - rx2) * cosa;
  M[1][0] = rxry1mcosa - rz * sina;
  M[2][0] = rxrz1mcosa + ry * sina;
  
  M[0][1] = rxry1mcosa + rz * sina;
  M[1][1] = ry2 + (1. - ry2) * cosa;
  M[2][1] = ryrz1mcosa - rx * sina;

  M[0][2] = rxrz1mcosa - ry * sina;
  M[1][2] = ryrz1mcosa + rx * sina;
  M[2][2] = rz2 + (1. - rz2) * cosa;

  M[3][0] =  A[0] * (1 - M[0][0]) - A[1] * M[1][0] - A[2] * M[2][0];
  M[3][1] = -A[0] * M[0][1] + A[1] * (1 - M[1][1]) - A[2] * M[2][1];
  M[3][2] = -A[0] * M[0][2] - A[1] * M[1][2] + A[2] * (1 - M[2][2]);

  M[0][3] = M[1][3] = M[2][3] = 0.;
  M[3][3] = 1.;
}

int main(int argc, char *argv[]) {
  DtMatrix4x4 M;
  DtFloat angle;
  DtPoint3 c1, c2;

  fscanf(stdin,"%lf%lf%lf",&c1[0],&c1[1],&c1[2]);
  fscanf(stdin,"%lf%lf%lf",&c2[0],&c2[1],&c2[2]);
  fscanf(stdin,"%lf",&angle);

  MkArbitraryAxisRotMat4x4(c1,c2,angle,M);

  fprintf(stdout,"%lf %lf %lf %lf\n",M[0][0], M[0][1], M[0][2], M[0][3]);
  fprintf(stdout,"%lf %lf %lf %lf\n",M[1][0], M[1][1], M[1][2], M[1][3]);
  fprintf(stdout,"%lf %lf %lf %lf\n",M[2][0], M[2][1], M[2][2], M[2][3]);
  fprintf(stdout,"%lf %lf %lf %lf\n",M[3][0], M[3][1], M[3][2], M[3][3]);
}
  
