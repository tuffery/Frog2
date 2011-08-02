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

#define DtFloat         double          /* Le type flottant.            */
typedef DtFloat         DtPoint3[3];
typedef DtFloat         DtPoint4[4];
typedef DtFloat         DtMatrix3x3[3][3];
typedef DtFloat         DtMatrix4x4[4][4];


#define boucle(i,min,max) for((i)=(min);(i)<(max);(i)++)
#define max(a,b) (a<b?b:a)
#define min(a,b) (a<b?a:b)
#define sign(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))
#define EPS 1.e-10


/* ===============================================================
 * ---------- Transformation d'une coordonnee par rmat. ----------
 * =============================================================== */

void crdRotate(DtPoint3 p,DtMatrix4x4 rmat)
{
  double x,y,z;

  x = p[0] * rmat[0][0] + p[1] * rmat[1][0] + p[2] * rmat[2][0] + rmat[3][0];
  y = p[0] * rmat[0][1] + p[1] * rmat[1][1] + p[2] * rmat[2][1] + rmat[3][1];
  z = p[0] * rmat[0][2] + p[1] * rmat[1][2] + p[2] * rmat[2][2] + rmat[3][2];
  p[0] = x;
  p[1] = y;
  p[2] = z;
}

int main(int argc, char *argv[])
{
  DtPoint3 *c1, *c2;
  DtMatrix4x4 M;
  double squared_rmsd;
  int i,j;
  int n;

  fscanf(stdin,"%d",&n);
  c1 = calloc(n,sizeof(DtPoint3));

  for (i=0; i<n; i++) {
    j = i*3 + 1;
    fscanf(stdin,"%lf%lf%lf",&c1[i][0],&c1[i][1],&c1[i][2]);
  }
  fscanf(stdin,"%lf%lf%lf%lf",&M[0][0],&M[0][1],&M[0][2],&M[0][3]);
  fscanf(stdin,"%lf%lf%lf%lf",&M[1][0],&M[1][1],&M[1][2],&M[1][3]);
  fscanf(stdin,"%lf%lf%lf%lf",&M[2][0],&M[2][1],&M[2][2],&M[2][3]);
  fscanf(stdin,"%lf%lf%lf%lf",&M[3][0],&M[3][1],&M[3][2],&M[3][3]);

  for (i=0; i<n; i++) {
    crdRotate(c1[i],M);
    fprintf(stdout,"%10.3lf%10.3lf%10.3lf\n",c1[i][0],c1[i][1],c1[i][2]);
  }

/*   for (i=0; i<n; i++) { */
/*     fprintf(stdout,"ATOM  %5d  %.3s%c%.3s %c%4d%c   %8.3lf%8.3lf%8.3lf\n", */
/*          i, "CA", ' ', "ALA", ' ',i,' ',c2[i][0],c2[i][1],c2[i][2]); */
/*   } */

}
