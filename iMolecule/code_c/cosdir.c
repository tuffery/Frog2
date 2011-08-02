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

/* ===============================
 * Matrice de passage ???
 * ===============================
 */
static void cosdir(DtMacRec *M, 
	    int i, int j, int k, 
	    DtPoint3 A, DtPoint3 B, DtPoint3 C)
{
  /* Local variables */
  static DtFloat a, b, c, pm;

  /* Vecteur ij norme */
  A[0] = M->corm[j][0] - M->corm[i][0];
  B[0] = M->corm[j][1] - M->corm[i][1];
  C[0] = M->corm[j][2] - M->corm[i][2];
  pm = sqrt(A[0] * A[0] + B[0] * B[0] + C[0] * C[0]);
  A[0] /= pm;
  B[0] /= pm;
  C[0] /= pm;

  /* Vecteur jk */
  a = M->corm[k][0] - M->corm[j][0];
  b = M->corm[k][1] - M->corm[j][1];
  c = M->corm[k][2] - M->corm[j][2];

  /* Vecteur ?? norme */
  A[2] = B[0] * c - C[0] * b;
  B[2] = C[0] * a - A[0] * c;
  C[2] = A[0] * b - B[0] * a;
  pm = sqrt(A[2] * A[2] + B[2] * B[2] + C[2] * C[2]);
  A[2] /= pm;
  B[2] /= pm;
  C[2] /= pm;

  /* Vecteur ?? norme */
  A[1] = B[2] * C[0] - C[2] * B[0];
  B[1] = C[2] * A[0] - A[2] * C[0];
  C[1] = A[2] * B[0] - B[2] * A[0];
  pm = sqrt(A[1] * A[1] + B[1] * B[1] + C[1] * C[1]);
  A[1] /= pm;
  B[1] /= pm;
  C[1] /= pm;
}
