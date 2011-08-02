/*----------------------------------------------------
 * gene_mesh.h
 * Interface of gene_mesh.cc
 * Random generation of isomorphic ARGraphs which are
 * almost meshes.
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: gene_mesh.h,v 1.1 2003/03/19 20:58:52 wc2so1 Exp $
 ----------------------------------------------------*/

/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: gene_mesh.h,v $
 *   Revision 1.1  2003/03/19 20:58:52  wc2so1
 *   Added extensions directory for extensions to frowns (they can now
 *   all be built together)
 *
 *   Seperated MoleculeDrawer and TkMoleculeDrawer
 *
 *   Updated perception code and removed a couple of bugs.
 *
 *   Revision 1.2  1999/01/06 20:03:53  foggia
 *   Added the option to generate subgraphs
 *
 *   Revision 1.1  1998/12/26 17:28:14  foggia
 *   Initial revision
 *
 *   Revision 1.2  1998/12/08 13:31:01  foggia
 *   Minor changes
 *
 *---------------------------------------------------------------*/

#ifndef GENE_MESH_H

#include "argraph.h"



void GenerateMesh(int nodes, int extra_edges, Graph **g1, Graph **g2,
                  int sub_nodes=-1);

#endif
