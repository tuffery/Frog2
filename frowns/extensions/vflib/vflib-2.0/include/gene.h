/*----------------------------------------------------
 * gene.h
 * Interface of gene.cc
 * Random generation of isomorphic ARGraphs
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: gene.h,v 1.1 2003/03/19 20:58:52 wc2so1 Exp $
 ----------------------------------------------------*/

/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: gene.h,v $
 *   Revision 1.1  2003/03/19 20:58:52  wc2so1
 *   Added extensions directory for extensions to frowns (they can now
 *   all be built together)
 *
 *   Seperated MoleculeDrawer and TkMoleculeDrawer
 *
 *   Updated perception code and removed a couple of bugs.
 *
 *   Revision 1.2  1998/12/08 13:31:01  foggia
 *   Minor changes
 *
 *---------------------------------------------------------------*/

#ifndef GENE_H

#include "argraph.h"



void Generate(int nodes, int edges, Graph **g1, Graph **g2, 
              bool connected=true);

#endif
