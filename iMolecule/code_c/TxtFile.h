/* ========================================================================== */
/*                                                                            */
/*  LigCSRre (C) 2008 - Small compounds 3D molecular similarity screening by: */
/*                                                                            */
/* Pierre Tuffery(1), Flavien Quintus(1), J. Grynberg(1), O. Sperandio(1),    */
/* M. Petitjean(2)                                                            */
/*  (1): MTi, INSERM UMR-S 973, Université Paris Diderot - Paris 7,           */
/*     	 F75013, Paris, France.                                               */
/* 	 (http://www.mti.univ-paris-diderot.fr)                               */
/* (2): CEA/DSV/iBiTec-S/SB2SM (CNRS URA 2096),                               */
/*       F91191, Gif-sur-Yvette, France.                                      */
/*                                                                            */
/* Grants to use: See the LICENSE.TXT file that comes with the package.       */
/*                                                                            */
/* ========================================================================== */
#ifndef __TXTFILE_H__
#define __TXTFILE_H__


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

extern char ** txtFileRead(char *fname, int *lSze);
extern char ** txtFileFree(char **l, int lSze);
extern char ** stdinRead(int *lSze);

#endif
