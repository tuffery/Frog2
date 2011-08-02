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

       #include <stdio.h>
       #include <stdlib.h>
       #include <string.h>

       int
       main(int argc, char *argv[])
       {
           char *str1, *str2, *token, *subtoken;
           char *saveptr1, *saveptr2;
           int j;

           if (argc != 4) {
               fprintf(stderr, "Usage: %s string delim subdelim\n",
                       argv[0]);
               exit(EXIT_FAILURE);
           }

           for (j = 1, str1 = argv[1]; ; j++, str1 = NULL) {
               token = strtok_r(str1, argv[2], &saveptr1);
               if (token == NULL)
                   break;
               printf("%d: %s", j, token);

               for (str2 = token; ; str2 = NULL) {
                   subtoken = strtok_r(str2, argv[3], &saveptr2);
                   if (subtoken == NULL)
                       break;
                   printf(" --> %s0, subtoken);
               }
           }

           exit(EXIT_SUCCESS);
       } /* main */

