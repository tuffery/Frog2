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
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>


/* ===========================================================
 * Renvoie les lignes d'un fichier et leur nombre de lignes.
 * ===========================================================
 */

char ** txtFileRead(char *fname, int *lSze)
{
  struct stat buf;
  char *tmp;
  char **line;
  int fd;
  int fSze, done;
  int nLines;
  int i,j;

  if ((fd = open(fname,O_RDONLY)) < 0) {
    return NULL;
  }

  if (fstat(fd,&buf)) {
    close(fd);
    return NULL;
  }

  if (!S_ISREG (buf.st_mode)) {
    close(fd);
    return NULL;
  }

  if (!(fSze = buf.st_size)) {
    close(fd);
    return NULL;
  }

  /* La taille du buffer ne dépasse-t-elle pas les capacités de fSze ? */
  if (fSze != buf.st_size) {
    fprintf(stderr, "%s: File too heavy to be loaded.\n", fname);
    close(fd);
    return NULL;
  }

  if (!(tmp = malloc((fSze+1) * sizeof(char)))) {
    close(fd);
    return NULL;
  }

  done = 0;
  while (done < fSze) {
    int r = read(fd, &(tmp[done]), fSze - done);
    if (r < 0) {
      free(tmp);
      close(fd);
      return NULL;
    }
    done += r;
  }
  close(fd);

  nLines = 0;
  for (i = 0; i < fSze; i++) {
    if (tmp[i] == '\n') nLines++;
  }
  if (tmp[fSze-1] != '\n') nLines++;

  if (!(line = malloc((nLines+2) * sizeof(char *)))) {
      free(tmp);
      return NULL;
  }

  line[0] = tmp;
  j = 1;
  for (i = 0; i < fSze; i++) {
    if (tmp[i] == '\n') {
      line[j++] = &tmp[i+1];
      tmp[i] = '\000';
    }
  }

  line[j] = NULL;
  *lSze = nLines;

  return (line);
}

/* ===============================================
 * Réalloue la mémoire occupée par txtFileRead.
 * ===============================================
 */

char ** txtFileFree(char **l, int lSze)
{
  /* Attention: on ne libere que la premiere ligne,
     car c'est un seul malloc pour l'ensemble des lignes
  */
  free(l[0]);
  free(l);
  return (NULL);
}

/* ===============================================
 * Grab from stdin
 * ===============================================
 */
char ** stdinRead(int *lSze)
{
  char  buff[BUFSIZ];
  char **line;
  char *tmp = NULL;
  int   nLines;
  int   fSze = 0;
  FILE *fd = stdin;
  int   i, j;

  tmp = malloc(1000 * sizeof(char));
  tmp[0] = '\0';
  while ((fgets(buff, BUFSIZ, stdin)) != NULL) {
    fSze += strlen(buff);
    if (!(tmp = realloc(tmp, (fSze+1) * sizeof(char)))) {
      return NULL;
    }
    strcat(tmp, buff);
  }

  // Number of lines
  nLines = 0;
  for (i = 0; i < fSze; i++) {
    if (tmp[i] == '\n') nLines++;
  }
  if (tmp[fSze-1] != '\n') nLines++;

  // Array of lines
  if (!(line = malloc((nLines+2) * sizeof(char *)))) {
      free(tmp);
      return NULL;
  }

  // the lines
  line[0] = tmp;
  j = 1;
  for (i = 0; i < fSze; i++) {
    if (tmp[i] == '\n') {
      line[j++] = &tmp[i+1];
      tmp[i] = '\000';
    }
  }

  line[j] = NULL;
  *lSze = nLines;

  return (line);

}
