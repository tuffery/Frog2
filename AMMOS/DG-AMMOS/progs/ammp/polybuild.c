/* polybuild.c
*
* write a united atom /residue ammp script for building 
* a polymer.
* read a three letter code list of residues
*
*  tsearch residue at a time until built
*  optimize against bond angle torsion and tether
*/
#define MAXATOM 3000
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
void main()
{

    int ires[MAXATOM];
    char aid[MAXATOM][10];
    FILE *out;
    int iatom,jj;
    char res[4];

    iatom = 0;
    while ( scanf("%s",&res[0]) != EOF )
    {
        ires[iatom] = iatom*100 + 100; /* q+d */
        for( jj=0; jj< 3; jj++)
            aid[iatom][jj] = res[jj];
        aid[iatom][3] = '.';
        aid[iatom][4] = 'c';
        aid[iatom][5] = 'a';
        aid[iatom][6] = '\0';
        iatom ++;
    }

    /* default for q+d */
    out = stdout;
    /* read the stuff in */

    /*  now we should write the file with the torsions */
    fprintf(out," echo off; setf qu 0.;\n setf au 20.;\n setf bu 6000.; \n setf mu 40.; \n");
    fprintf(out," setf kbu 100.;\n setf kau 100.;\n setf ang 109.5;\n setf tku 10.;\n");
    fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
            0.,0.,0.,ires[0],&aid[0][0]);
    fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
            0.,0.,0.,ires[1],&aid[1][0]);
    fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
            0.,0.,0.,ires[2],&aid[2][0]);
    fprintf( out, "bond  %d %d 3.8 kbu;\n",ires[0],ires[1]);
    fprintf(out, "bond  %d %d 3.8 kbu;\n",ires[1],ires[2]);
    fprintf(out,"angle %d %d %d kau 109.5;\n",
            ires[0],ires[1],ires[2]);
    for( jj = 0; jj < iatom -3 ; jj++)
    {
        fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
                0.,0.,1000. + 10.*jj,ires[jj+3],&aid[jj+3][0]);
        fprintf( out," tailor exclude %d %d; \n",
                 ires[jj+2],ires[jj+3]);
    }
    fprintf(out," # tethers go here ; \n");
    fprintf(out," use none tether; cngdel 1000 1000 0; echo; \n");
    fprintf( out," use none bond angle torsion tether;\n");
    for( jj = 0; jj < iatom -3 ; jj++)
    {
        fprintf( out," mov x %d.x; add x 3.8;\n", ires[jj+2]);
        fprintf( out," mov y %d.y; ", ires[jj+2]);
        fprintf( out," mov z %d.z; ", ires[jj+2]);
        fprintf( out,"atom x y z %d %s qu au bu mu;\n",
                 ires[jj+3],&aid[jj+3][0]);

        /*
        	fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
        		0.,0.,0.,ires[jj+3],&aid[jj+3][0]);
        */
        fprintf(out," use none tether; cngdel 10 10 0;\n");
        fprintf(out, "bond  %d %d 3.8 kbu;\n",ires[jj+2],ires[jj+3]);
        /*	fprintf(out,"angle %d %d %d kau 109.5;\n", */
        fprintf(out,"angle %d %d %d kau 100.;\n",
                ires[jj+1],ires[jj+2],ires[jj+3]);
        fprintf(out,"torsion %d %d %d %d  47.81527 1 286.7432;\n"
                ,ires[jj],ires[jj+1],ires[jj+2],ires[jj+3]);
        fprintf(out,"torsion %d %d %d %d  22.79931 2 413.2828;\n"
                ,ires[jj],ires[jj+1],ires[jj+2],ires[jj+3]);
        fprintf(out,"torsion %d %d %d %d  41.48926 3 346.9368;\n"
                ,ires[jj],ires[jj+1],ires[jj+2],ires[jj+3]);
        fprintf( out," use none bond angle  tether;\n");
        fprintf( out," cngdel 20 20 0; v_maxwell 300; pac 20 .00001;\n");
        fprintf( out," cngdel 30 30 0;\n");
        fprintf(out," tgroup 1 %d %d %d %d ; tsearch 1;\n",
                ires[jj],ires[jj+1],ires[jj+2],ires[jj+3]);
        fprintf( out," use none bond angle nonbon torsion tether;\n");
        fprintf( out," cngdel 30 30 0;\n");
    }
    exit(0);
}

