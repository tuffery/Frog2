/* this program reads
* an AMMP file
* it extracts the ca positions and numbers
* and then generates a linked ca ammp file
* the purpose of this is to generate alpha
* carbons for missing atoms based on a group
* search
*/
#define MAXATOM 3000

#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

void main()
{
    char file1[80],file2[80],work[1000];
    FILE *in1,*out,*fopen();
    int i,ifile,ii,jj,kk,iatom;
    float x[MAXATOM];
    float y[MAXATOM];
    float z[MAXATOM];
    char aid[MAXATOM][10];
    int  ires[MAXATOM];
    char aname[10],atype[10],*wp,*cp;
    float a;

    /* strip out the file name and put it in the right place for opening */
    fgets( work,80,stdin );
    for(i= 0; i<80; i++)
        if( work[i] != ' ') break;
    for( ifile = i; ifile < 80 ; ifile++)
    {
        if( work[ifile] == ' ' ) break;
        if( work[ifile] == '\0' ) break;
        if( work[ifile] == '\n' ) break;
        file1[ifile -i ] = work[ifile];
        file1[ifile -i +1 ] = '\0';
    }
    fgets( work,80,stdin );
    for(i= 0; i<80; i++)
        if( work[i] != ' ') break;
    for( ifile = i; ifile < 80 ; ifile++)
    {
        if( work[ifile] == ' ' ) break;
        if( work[ifile] == '\0' ) break;
        if( work[ifile] == '\n' ) break;
        file2[ifile -i ] = work[ifile];
        file2[ifile -i +1 ] = '\0';
    }
    in1 = fopen( file1,"r" );
    out = fopen( file2,"w" );
if( out == NULL) { fprintf(stderr,"cannot open %s\n",file2); exit(0);}
    if( in1 == NULL) { fprintf(stderr,"cannot open %s\n",file1); exit(0);}


    /* now lets read the ammp file and for each atom record check for xxx.ca
      */
    ii = 0;
    iatom = 0;
    while( (work[ii] = fgetc( in1 )) != EOF && ii < 1000 )
    {
        if( isupper(work[ii]))
            work[ii] = tolower(work[ii]);
        /* strip leading blanks and stuff */
        if( ii == 0 && (work[ii] == ' '||
                        work[ii] == '\n' || work[ii] == '\0')) ii = -1;
        /* skip linefeeds and eolines */
        if( work[ii] == '\n' || work[ii] == '\0') ii = ii-1;

        if( work[ii] == ';')
        {
            if( work[0] == 'a' && work[1] == 't' && work[2] == 'o' &&
                    work[3] == 'm')
            { /*then its an atom command */

                wp = &work[4];
                /* skip x,y,z */
                while( *wp == ' ') {wp++; if( *wp == ';'){ wp--; break;} }
                while( *wp != ' ') {wp++; if( *wp == ';'){ wp--; break;} }
                while( *wp == ' ') {wp++; if( *wp == ';'){ wp--; break;} }
                while( *wp != ' ') {wp++; if( *wp == ';'){ wp--; break;} }
                while( *wp == ' ') {wp++; if( *wp == ';'){ wp--; break;} }
                while( *wp != ' ') {wp++; if( *wp == ';'){ wp--; break;} }
                sscanf(wp,"%d %s", &jj,aname);
                cp = &aname[0];
                while( *cp != '.' && *cp != '\0') cp++;
                if( *cp != '\0') cp++;
                if( strcmp("ca",cp) == 0)
                { /* we've got a ca if here */
                    wp = &work[4];
                    sscanf( wp,"%f %f %f",&x[iatom],&y[iatom],&z[iatom]);
                    ires[iatom] = jj;
                    for( jj=0; jj< 10; jj++)
                        aid[iatom][jj] = aname[jj];
                    iatom += 1;
                }
            }
            ii = -1; /* put the line (after increment) in the right place */
        }/* end of end of command if */
        ii++; /* increment ii for reading the next character */
    }/* end of while */
    /*  now we should write the file with the torsions */
    fprintf(out," setf qu 0.;\n setf au 20.;\n setf bu 6000.; \n setf mu 40.; \n");
    fprintf(out," setf kbu 100.;\n setf kau 100.;\n setf ang 109.5;\n setf tku 10.;\n");
    fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
            x[0],y[0],z[0],ires[0],&aid[0][0]);
    fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
            x[1],y[1],z[1],ires[1],&aid[1][0]);
    fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
            x[2],y[3],z[2],ires[2],&aid[2][0]);
    fprintf( out, "bond  %d %d 3.8 kbu;\n",ires[0],ires[1]);
    fprintf(out, "bond  %d %d 3.8 kbu;\n",ires[1],ires[2]);
    fprintf(out,"angle %d %d %d kau 109.5;\n",
            ires[0],ires[1],ires[2]);
    for( jj = 0; jj < iatom -3 ; jj++)
    {
        fprintf(out," atom %f %f %f %d %s qu au bu mu;\n",
                x[jj+3],y[jj+3],z[jj+3],ires[jj+3],&aid[jj+3][0]);
        fprintf(out, "bond  %d %d 3.8 kbu;\n",ires[jj+2],ires[jj+3]);
        fprintf(out,"angle %d %d %d kau 109.5;\n",
                ires[jj+1],ires[jj+2],ires[jj+3]);
        fprintf(out,"torsion %d %d %d %d  0.4781527 1 286.7432;\n"
                ,ires[jj],ires[jj+1],ires[jj+2],ires[jj+3]);
        fprintf(out,"torsion %d %d %d %d  0.2279931 2 413.2828;\n"
                ,ires[jj],ires[jj+1],ires[jj+2],ires[jj+3]);
        fprintf(out,"torsion %d %d %d %d  0.4148926 3 346.9368;\n"
                ,ires[jj],ires[jj+1],ires[jj+2],ires[jj+3]);
    }
    /*
    	for( jj=0; jj< iatom; jj++)
    	{
    	if( x[jj] != 0.)
    	if( y[jj] != 0.)
    	if( z[jj] != 0.)
    	fprintf(out,"tether 100 %d %f %f %f;\n", ires[jj],
    		x[jj],y[jj],z[jj]);
    	}
    */
    exit(0);
}
