///////////////////////////////////////////////////
///                                             ///
///                AMMOS_SmallMol               ///
///          written by T. Pencheva, M.Miteva   ///
///                 2006-2007                   ///
///                                             ///
///////////////////////////////////////////////////


// ----------------------------------------------------------------
// Save energies from two files - before_opt.txt and after_opt.txt 
// and save only needed information in one new file energy.txt


#include <stdio.h>		
#include <string.h>	 
#include <stdlib.h>

#define max_numb_atoms 1000
#define max_name_len 10


main()
{

int	i = 0, j = 0, k = 0,  		// Counters  
	line_length = 1000, 		// Length of the line  
	atom_numb;			// For the number of the atom

float	numb_1, numb_2, numb_3,		// For the float numbers in the information from analyze
	bef_int, bef_ext, bef_total,	// For the different energies
	aft_int, aft_ext, aft_total,
	delta_energy;

char	buffer [line_length],		// Buffer with line_length number of chars - max width of one row  

	*word_1 [max_numb_atoms],	// For the first word in the line  
	*word_2 [max_numb_atoms],	
	*word_3 [max_numb_atoms],	
	*word_4 [max_numb_atoms],	
	*word_5 [max_numb_atoms],	
	*word_6 [max_numb_atoms];	
		

// Initialization of pointers for atoms names and types and for bonds types
for (i = 0; i < max_numb_atoms; i++) 
{
	word_1 [i] = (char *) malloc (line_length);
	word_2 [i] = (char *) malloc (line_length);
	word_3 [i] = (char *) malloc (line_length);
	word_4 [i] = (char *) malloc (line_length);
	word_5 [i] = (char *) malloc (line_length);
	word_6 [i] = (char *) malloc (line_length);
}

FILE 	*before, *after, *final, *fopen(); 

// Open working files for reading and writting  
	before = fopen ("before_opt.txt", "r"); 
	after = fopen ("after_opt.txt", "r" );
	final = fopen ("energy.txt", "w" );


// Reading from file before_opt.txt 
for (j = 0; j < max_numb_atoms; j++)
{
	fgets (buffer, line_length, before);	
	sscanf (buffer, "%s %s %s %i %s %f %s %f %s %f", word_1 [j], word_2 [j], word_3 [j], 
		&atom_numb, word_4 [j], &numb_1, word_5 [j], &numb_2, word_6 [j], &numb_3);	
	if (strncmp (word_3 [j], "internal", 8) == 0) break;
	else continue;
}	

	sscanf (buffer, "%s %s %s %f", word_1 [j], word_2 [j], word_3 [j], &bef_int);
	fgets (buffer, line_length, before);	
	sscanf (buffer, "%s %s %s %f", word_1 [j], word_2 [j], word_3 [j], &bef_ext);
	fgets (buffer, line_length, before);	
	sscanf (buffer, "%s %s %f", word_1 [j], word_2 [j], &bef_total);

// Reading from file after_opt.txt 
for (j = 0; j < max_numb_atoms; j++)
{
	fgets (buffer, line_length, after);	
	sscanf (buffer, "%s %s %s %i %s %f %s %f %s %f", word_1 [j], word_2 [j], word_3 [j], 
		&atom_numb, word_4 [j], &numb_1, word_5 [j], &numb_2, word_6 [j], &numb_3);	
	if (strncmp (word_3 [j], "internal", 8) == 0) break;
	else continue;
}	
	sscanf (buffer, "%s %s %s %f", word_1 [j], word_2 [j], word_3 [j], &aft_int);
	fgets (buffer, line_length, after);	
	sscanf (buffer, "%s %s %s %f", word_1 [j], word_2 [j], word_3 [j], &aft_ext);
	fgets (buffer, line_length, after);	
	sscanf (buffer, "%s %s %f", word_1 [j], word_2 [j], &aft_total);

	delta_energy = bef_int - aft_int;

// Writting in final file

	fprintf (final, " %f %s %f %s %f", 
		 bef_int, ":", aft_int, ":", delta_energy);


// Before termination, free up space and close files
// Free located for atom_name and type1 memory
for (i = 0; i < max_numb_atoms; ++i)
{
	free (word_1 [i]); 
	free (word_2 [i]);
	free (word_3 [i]);
	free (word_4 [i]); 
	free (word_5 [i]);
	free (word_6 [i]);
}


// Close working files
fclose (before);	// Close file before_opt.txt
fclose (after);		// Close file after_opt.txt   
fclose (final);		// Close file energy.txt  
// remove ("before_opt.txt");	// Remove temporary files
// remove ("after_opt.txt");	// Remove temporary files
}

