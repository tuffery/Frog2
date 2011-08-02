///////////////////////////////////////////////////
///                                             ///
/// Read data from output file for protein from ammp   ///
///          written by T. Pencheva, M.Miteva          ///
///                 2007-2008                   ///
///                                             ///
///////////////////////////////////////////////////


// ---------------------------------------------------------------------------------------------
// Read data from output file for protein from ammp and replace coordinates in the initial *.pdb


#include <stdio.h>		
#include <string.h>	 
#include <stdlib.h> 

#define max_numb_atoms 200000
#define max_name_length 10

// Function lower converts characters in string a to lower case
void lower (char *str, int len)
{
	int i;
	for (i = 0; i < len; i++) str [i] = tolower (str[i]); 
}

void upper (char *str, int len)
{
	int i;
	for (i = 0; i < len; i++) str [i] = toupper (str[i]); 
}

main()
{

int	i = 0, l = 0, p = 0,  	// Counters  
	ammp_numb_at [max_numb_atoms],		// Current number of atom in *.ammp
	pdb_numb_at,				// Current number of atom in *.pdb
	pdb_numb_at_output,
	pdb_numb_res [max_numb_atoms],		// Number of residue in *.pdb
	line_length = 1000; 			// Length of the line  

float	pdb_coordx, pdb_coordy, pdb_coordz,	// Atoms coordinates in the initial *.pdb
	ammp_coordx [max_numb_atoms], 
	ammp_coordy [max_numb_atoms], 
	ammp_coordz [max_numb_atoms],		// Atoms coordinates read from *.ammp  
	pdb_n1, pdb_n2,				// Used when read from *.pdb
	ammp_n1, ammp_n2, ammp_n3, ammp_n4;	// Used when read from *.ammp


char	*ammp_name [max_numb_atoms],		// Name of the atom in *.ammp
	*pdb_name [max_numb_atoms],		// Name of the atom in *.pdb
	*pdb_name_comp [max_numb_atoms],	// Name of the residue in *.pdb
	*pdb_res [max_numb_atoms],		// Name of the residue in *.pdb
	pdb_chain  [max_name_length],		// Name of the chain in *.pdb
	pdb_atom_type [max_name_length],	// Atom type in *.pdb
	output_word [max_name_length], 		// For "ATOM" when read from *.ammp
	pdb_word [max_name_length], 		// For "ATOM" when read from *.pdb
	output_buffer [line_length],		// Buffer for *.ammp file
	pdb_buffer [line_length];		// Buffer for *.pdb file


// Initialization of pointers for word atom and atoms names
for (i = 0; i < max_numb_atoms; i++) 
{
	ammp_name [i] = (char *) malloc (max_name_length);
	pdb_res [i] = (char *) malloc (max_name_length);
	pdb_name [i] = (char *) malloc (max_name_length);
	pdb_name_comp [i] = (char *) malloc (max_name_length);
}


FILE 	*ammp, *pdb_in, *pdb_out, *fopen(); 

// Open working files for reading and writting  
ammp = fopen ("output.ammp", "r");		// File contains optimized structure of protein after running AMMP
pdb_in = fopen ("input_protein.pdb", "r"); 	// Reading necessary information from *.pdb
pdb_out = fopen ("output_protein.pdb", "w"); 	// Writting new coordinates from *.ammp into new *.pdb file 


// Reading data for atoms from output.ammp - only new coordinates of active atoms of the protein are needed
for (i = 0; i < max_numb_atoms; i++)
{
	fgets (output_buffer, line_length, ammp);
	sscanf (output_buffer, "%s", output_word); 

	if (!strncmp (output_word, "bond", 4)) break;

	if (!strncmp (output_word, "mompar", 6)) 
	{
		fgets (output_buffer, line_length, ammp); 
		sscanf (output_buffer, "%s", output_word); 
	}

	if (!strncmp (output_word, "atom", 4)) 
	{
		sscanf (output_buffer, "%s %f %f %f %i %s %f %f %f %f", output_word, 
			&ammp_coordx [i], &ammp_coordy [i], &ammp_coordz [i], &ammp_numb_at [i],
			ammp_name [i], &ammp_n1, &ammp_n2, &ammp_n3, &ammp_n4);
		p++;
	}
}


// Reading data from input_protein.pdb before starting ATOM information
fgets (pdb_buffer, line_length, pdb_in);			
sscanf (pdb_buffer, "%s", pdb_word);
while (strncmp (pdb_word, "ATOM", 4) != 0 ) 
{	
	fgets (pdb_buffer, line_length, pdb_in);
	sscanf (pdb_buffer, "%s", pdb_word);
}

l = 0;
pdb_numb_at_output = 1;
for (i = 0; i < p; i++) 
{
	sscanf (pdb_buffer, "%s %i %s %s %s %i %f %f %f %f %f %s", pdb_word, &pdb_numb_at, pdb_name [i], pdb_res [i], 	
		pdb_chain, &pdb_numb_res [i], &pdb_coordx, &pdb_coordy, &pdb_coordz, &pdb_n1, &pdb_n2, pdb_atom_type);

	if (strncmp (pdb_word, "END", 3) == 0) 
	{
		if (l == p) break;
		else
		{
			while (l!= p)
			{
				upper (ammp_name [l], 10); 
				fprintf (pdb_out, "%4s%7i%5s%4s%2s%5i%11.3f%8.3f%8.3f%6s%6s\n", "ATOM", pdb_numb_at_output, 
					 &ammp_name [l][4], pdb_res [i - 1], pdb_chain,pdb_numb_res [i - 1], 
					 ammp_coordx [l], ammp_coordy [l], ammp_coordz [l], "1.00", "99.99");
				pdb_numb_at_output++;
				l++;
			}
		continue;
		}
	}

	else
	{
		strcat (pdb_name_comp [i], pdb_res [i]);
		strcat (pdb_name_comp [i], ".");
		strcat (pdb_name_comp [i], pdb_name [i]);
		lower (pdb_name_comp [i], 10); 

		if (ammp_numb_at [l] >= 1000000) 
		{
			l++;
			continue;
		}

		if (strcmp (pdb_name_comp [i], ammp_name [l]) == 0)
		{
			if (strncmp (pdb_word, "HETATM", 6) == 0) 
			{	
				if (strncmp (pdb_chain, "", 0) == 0) strcpy (pdb_chain, "A");
								
				fprintf (pdb_out, "%6s%5i%5s%4s%2s%5i%11.3f%8.3f%8.3f%6.2f%6.2f%12s\n", 
					 pdb_word, pdb_numb_at_output, pdb_name [i], pdb_res [i], pdb_chain, ammp_numb_at [l]/100, 	
					 ammp_coordx [l], ammp_coordy [l], ammp_coordz [l], pdb_n1, pdb_n2, pdb_atom_type);
				pdb_numb_at_output++;
				l++;
			}
			
			else
			{
				fprintf (pdb_out, "%4s%7i%5s%4s%2s%5i%11.3f%8.3f%8.3f%6.2f%6.2f%12s\n", 
					 pdb_word, pdb_numb_at_output, pdb_name [i], pdb_res [i], pdb_chain, pdb_numb_res [i], 	
					 ammp_coordx [l], ammp_coordy [l], ammp_coordz [l], pdb_n1, pdb_n2, pdb_atom_type);
				pdb_numb_at_output++;
				l++;
			}
		}
	
		else 
		{
			while ((strcmp (pdb_name_comp [i], ammp_name [l]) != 0))
			{
				upper (ammp_name [l], 10); 
				fprintf (pdb_out, "%4s%7i%5s%4s%2s%5i%11.3f%8.3f%8.3f%6s%6s\n", pdb_word, pdb_numb_at_output, 
					 &ammp_name [l][4], pdb_res [i - 1], pdb_chain, pdb_numb_res [i - 1], 
					 ammp_coordx [l], ammp_coordy [l], ammp_coordz [l], "1.00", "99.99");
				pdb_numb_at_output++;
				l++;
			}
			continue;
		}
	}
	fgets (pdb_buffer, line_length, pdb_in);	
}

fprintf (pdb_out, "%s", "END"); 


// Free located memory for atom_name
for (i = 0; i < max_numb_atoms; ++i) 
{
	free (ammp_name [i]);
	free (pdb_res [i]);
	free (pdb_name [i]);
	free (pdb_name_comp [i]);
}

// Close working files
fclose (pdb_in);		// Close *.pdb file 
fclose (pdb_out);		// Close new *.pdb file, with results from ammp
fclose (ammp);			// Close output.ammp file  

}

