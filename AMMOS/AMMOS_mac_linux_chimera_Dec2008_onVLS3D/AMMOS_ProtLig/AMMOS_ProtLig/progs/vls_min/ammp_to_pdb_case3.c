///////////////////////////////////////////////////
///                                             ///
///  Read data from output file for protein from ammp   ///
///          written by T. Pencheva, M.Miteva   ///
///                 2007-2008                   ///
///                                             ///
///////////////////////////////////////////////////


// ---------------------------------------------------------------------------------------------
// Read data from output file for protein from ammp and replace coordinates in the initial *.pdb


#include <stdio.h>		
#include <string.h>	 
#include <stdlib.h>	

#define max_numb_atoms 120000
#define max_name_length 10

// Function upper converts characters in string a to upper case
void upper (char *str, int len)
{
	int i;
	for (i = 0; i < len; i++) str [i] = toupper (str[i]); 
}

main()
{

int	a = 0, i = 0, j = 0, 
	k = 0, l = 0, p = 0, 			// Counters  
	ammp_numb_at [max_numb_atoms],		// Current number of atom in *.ammp
	ammp_numb_curr,				// Current number of atom in *.ammp
	ammp_res_numb [max_numb_atoms],		// Current residue in *.ammp
	pdb_numb_at,				// Current number of atom in *.pdb
	pdb_numb_at_output,
	pdb_numb_res [max_numb_atoms],		// Number of residue in *.pdb
	pdb_numb_res_curr,			// Number of residue in *.pdb
	line_length = 1000; 			// Length of the line  

float	pdb_coordx, pdb_coordy, pdb_coordz,	// Atoms coordinates in the initial *.pdb
	ammp_currx, ammp_curry, ammp_currz,	// Atoms coordinates read from *.ammp  
	ammp_coordx [max_numb_atoms], 
	ammp_coordy [max_numb_atoms], 
	ammp_coordz [max_numb_atoms],		
	pdb_n1_curr, pdb_n2_curr,		// Used when read from *.pdb
	pdb_n1 [max_numb_atoms],
	pdb_n2 [max_numb_atoms],
	ammp_n1, ammp_n2, ammp_n3, ammp_n4;	// Used when read from *.ammp


char	*ammp_name [max_numb_atoms],		// Name of the atom in *.ammp
	ammp_name_curr [max_name_length],	
	*ammp_res [max_numb_atoms],		// Name of the residue in *.ammp
	*pdb_name [max_numb_atoms],	
	*pdb_name_comp [max_numb_atoms],	// Name of the atom in *.pdb
	pdb_name_curr [max_name_length],	
	*pdb_res [max_numb_atoms],		// Name of the residue in *.pdb
	pdb_res_curr [max_name_length],		
	*pdb_chain [max_numb_atoms],		// Name of the chain in *.pdb
	pdb_chain_curr [max_name_length],	
	pdb_atom_type_curr [max_name_length],	// Atom type in *.pdb
	*pdb_atom_type [max_numb_atoms],		// Atom type in *.pdb
	output_word [max_name_length], 		// For "ATOM" when read from *.ammp
	pdb_word [max_name_length], 		// For "ATOM" when read from *.pdb
	output_buffer [line_length],		// Buffer for *.ammp file
	pdb_buffer [line_length];		// Buffer for *.pdb file


// Initialization of pointers for word atom and atoms names
for (i = 0; i < max_numb_atoms; i++) 
{
	ammp_name [i] = (char *) malloc (max_name_length);
	ammp_res [i] = (char *) malloc (max_name_length);
	pdb_res [i] = (char *) malloc (max_name_length);
	pdb_name [i] = (char *) malloc (max_name_length);
	pdb_chain [i] = (char *) malloc (max_name_length);
	pdb_name_comp [i] = (char *) malloc (max_name_length);
	pdb_atom_type [i] = (char *) malloc (max_name_length);
}


FILE 	*ammp, *pdb_in, *pdb_out, *fopen(); 

// Open working files for reading and writting  
ammp = fopen ("output.ammp", "r");		// File contains optimized structure of protein after running AMMP
pdb_in = fopen ("input_protein.pdb", "r"); 	// Reading necessary information from *.pdb
pdb_out = fopen ("output_protein.pdb", "w"); 	// Writting new coordinates from *.ammp into new *.pdb file 

// Reading data from input_protein.pdb before starting ATOM information
fgets (pdb_buffer, line_length, pdb_in);			
sscanf (pdb_buffer, "%s", pdb_word);
while (strncmp (pdb_word, "ATOM", 4) != 0 ) 
{	
	fgets (pdb_buffer, line_length, pdb_in);
	sscanf (pdb_buffer, "%s", pdb_word);
}

a = 0;
for (j = 0; j < max_numb_atoms; j++) 
{
	sscanf (pdb_buffer, "%s %i %s %s %s %i %f %f %f %f %f %s", pdb_word, &pdb_numb_at, pdb_name_curr, pdb_res_curr, 	
		pdb_chain_curr, &pdb_numb_res_curr, &pdb_coordx, &pdb_coordy, &pdb_coordz, &pdb_n1_curr, &pdb_n2_curr, pdb_atom_type_curr);


// printf ("%s %i %s %s %s %i %f %f %f %f %f %s\n", pdb_word, pdb_numb_at, pdb_name_curr, pdb_res_curr, 	
//	pdb_chain_curr, pdb_numb_res_curr, pdb_coordx, pdb_coordy, pdb_coordz, pdb_n1_curr, pdb_n2_curr, pdb_atom_type_curr);

	if (!strncmp (pdb_word, "END", 3)) break;

	else 
	{
		strcpy (pdb_name [a], pdb_name_curr);
		strcpy (pdb_res [a], pdb_res_curr);
		strcpy (pdb_chain [a], pdb_chain_curr);
		strcpy (pdb_atom_type [a], pdb_atom_type_curr);
		strcat (pdb_name_comp [a], pdb_res_curr);
		strcat (pdb_name_comp [a], ".");
		strcat (pdb_name_comp [a], pdb_name_curr);
		pdb_numb_res [a] = pdb_numb_res_curr;
		pdb_n1 [a] = pdb_n1_curr;
		pdb_n2 [a] = pdb_n2_curr;
		a++;
	}

fgets (pdb_buffer, line_length, pdb_in);
}


// Reading data for atoms from output.ammp - only new coordinates of active atoms of the protein are needed
for (i = 0; i < max_numb_atoms; i++)
{
	if (!strncmp (output_word, "atom", 4)) 
	{
		sscanf (output_buffer, "%s %f %f %f %i %s %f %f %f %f", output_word, 
			&ammp_currx, &ammp_curry, &ammp_currz, &ammp_numb_curr,
			ammp_name_curr, &ammp_n1, &ammp_n2, &ammp_n3, &ammp_n4);

		if (ammp_numb_curr >= 1000000) 
		{
			l = 1;
			fgets (output_buffer, line_length, ammp); 
			sscanf (output_buffer, "%s", output_word); 
			continue;
		}	
		else l = 0;
	}		

	fgets (output_buffer, line_length, ammp); 
	sscanf (output_buffer, "%s", output_word); 

	if (!strncmp (output_word, "mompar", 6)) 
	{
		fgets (output_buffer, line_length, ammp); 
		sscanf (output_buffer, "%s", output_word); 
	}

	if (!strncmp (output_word, "inactive", 8)) 
	{
		fgets (output_buffer, line_length, ammp); 
		sscanf (output_buffer, "%s", output_word); 
		if (!strncmp (output_word, "bond", 4)) break;
		continue;
	}
		
	else
	{
		if (l == 1) continue;
		else
		{
			ammp_coordx [p] = ammp_currx;
			ammp_coordy [p] = ammp_curry;
			ammp_coordz [p] = ammp_currz;
			ammp_numb_at [p] = ammp_numb_curr;
			div_t numb_res;
			numb_res = div (ammp_numb_curr, 100);
			ammp_res_numb [p] = numb_res.quot;
			strcpy (ammp_name [p], ammp_name_curr);
			upper (ammp_name [p], 10); 
			strncpy (ammp_res [p], ammp_name [p], 3);
			p++;
			l = 0;
		}	
	}

	if (!strncmp (output_word, "atom", 4)) continue;

	if (!strncmp (output_word, "bond", 4)) 
	{
			ammp_coordx [p] = ammp_currx;
			ammp_coordy [p] = ammp_curry;
			ammp_coordz [p] = ammp_currz;
			ammp_numb_at [p] = ammp_numb_curr;
			strcpy (ammp_name [p], ammp_name_curr);
			upper (ammp_name [p], 10); 
			strncpy (ammp_res [p], ammp_name [p], 3);
			break;
	}

	fgets (output_buffer, line_length, ammp);
	sscanf (output_buffer, "%s", output_word); 
}

pdb_numb_at_output = 1;
l = 0;
for (i = 0; i < a; i++)
{
	if ( l == p - 1) break;

	if (pdb_numb_res [i] == ammp_res_numb [l+1] && strcmp (pdb_name_comp [i], ammp_name [l+1]) == 0)
	{
		fprintf (pdb_out, "%4s%7i%5s%4s%2s%5i%11.3f%8.3f%8.3f%6.2f%6.2f%12s\n", 
			 "ATOM", pdb_numb_at_output, pdb_name [i], pdb_res [i], pdb_chain [i], pdb_numb_res [i], 	
			 ammp_coordx [l+1], ammp_coordy [l+1], ammp_coordz [l+1], pdb_n1 [i], pdb_n2 [i], pdb_atom_type [i]);
		pdb_numb_at_output++;
		l++;
		continue;
	}

	if (pdb_numb_res [i] != ammp_res_numb [l+1] && ammp_res_numb [l+1] == ammp_res_numb [l]) 
	{
		fprintf (pdb_out, "%4s%7i%5s%4s%2s%5i%11.3f%8.3f%8.3f%6s%6s\n", 
			 "ATOM", pdb_numb_at_output, &ammp_name [l+1][4], ammp_res [l+1], pdb_chain [i-1], ammp_res_numb [l+1], 
			 ammp_coordx [l+1], ammp_coordy [l+1], ammp_coordz [l+1], "1.00", "99.99");
		pdb_numb_at_output++;
		l++;
		i--;
		continue;
	}

	else continue;

}
fprintf (pdb_out, "%s", "END");


// Free located memory for atom_name
for (i = 0; i < max_numb_atoms; ++i) 
{
	free (ammp_name [i]);
	free (ammp_res [i]);
	free (pdb_res [i]);
	free (pdb_name [i]);
	free (pdb_chain [i]);
	free (pdb_name_comp [i]);
	free (pdb_atom_type [i]);
}

// Close working files
fclose (pdb_in);		 
fclose (pdb_out);		
fclose (ammp);		  

}

