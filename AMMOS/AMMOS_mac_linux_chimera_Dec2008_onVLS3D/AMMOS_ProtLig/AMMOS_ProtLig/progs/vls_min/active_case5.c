///////////////////////////////////////////////////
///                                             ///
///         Definition of active atoms          ///
///         written by T. Pencheva, M.Miteva    ///
///                 2007-2008                   ///
///                                             ///
///////////////////////////////////////////////////


// ---------------------------------------------------------------------
// Definition of active atoms in the sphere around the ligand
// The whole protein is inactive, only the ligand is active


#include <stdio.h>		
#include <stdlib.h>
#include <math.h>	

#define max_numb_atoms 250000
#define max_name_length 15

main()
{
int	i = 0, j = 0, k = 0, n = 0, 
	p = 0, l = 0, flag = 0,		// Counters 
	end_atom = 0, 			// For number of atoms in the protein
	line_length = 1000, 		// Length of the line  
	res_atom_number = 0,		// Counter for current atom in the residue
	lig_atom_number = 0,		// Counter for current atom in the ligand
	all_active [max_numb_atoms],	// Array for active atoms in the sphere around the ligand
	active [max_numb_atoms],	// Array for active atoms after remove of repeating
	lig_atoms [max_numb_atoms],	// Array for all atoms of the ligand
	res_atoms [max_numb_atoms];	// Array for all atoms of the protein

float	radius = 0,					// Radius of the sphere around the ligand
	dist = 0,						
	res_coordx, res_coordy, res_coordz, 		  
	lig_coordx, lig_coordy, lig_coordz, 		  
	res_numb1, res_numb2, res_numb3, res_numb4,	
	lig_numb1, lig_numb2, lig_numb3, lig_numb4;	 

char	res_buf [line_length],		// Buffer with line_length number of chars 
	lig_buf [line_length],		// Define string sufficiently large to store a line of input  
	res_word [max_name_length],	 
	lig_word [max_name_length],	   
	*res_name [max_numb_atoms],	  
	*lig_name [max_numb_atoms];	   	

// Initialization of pointers for word atom and atoms names
	for (j = 0; j < max_numb_atoms; j++) 
	{
	    res_name [j] = (char *) malloc (max_name_length);
	    lig_name [j] = (char *) malloc (max_name_length);
	}

FILE 	*protein, *ligand, *numb_act, *fopen(); 

// Open working files for reading and writting  
	ligand = fopen ("input_ligand.ammp", "r"); 	
	protein = fopen ("protein.ammp", "r"); 	
	numb_act = fopen ("active.ammp", "w");		


// Reading from *.ammp file for ligand

for (i = 0; i < max_numb_atoms; i++)
{
fgets (lig_buf, line_length, ligand);	
sscanf (lig_buf, "%s", lig_word);	

	if (!strncmp (lig_word, "bond", 4)) break; 
	if (!strncmp (lig_word, "mompar", 6)) 
	{
		fgets (lig_buf, line_length, ligand); 
		sscanf (lig_buf, "%s", lig_word); 
	}
	if (!strncmp (lig_word, "inactive", 8)) 
	{
		fgets (lig_buf, line_length, ligand); 
		sscanf (lig_buf, "%s", lig_word); 
	}
	if (!strncmp (lig_word, "atom", 4)) 
	{
		sscanf (lig_buf, "%s %f %f %f %i %s %f %f %f %f", lig_word, &lig_coordx, &lig_coordy, &lig_coordz, 
			&lig_atom_number, lig_name [k], &lig_numb1, &lig_numb2, &lig_numb3, &lig_numb4);
	
		lig_atoms [l] = lig_atom_number;
		l++;
	}
}

// Reading from *.ammp file for protein

for (j = 0; j < max_numb_atoms; j++)
{
fgets (res_buf, line_length, protein);	
sscanf (res_buf, "%s", res_word);	

	if (!strncmp (res_word, "bond", 4)) break;

	if (!strncmp (res_word, "mompar", 6)) 
	{	
		fgets (res_buf, line_length, protein); 
		sscanf (res_buf, "%s", res_word); 
	}
	if (!strncmp (res_word, "atom", 4)) 
	{
		sscanf (res_buf, "%s %f %f %f %i %s %f %f %f %f", res_word, &res_coordx, &res_coordy, &res_coordz, 
			&res_atom_number, res_name [i], &res_numb1, &res_numb2, &res_numb3, &res_numb4);

			res_atoms [p] = res_atom_number;
			p ++;
	}
}


// Print of active atoms
for (j = 0; j < l; j ++) fprintf (numb_act, "%s %i %i%s\n", "active", lig_atoms [j], lig_atoms [j], ";");


// Free located memory for atom_name
for (i = 0; i < max_numb_atoms; ++i) 
{
	free (res_name [i]);
	free (lig_name [i]);
}


// Close working files
fclose (protein);	
fclose (ligand);
fclose (numb_act);	
}

