///////////////////////////////////////////////////
///                                             ///
///                 VLS_AMMP                    ///
///          written by T. Pencheva             ///
///                 2006-2007                   ///
///                                             ///
///////////////////////////////////////////////////



// -------------------------------------------------------------
// Read data from output file from ammp and replace coordinates in the initial *.mol2


#include <stdio.h>	// Allows the program to interact with the screen, keyboard and filesystem of the computer	
#include <string.h>	// For concatenation of two strings 
#include <stdlib.h>

#define max_numb_atoms 1000
#define max_name_len 10


main()
{

int	i = 0, j = 0, k = 0, 
	l = 0, p = 0,		  	// Counters  
	atoms_number, 			// Number of atoms
	bonds_number, 			// Number of bonds
	number,				// Current number of atom 
	atom1, atom2,			// Used when read bonds
	line_length = 1000; 		// Length of the line  

float	old_coordx, old_coordy, old_coordz,	// Atoms coordinates in the initial *.mol2
	new_coordx [max_numb_atoms], 	// Atoms coordinates read from ammp  
	new_coordy [max_numb_atoms], 	// Atoms coordinates read from ammp  
	new_coordz [max_numb_atoms], 	// Atoms coordinates read from ammp 
	numb1, numb2, numb3, numb4,	// Used when read from *.ammp
 	charge;				// Atoms charge


char	type1 [max_name_len],		// Type of the atom - column 1  
	type2 [max_name_len],		// Type of the atom - column 2  
	type3 [max_name_len],		// Type of the atom - column 3
	word [max_name_len], 		// For "ATOM" when read from *.ammp
	name_mol2 [max_name_len],	// Name of the atom in *.mol2 file
	name_ammp [max_name_len],	// Name of the atom in *.ammp file  
	bond [max_name_len],		// For bond
	mol2_buffer [line_length],	// Buffer for *.mol2 fiel with line_length number of chars - max width of one row  
					// Define string sufficiently large to store a line of input  
	ammp_buffer [line_length];	// Buffer for *.ammp file



FILE 	*ammp, *mol2, *mol2_out, *fopen(); 

// Open working files for reading and writting  
	ammp = fopen ("output.ammp", "r");
	mol2 = fopen ("input.mol2", "r"); 	// Reading necessary information from *.mol2
	mol2_out = fopen ("output.mol2", "w"); 	// Writting new coordinates from *.ammp into new *.mol2 file 

// Reading from *.mol2 file and writing first rows (before "@<TRIPOS>ATOM") in new *.mol2 file
	fgets (mol2_buffer, line_length, mol2);			
	fprintf (mol2_out, "%s", mol2_buffer);
	while (strncmp (mol2_buffer, "@<TRIPOS>MOLECULE", 17) != 0 ) 
	{
		fgets (mol2_buffer, line_length, mol2);
		fprintf (mol2_out, "%s", mol2_buffer);
	}
	fgets (mol2_buffer, line_length, mol2);	
	fprintf (mol2_out, "%s", mol2_buffer); 
	fgets (mol2_buffer, line_length, mol2);	
	fprintf (mol2_out, "%s", mol2_buffer);
	sscanf (mol2_buffer, "%i %i", &atoms_number, &bonds_number);	// Reading atoms and bonds numbers 
	fgets (mol2_buffer, line_length, mol2);
	fprintf (mol2_out, "%s", mol2_buffer);
	while (strncmp (mol2_buffer, "@<TRIPOS>ATOM", 13) != 0 ) 	
	{
		fgets (mol2_buffer, line_length, mol2);
    		fprintf (mol2_out, "%s", mol2_buffer);
	}

// Reading data for atoms from output.ammp - only new coordinates of the ligand are needed
for (i = 0; i < max_numb_atoms; i++)
{
	if (number == atoms_number - 1) break;

	if (strncmp (word, "atom", 4) || strncmp (word, "mompar", 6)) 
	{
		fgets (ammp_buffer, line_length, ammp);
		sscanf (ammp_buffer, "%s", word); 
	}

	if (!strncmp (word, "mompar", 6)) 
	{
		fgets (ammp_buffer, line_length, ammp); 
		sscanf (ammp_buffer, "%s", word); 
	}
 
	if (!strncmp (word, "atom", 4)) 
	{
		sscanf (ammp_buffer, "%s %f %f %f %i %s %f %f %f %f", word, &new_coordx [i], &new_coordy [i], &new_coordz [i], 
			&number, name_ammp , &numb1, &numb2, &numb3, &numb4);
// printf("%f %f %f\n", new_coordx [i], new_coordy [i], new_coordz [i]);
	}
}

// Writing data for atoms in new *.mol2
	for (k = 0; k < atoms_number; k++)
	{
		fgets (mol2_buffer, line_length, mol2);	
		sscanf (mol2_buffer, "%i %s %f %f %f %s %s %s %f", &l, name_mol2, 
			&old_coordx, &old_coordy, &old_coordz, type1, type2, type3, &charge);
		
		fprintf (mol2_out, "%7i %-5s%13.4f%10.4f%10.4f %-5s%6s%4s%15.4f\n",  
			 l, name_mol2, new_coordx [k], new_coordy [k], new_coordz [k], type1, type2, type3, charge);
	}

// Reading data for bonds from old *.mol2 and write it in new *.mol2
	fgets (mol2_buffer, line_length, mol2);
 	fprintf (mol2_out, "%s", mol2_buffer); 	// For "@<TRIPOS>BOND"
	for (j = 0; j < bonds_number; j++)
	{
		fgets (mol2_buffer, line_length, mol2);
		sscanf (mol2_buffer, "%i %i %i %s", &l, &atom1, &atom2, bond);
		fprintf (mol2_out, "%6i%5i%5i %s\n", l, atom1, atom2, bond);  
	}

// Close working files
fclose (mol2);		// Close *.mol2 file 
fclose (mol2_out);	// Close new *.mol2 file, with results from ammp
fclose (ammp);		// Close output.ammp file  

}

