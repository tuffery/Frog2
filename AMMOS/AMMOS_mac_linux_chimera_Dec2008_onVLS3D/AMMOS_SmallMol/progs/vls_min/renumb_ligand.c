///////////////////////////////////////////////////
///                                             ///
///                 AMMOS_SmallMol              ///
///          written by T. Pencheva, M.Miteva   ///
///                 2007-2008                   ///
///                                             ///
///////////////////////////////////////////////////


// ------------------------------------------------------------------------------------------------------------
// Re-numbering of molecule atoms after preammp 
// Add value of "plus" to the atom numbers


#include <stdio.h>		
#include <stdlib.h>	

#define max_numb_atoms 10000
#define max_name_len 10
#define plus 0 

main()
{

int	i = 0, j = 0, k = 0, 			// Counters  
	atom_number,				// Used when read for "atom" from *.ammp
	mompar_number,				// Used when read for "mompar" from *.ammp
	bond_atom1, bond_atom2,			// Used when read for "bond" from *.ammp
	angle_atom1, angle_atom2, angle_atom3,	// Used when read for "angle" from *.ammp
	hybrid_atom1, hybrid_atom2, 
	hybrid_atom3, hybrid_atom4,		// Used when read for "hybrid" from *.ammp
	torsion_atom1 [max_numb_atoms],   
	torsion_atom2 [max_numb_atoms],
	torsion_atom3 [max_numb_atoms], 
	torsion_atom4 [max_numb_atoms], 
	torsion_n2,				// Used when read for "torsion" from *.ammp
	line_length = 1000; 			// Length of the line  

float	atom_coordx, atom_coordy, atom_coordz,	// Used when read for "atom" from *.ammp 
	atom_n1, atom_n2, atom_n3, atom_n4,	// Used when read for "atom" from *.ammp
	mompar_n1, mompar_n2,			// Used when read for "mompar" from *.ammp
	bond_n1, bond_n2,			// Used when read for "bond" from *.ammp
	angle_n1, angle_n2,			// Used when read for "angle" from *.ammp
	hybrid_n1, hybrid_n2,			// Used when read for "hybrid" from *.ammp
	torsion_n1, torsion_n3;			// Used when read for "torsion" from *.ammp


char	word [max_name_len], 		// For the first word when read from *.ammp
	atom_name [max_name_len],	// Name of the atom when read for "atom" from *.ammp 
	buffer [line_length];		// Buffer for *.ammp file


FILE 	*ammp, *ammp_curr, *fopen(); 

// Open working files for reading and writting  
	ammp = fopen ("input_ligand.ammp", "r");	// File for ligand after running preammp
	ammp_curr = fopen ("ammp_curr.ammp", "w");	// Temporary file

// Reading data from *.ammp
for (i = 0; i < max_numb_atoms; i++)
{
	fgets (buffer, line_length, ammp); 
	sscanf (buffer, "%s", word); 

	if (!strncmp (word, "atom", 4)) 
	{
		sscanf (buffer, "%s %f %f %f %i %s %f %f %f %f", word, 
			&atom_coordx, &atom_coordy, &atom_coordz, &atom_number,
			atom_name, &atom_n1, &atom_n2, &atom_n3, &atom_n4);
		atom_number = atom_number + plus;
		fprintf (ammp_curr, "%s %f %f %f %i %s %f %f %f %f %s\n", word, 
			atom_coordx, atom_coordy, atom_coordz, atom_number,
			atom_name, atom_n1, atom_n2, atom_n3, atom_n4, ";");
	}

	if (!strncmp (word, "mompar", 6)) 
	{
		sscanf (buffer, "%s %i %f %f", word, &mompar_number, &mompar_n1, &mompar_n2);
		mompar_number = mompar_number + plus;
		fprintf (ammp_curr, "%s %i %f %f %s\n", word, mompar_number, mompar_n1, mompar_n2, ";");
	}

	if (!strncmp (word, "bond", 4)) 
	{
		sscanf (buffer, "%s %i %i %f %f", word, &bond_atom1, &bond_atom2, &bond_n1, &bond_n2);
		bond_atom1 = bond_atom1 + plus;
		bond_atom2 = bond_atom2 + plus;
		fprintf (ammp_curr, "%s %i %i %f %f %s\n", word, bond_atom1, bond_atom2, bond_n1, bond_n2, ";");
	}

	if (!strncmp (word, "angle", 5)) 
	{
		sscanf (buffer, "%s %i %i %i %f %f", word, &angle_atom1, &angle_atom2, &angle_atom3, &angle_n1, &angle_n2);
		angle_atom1 = angle_atom1 + plus;
		angle_atom2 = angle_atom2 + plus;
		angle_atom3 = angle_atom3 + plus;
		fprintf (ammp_curr, "%s %i %i %i %f %f %s\n", word, angle_atom1, angle_atom2, angle_atom3, angle_n1, angle_n2, ";");
	}

	if (!strncmp (word, "hybrid", 6)) 
	{
		sscanf (buffer, "%s %i %i %i %i %f %f", word, &hybrid_atom1, &hybrid_atom2, &hybrid_atom3, &hybrid_atom4, 
			&hybrid_n1, &hybrid_n2);
		hybrid_atom1 = hybrid_atom1 + plus;
		hybrid_atom2 = hybrid_atom2 + plus;
		hybrid_atom3 = hybrid_atom3 + plus;
		hybrid_atom4 = hybrid_atom4 + plus;
		fprintf (ammp_curr, "%s %i %i %i %i %f %f %s\n", word, hybrid_atom1, hybrid_atom2, 
			 hybrid_atom3, hybrid_atom4, hybrid_n1, hybrid_n2, ";");
	}

	if (!strncmp (word, "torsion", 7)) 
	{
		sscanf (buffer, "%s %i %i %i %i %f %i %f", word, &torsion_atom1 [i], &torsion_atom2[i], 
			&torsion_atom3 [i], &torsion_atom4 [i], &torsion_n1, &torsion_n2, &torsion_n3);
		torsion_atom1 [i] = torsion_atom1 [i] + plus;
		torsion_atom2 [i] = torsion_atom2 [i] + plus;
		torsion_atom3 [i] = torsion_atom3 [i] + plus;
		torsion_atom4 [i] = torsion_atom4 [i] + plus;

		if (i > 1 &&  torsion_atom1 [i] == torsion_atom1 [i-1] &&
			      torsion_atom2 [i] == torsion_atom2 [i-1] &&
			      torsion_atom3 [i] == torsion_atom3 [i-1] &&
			      torsion_atom4 [i] == torsion_atom4 [i-1]) break;

		fprintf (ammp_curr, "%s %i %i %i %i %f %i %f %s\n", word, torsion_atom1 [i], torsion_atom2 [i], 
			 torsion_atom3 [i], torsion_atom4 [i], torsion_n1, torsion_n2, torsion_n3, ";");
	}

}


// Close working files
fclose (ammp);			  
fclose (ammp_curr);		// Close temporary *.ammp file  
rename ("ammp_curr.ammp", "input_ligand.ammp"); // Save re-numbering with the same file name
remove ("ammp_curr.ammp");	// Remove temporary file ammp_active

}

