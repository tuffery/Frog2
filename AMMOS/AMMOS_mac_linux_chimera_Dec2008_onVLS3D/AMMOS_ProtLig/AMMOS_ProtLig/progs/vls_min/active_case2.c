///////////////////////////////////////////////////
///                                             ///
///        Definition of inactive atoms         ///
///          written by T. Pencheva, M.Miteva   ///
///                 2007-2008                   ///
///                                             ///
///////////////////////////////////////////////////


// ---------------------------------------------------------------------
// Definition of inactive atoms in the proteins - the atoms from backbone

#include <stdio.h>		
#include <stdlib.h>

#define max_numb_atoms 400000
#define max_name_length 15

main()
{
int	i = 0,	j = 0,			// Counters 
	line_length = 1000, 		// Length of the line  
	res_atom_number = 0;		// Counter for current atom in the residue
	
float	coordx, coordy, coordz, 	// Atoms coordinates  
	numb1, numb2, numb3, numb4; 	// Four fields after atoms (residue) name
	
char	buffer [line_length],		  
					  
	word [max_name_length],		  
	*name [max_numb_atoms];		  
 	
// Initialization of pointers for word atom and atoms names
	for (j = 0; j < max_numb_atoms; j++)
	    name [j] = (char *) malloc (max_name_length);


FILE 	*protein, *numb_inact, *fopen(); 

// Open working files for reading and writting  
	protein = fopen ("protein.ammp", "r"); 	// Open working *.ammp file for reading
	numb_inact = fopen ("inactive.ammp", "w");	// Create .txt file for writting

	

// Reading from *.ammp file for protein

while (strncmp (word, "bond", 4))
{
	fgets (buffer, line_length, protein);	
	sscanf (buffer, "%s", word);	

	if (!strncmp (word, "mompar", 6)) 
	{	
		fgets (buffer, line_length, protein); 
		sscanf (buffer, "%s", word); 
	}
	if (!strncmp (word, "atom", 4)) 
	{
		sscanf (buffer, "%s %f %f %f %i %s %f %f %f %f", word, &coordx, &coordy, &coordz, 
			&res_atom_number, name [i], &numb1, &numb2, &numb3, &numb4);

// Check for the backbone atoms - N, C alfa, C, O and H at N
 		if (!strcmp (&name [i][4], "n") && !strcmp (&name [i][5], "")) 
		{	
			fprintf (numb_inact, "%s %i %i%s\n", "inactive", res_atom_number, res_atom_number, ";");
			continue;
		}

		if (!strcmp (&name [i][4], "ca") && !strcmp (&name [i][6], "")) 
		{
			fprintf (numb_inact, "%s %i %i%s\n", "inactive", res_atom_number, res_atom_number, ";");
			continue;
		}

		if (!strcmp (&name [i][4], "c") && !strcmp (&name [i][5], "")) 
		{
			fprintf (numb_inact, "%s %i %i%s\n", "inactive", res_atom_number, res_atom_number, ";");
			continue;
		}

		if (!strcmp (&name [i][4], "o") && !strcmp (&name [i][5], "")) 
		{
			fprintf (numb_inact, "%s %i %i%s\n", "inactive", res_atom_number, res_atom_number, ";");
			continue;
		}

		if (!strcmp (&name [i][4], "h") && !strcmp (&name [i][5], "")) 
		{	
			fprintf (numb_inact, "%s %i %i%s\n", "inactive", res_atom_number, res_atom_number, ";");
			continue;
		}
		if (!strcmp (&name [i][4], "ha") && !strcmp (&name [i][6], "")) 
		{	
			fprintf (numb_inact, "%s %i %i%s\n", "inactive", res_atom_number, res_atom_number, ";");
			continue;
		}
		if (!strcmp (&name [i][4], "1ha") && !strcmp (&name [i][7], "")) 
		{	
			fprintf (numb_inact, "%s %i %i%s\n", "inactive", res_atom_number, res_atom_number, ";");
			continue;
		}
		if (!strcmp (&name [i][4], "2ha") && !strcmp (&name [i][7], "")) 
		{	
			fprintf (numb_inact, "%s %i %i%s\n", "inactive", res_atom_number, res_atom_number, ";");
			continue;
		}
	}
}


// Free located memory for atom_name
for (i = 0; i < max_numb_atoms; ++i) free (name [i]);


// Close working files
fclose (protein);		 
fclose (numb_inact);		   

}

