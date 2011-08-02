///////////////////////////////////////////////////
///                                             ///
///  Definition of active atoms in the sphere   ///
///          written by T. Pencheva, M.Miteva   ///
///                 2007-2008                   ///
///                                             ///
///////////////////////////////////////////////////



// ---------------------------------------------------------------------
// Definition of active atoms in the sphere around the ligand
// All atoms in the sphere will be active, the whole rest protein will be inactive 
// Active will be not only the atoms, failed in the sphere, but also all atoms from the touched residues, 
// excepts the backbone of the protein in the sphere

#include <stdio.h>		
#include <stdlib.h>
#include <math.h>	

#define max_numb_atoms 160000
#define max_name_length 15
#define plus 1000000

// Perfoming of the bubble algorithm to sort an integer array
void bubble_sort (int *arr, int n)
{
	int     i, j;
  	int     temp;

  	/* Let i run over all entries (except the last). */
  	for (i = 0; i < n - 1; i++)
  	{
    		/* Let j run over all entries from (i+1) until the last. */
    		for (j = i + 1; j < n; j++)
    		{
      			/* If arr[i] > arr[j], swap the two entries. */
      			if (arr[i] > arr[j])
      			{
        			temp = arr[i];
        			arr[i] = arr[j];
        			arr[j] = temp;
      			}
    		}
  	}

return;
}

main()
{
int	i = 0, j = 0, k = 0, l = 0, m = 0, a = 0, b = 0, s = 0,
	n = 0, q = 0, p = 0, flag = 0,	// Counters 
	end_atom = 0, 			// For number of atoms in the protein
	line_length = 1000, 		// Length of the line  
	res_atom_number = 0,		// Counter for current atom in the residue
	lig_atom_number = 0,		// Counter for current atom in the ligand
	all_active [max_numb_atoms],	// Array for active atoms in the sphere around the ligand
	active [max_numb_atoms],	// Array for active atoms after remove of repeating
	inact_backbone [max_numb_atoms],// Array for inactive atoms in the backbone	
	backbone [max_numb_atoms],	// Array for inactive atoms in the backbone
	inactive [max_numb_atoms],	// Array for inactive atoms	
	res_atoms [max_numb_atoms],	// Array for all atoms of the protein
	res_pointers [max_numb_atoms],	
	pointers [max_numb_atoms],	
	all_sphere [max_numb_atoms];	
	

float	radius = 0,					// Radius of the sphere around the ligand
	dist = 0,					// Need for determination of the sphere around the ligand	
	res_coordx, res_coordy, res_coordz, 		// Residue coordinates  
	lig_coordx, lig_coordy, lig_coordz, 		// Residue coordinates  
	res_numb1, res_numb2, res_numb3, res_numb4,	// Four fields after atoms name in residue
	lig_numb1, lig_numb2, lig_numb3, lig_numb4;	// Four fields after atoms name in ligand 

char	res_buf [line_length],		// Buffer with line_length number of chars - max width of one row
	lig_buf [line_length],		// Define string sufficiently large to store a line of input  
	param_buf [line_length],	// Define string sufficiently large to store a line of input  
	res_word [max_name_length],	 
	lig_word [max_name_length],	   
	param_word [max_name_length],	   
	*res_name [max_numb_atoms],	 
	*lig_name [max_numb_atoms];	// Name of the atom   	

// Initialization of pointers for word atom and atoms names
for (j = 0; j < max_numb_atoms; j++) 
{
    res_name [j] = (char *) malloc (max_name_length);
    lig_name [j] = (char *) malloc (max_name_length);
}

FILE 	*protein, *ligand, *numb_act, *param, *fopen(); 

// Open working files for reading and writting  
// protein.ammp will be open later in order to be open each time for each ligand atom
ligand = fopen ("input_ligand.ammp", "r"); 	// Open working *.ammp file for ligand for reading
numb_act = fopen ("active.ammp", "w");		// Create .txt file for writting of active atoms
param = fopen ("param.temp", "r"); 		// Open ammp.param file, created by user for reading of desired value of the radius
	
// Input from the keyboard of the desired value of sphere radius 
// printf ("Input a desired value of the sphere radius: ");
// scanf ( "%f", &radius);

// Input the desired value of sphere radius from the parameter file created by user
for (i = 0; i < 5; i++) fgets (param_buf, line_length, param);
sscanf (param_buf, "%s %f", param_word, &radius);

// Reading from *.ammp file for ligand
for (i = 0; i < 200; i++)
{
	flag ++;
	protein = fopen ("protein.ammp", "r"); 		// Open working *.ammp file for protein for reading
	fgets (lig_buf, line_length, ligand);	
	sscanf (lig_buf, "%s", lig_word);	

	if (!strncmp (lig_word, "bond", 4)) break; 
	if (!strncmp (lig_word, "mompar", 6)) 
	{
		fgets (lig_buf, line_length, ligand); 
		sscanf (lig_buf, "%s", lig_word); 
	}
	if (!strncmp (lig_word, "atom", 4)) 
	{
		sscanf (lig_buf, "%s %f %f %f %i %s %f %f %f %f", lig_word, &lig_coordx, &lig_coordy, &lig_coordz, 
			&lig_atom_number, lig_name [k], &lig_numb1, &lig_numb2, &lig_numb3, &lig_numb4);
		q++;
	}


// Reading from *.ammp file for protein
	for (j = 0; j < max_numb_atoms; j++)
	{
		fgets (res_buf, line_length, protein);	
		sscanf (res_buf, "%s", res_word);	

		if (!strncmp (res_word, "bond", 4)) 
		{
// Determination of the number of atoms in the protein (end_atoms)
			if (flag == 1) 
			{
				end_atom = p;
				break; 
			}
			else break; 
		}
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

// Check for the backbone atoms - N, C alfa, C, O and H at N			
			if (!strcmp (&res_name [i][4], "n") && !strcmp (&res_name [i][5], ""))
			{	
				inact_backbone [b] =  res_atom_number;
				b++;
			}
			if (!strcmp (&res_name [i][4], "ca") && !strcmp (&res_name [i][6], "")) 
			{	
				inact_backbone [b] =  res_atom_number;
				b++;
			}
			if (!strcmp (&res_name [i][4], "c") && !strcmp (&res_name [i][5], "")) 
			{	
				inact_backbone [b] =  res_atom_number;
				b++;
			}
			if (!strcmp (&res_name [i][4], "o") && !strcmp (&res_name [i][5], "")) 
			{	
				inact_backbone [b] =  res_atom_number;
				b++;
			}
			if (!strcmp (&res_name [i][4], "h") && !strcmp (&res_name [i][5], "")) 
			{	
				inact_backbone [b] =  res_atom_number;
				b++;
			}
			if (!strcmp (&res_name [i][4], "ha") && !strcmp (&res_name [i][6], "")) 
			{	
				inact_backbone [b] =  res_atom_number;
				b++;
			}
			if (!strcmp (&res_name [i][4], "1ha") && !strcmp (&res_name [i][7], "")) 
			{	
				inact_backbone [b] =  res_atom_number;
				b++;
			}
			if (!strcmp (&res_name [i][4], "2ha") && !strcmp (&res_name [i][7], "")) 
			{	
				inact_backbone [b] =  res_atom_number;
				b++;
			}

// Check if the atom from the resudie "fails" into the sphere with determined radius
			dist = sqrt ((res_coordx - lig_coordx)*(res_coordx - lig_coordx) + 
			       (res_coordy - lig_coordy)*(res_coordy - lig_coordy) +	
		  	       (res_coordz - lig_coordz)*(res_coordz - lig_coordz));

			if (dist < radius) 
			{	   
// Print of all active atoms, with repeating ones
				all_active [n] = res_atom_number;
				n ++;
			}
			else continue;
		}
		continue;
	}
}


// Remove of repeating numbers from the all_active array
k = 0;	// Initialization of the internal counter

for (i = 0; i < n; )
{
	j = 0; 
	while (all_active [i] != active [j] && j <= k) j ++;
	if (j < k) 
	{	
		i ++;
		continue;
	}	
	active [k] = all_active [i];
	k ++;
}

// Determination of the residues, which atoms fails in the sphere
for (i = 0; i < k; i ++) 
{
// 	itoa (active [i], res_pointers [i], 10);
	div_t point;
	point = div (active [i], 100);
	res_pointers [i] =  point.quot;
}

// Remove of repeating numbers from the res_pointers array and sort them
k = 0;	// Initialization of the internal counter

for (i = 0; i < n; )
{
	j = 0; 
	while (res_pointers [i] != pointers [j] && j <= k) j ++;
	if (j < k) 
	{	
		i ++;
		continue;
	}
	pointers [k] = res_pointers [i];
	k ++;
}

l = k - 1;
bubble_sort (pointers, l);

// Remove of repeating numbers from the inact_backbone array
k = 0;	// Initialization of the internal counter

for (i = 0; i < n; )
{
	j = 0; 
	while (inact_backbone [i] != backbone [j] && j <= k) j ++;
	if (j < k) 
	{	
		i ++;
		continue;
	}	
	backbone [k] = inact_backbone [i];
	k ++;
}

p = k - 1;
bubble_sort (backbone, p);


// Taking those atoms from backbone, failed particularly in the sphere
s = 0;
for (i = 0; i < l; i ++)  
{
	for (j = 0; j < p ; j ++)  
	{
		div_t inact;
		inact = div (backbone [j], 100);
		if (pointers [i] == inact.quot)
		{
			inactive [s] = backbone [j];
			s++;
		}
		else continue;
	}
}

// Including of all atoms from the residues, failed particularly in the sphere, except the atoms from backbone
m = 0;
b = 0;
for (i = 0; i < l; i ++)  
{
	for (j = 0; j < end_atom ; j ++)  
	{
		div_t numb;
		numb = div (res_atoms [j], 100);
		if ( pointers [i] == numb.quot)
		{
			if (inactive [b] == res_atoms [j])
			{
				b++;
				continue;
			}
			all_sphere [m]	= res_atoms [j];
			m++;
		}
		else continue;
	}
}


// Print of active atoms for the ligand - add a value of "plus" to their number
for (j = 0; j < q; j ++) fprintf (numb_act, "%s %i %i%s\n", "active", plus + j, plus + j, ";");
// Print of active atoms of protein, without repeating ones
for (i = 0; i < m; i ++) fprintf (numb_act, "%s %i %i%s\n", "active", all_sphere [i], all_sphere [i], ";");


// Free located memory for atom_name
for (i = 0; i < max_numb_atoms; ++i) 
{
	free (res_name [i]);
	free (lig_name [i]);
}


// Close working files
fclose (protein);	// Close *.ammp file for protein
fclose (ligand);	// Close *.ammp file for ligand
fclose (numb_act);	// Close *.txt file for active atoms
fclose (param);		// Close ammp.param file
}

