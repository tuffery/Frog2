
CC = cc
GC = gcc
CFLAGS =   -O2 -c -g  
CLFLAGS =  -O2 
LFLAGS = -lm -o

all : mol2_to_templ_sp4.o mol2_to_templ_sp5.o renumb_ligand.o active_case2.o active_case3.o active_case4.o active_case5.o save_energy.o ammp_to_mol2.o ammp_to_pdb_case1.o ammp_to_pdb_case2.o ammp_to_pdb_case3.o ammp_to_pdb_case4.o
all : mol2_to_templ_sp4 mol2_to_templ_sp5 renumb_ligand active_case2 active_case3 active_case4 active_case5 save_energy ammp_to_mol2 ammp_to_pdb_case1 ammp_to_pdb_case2 ammp_to_pdb_case3 ammp_to_pdb_case4 

active_case2.o : active_case2.c 
	$(CC) $(CFLAGS) active_case2.c 

active_case2 : active_case2.o
	$(GC) $(CLFLAGS) active_case2.o $(LFLAGS) active_case2

active_case3.o : active_case3.c 
	$(CC) $(CFLAGS) active_case3.c 
	
active_case3 : active_case3.o 
	$(GC) $(CLFLAGS) active_case3.o $(LFLAGS) active_case3

active_case4.o : active_case4.c 
	$(CC) $(CFLAGS) active_case4.c 

active_case4 : active_case4.o 
	$(GC) $(CLFLAGS) active_case4.o $(LFLAGS) active_case4
	
active_case5.o : active_case5.c 
	$(CC) $(CFLAGS) active_case5.c 
	
active_case5 : active_case5.o
	$(GC) $(CLFLAGS) active_case5.o $(LFLAGS) active_case5
	
mol2_to_templ_sp4.o : mol2_to_templ_sp4.c 
	$(CC) $(CFLAGS) mol2_to_templ_sp4.c 
	
mol2_to_templ_sp4 : mol2_to_templ_sp4.o 
	$(GC) $(CLFLAGS) mol2_to_templ_sp4.o $(LFLAGS) mol2_to_templ_sp4

mol2_to_templ_sp5.o : mol2_to_templ_sp5.c 
	$(CC) $(CFLAGS) mol2_to_templ_sp5.c 

mol2_to_templ_sp5 : mol2_to_templ_sp5.o 
	$(GC) $(CLFLAGS) mol2_to_templ_sp5.o $(LFLAGS) mol2_to_templ_sp5
	
renumb_ligand.o : renumb_ligand.c 
	$(CC) $(CFLAGS) renumb_ligand.c 

renumb_ligand : renumb_ligand.o
	$(GC) $(CLFLAGS) renumb_ligand.o $(LFLAGS) renumb_ligand
	
save_energy.o : save_energy.c 
	$(CC) $(CFLAGS) save_energy.c 

save_energy : save_energy.o
	$(GC) $(CLFLAGS) save_energy.o $(LFLAGS) save_energy
	
ammp_to_mol2.o : ammp_to_mol2.c 
	$(CC) $(CFLAGS) ammp_to_mol2.c 
	
ammp_to_mol2 : ammp_to_mol2.o
	$(GC) $(CLFLAGS) ammp_to_mol2.o $(LFLAGS) ammp_to_mol2
	
ammp_to_pdb_case1.o : ammp_to_pdb_case1.c 
	$(CC) $(CFLAGS) ammp_to_pdb_case1.c 

ammp_to_pdb_case1 : ammp_to_pdb_case1.o
	$(GC) $(CLFLAGS) ammp_to_pdb_case1.o $(LFLAGS) ammp_to_pdb_case1

ammp_to_pdb_case2.o : ammp_to_pdb_case2.c 
	$(CC) $(CFLAGS) ammp_to_pdb_case2.c 

ammp_to_pdb_case2 : ammp_to_pdb_case2.o
	$(GC) $(CLFLAGS) ammp_to_pdb_case2.o $(LFLAGS) ammp_to_pdb_case2

ammp_to_pdb_case3.o : ammp_to_pdb_case3.c 
	$(CC) $(CFLAGS) ammp_to_pdb_case3.c 

ammp_to_pdb_case3 : ammp_to_pdb_case3.o
	$(GC) $(CLFLAGS) ammp_to_pdb_case3.o $(LFLAGS) ammp_to_pdb_case3

ammp_to_pdb_case4.o : ammp_to_pdb_case4.c 
	$(CC) $(CFLAGS) ammp_to_pdb_case4.c 

ammp_to_pdb_case4 : ammp_to_pdb_case4.o
	$(GC) $(CLFLAGS) ammp_to_pdb_case4.o $(LFLAGS) ammp_to_pdb_case4
