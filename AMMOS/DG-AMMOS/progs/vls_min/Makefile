STATIC= -static
CC = cc
GC = gcc
CFLAGS =   -O2 -c -g  
CLFLAGS =  -O2 
LFLAGS = $(STATIC) -lm -o

all : mol2_to_templ_sp4.o ammp_to_mol2.o save_energy.o 
all : mol2_to_templ_sp4 ammp_to_mol2 save_energy  


mol2_to_templ_sp4.o : mol2_to_templ_sp4.c 
	$(CC) $(CFLAGS) mol2_to_templ_sp4.c 

mol2_to_templ_sp4 : mol2_to_templ_sp4.o 
	$(GC) $(CLFLAGS) mol2_to_templ_sp4.o $(LFLAGS) mol2_to_templ_sp4

ammp_to_mol2.o : ammp_to_mol2.c 
	$(CC) $(CFLAGS) ammp_to_mol2.c 

ammp_to_mol2 : ammp_to_mol2.o
	$(GC) $(CLFLAGS) ammp_to_mol2.o $(LFLAGS) ammp_to_mol2

save_energy.o : save_energy.c 
	$(CC) $(CFLAGS) save_energy.c 

save_energy : save_energy.o
	$(GC) $(CLFLAGS) save_energy.o $(LFLAGS) save_energy
