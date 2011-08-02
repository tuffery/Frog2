# Installation of package AMMOS_ProtLig

EXE = mol2_to_templ_sp4 mol2_to_templ_sp5 renumb_ligand active_case2 active_case3 active_case4 active_case5 save_energy ammp_to_mol2 ammp_to_pdb_case1 ammp_to_pdb_case2 ammp_to_pdb_case3 ammp_to_pdb_case4
PY =  MinProteinLigands_sp4.py MinProteinLigands_sp5.py AMMOS_ProtLig_sp4.py AMMOS_ProtLig_sp5.py BankDivider.py

BIN = ../../bin

all :
	echo "Installation of AMMP program"
	cd ../progs/ammp; rm *.o; make ammp; cp ammp ammp_nongraph; mv ammp_nongraph $(BIN);
	echo ""; echo ""; echo ""

			
	echo "Installation of PREAMMP program"	
	cd ../progs/preammp; rm *.o; make presp4; cp presp4 preammp; mv preammp $(BIN);
	echo ""; echo ""; echo "";

	
	echo "Installation of package AMMOS_ProtLig"
	cd ../progs/vls_min; rm *.o; make all
	cd ../progs/vls_min; mv $(EXE) $(BIN);
	cd ../progs/Python_scripts; cp $(PY) $(BIN);
	
