# Installation of package VLS_AMMP_BuildMolecules.

EXE = mol2_to_templ_sp4 ammp_to_mol2 save_energy 
PY =  dg-ammos.py 
BIN = ../../bin

all :
	echo "Installation of AMMP program"
	cd ../progs/ammp; rm *.o; make -e ammp; cp ammp ammp_nongraph; mv ammp_nongraph $(BIN);
	echo ""; echo ""; echo ""


	echo "Installation of PREAMMP program"	
	cd ../progs/preammp; rm *.o; make -e presp4; cp presp4 preammp; mv preammp $(BIN);
	echo ""; echo ""; echo "";

	echo "Installation of package VLS_AMMP_BuildMolecules"
	cd ../progs/vls_min; rm *.o; make -e all
	cd ../progs/vls_min; cp $(EXE) $(BIN);
	cd ../progs/Python_scripts; cp $(PY) $(BIN);
