clean:
	\rm test_out.log XX5_multi_test.pdb
test_out.log:
	../www_iMolecule.py -isdf XX5.sdf -multi 5 -omol2 XX5_multi_test.mol2 -logFile test_out.log -wrkPath . -v
	../www_iMolecule.py -isdf XX5.sdf -multi 5 -opdb XX5_multi_test.pdb -logFile test_out.log -wrkPath . -v
tests: test_out.log	

