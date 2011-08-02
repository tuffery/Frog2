all: ammos frogc
mac: ammosmac frogcmac
ammos:
	cd AMMOS && make -e all
frogc:
	cd iMolecule/code_c && make -e all
clean:
	cd iMolecule/code_c && make clean
