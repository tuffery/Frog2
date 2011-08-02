#!/bin/sh

#This script is used for generate 3D coordinates of cycles from all library's given in argument. These library's
#must be in mol2 file format.
#Different variables need to be fixed in creaLibCycles.py.
#Options:
#  -m -M: Create a new library, without watch if a cycle already exists in the current library.




if test -e erreurs.txt; then
	`rm erreurs.txt`
fi

case $1 in
	"-M"|"-m") arg="True"; shift;;
	-) echo option $1 invalide; return;;
	*) arg="False";;
esac

for repertoire in $*; do
	echo "in $repertoire"
	for fichier in `ls $repertoire | grep mol2.gz`; do
		echo "   processing $fichier"
		`cp $repertoire/$fichier ./TMP.mol2.gz`
		echo "     unzipping library"
		`gunzip -q TMP.mol2.gz`
 		echo "         DONE"
		`python creaLibCycles.py $arg TMP.mol2`
		`rm TMP.mol2`
	done
done
