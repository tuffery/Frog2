STATIC=-static
OPT=-O2
all: TM QBestFit AxeRot monoconfV3 mol2clusterize
monoconfV3: monoconfV3.c argstr.c Mol2.c TxtFile.c Zuker.c Types.h Zuker.c RMSd.c 
	gcc $(STATIC) $(OPT) -o monoconfV3 monoconfV3.c argstr.c TxtFile.c Mol2.c Zuker.c RMSd.c -lm
mol2clusterize: mol2clusterize.c argstr.c Mol2.c TxtFile.c Zuker.c Types.h RMSd.c
	gcc $(STATIC) $(OPT) -o mol2clusterize mol2clusterize.c argstr.c Mol2.c TxtFile.c Zuker.c RMSd.c -lm
TM: TM.c
	gcc $(STATIC) $(OPT) -o TM TM.c -lm
QBestFit: QBestFit.c
	gcc $(STATIC) $(OPT) -o QBestFit QBestFit.c -lm
AxeRot: AxeRot.c
	gcc $(STATIC) $(OPT) -o AxeRot AxeRot.c -lm
clean:
	\rm TM QBestFit AxeRot mol2clusterize *~ monoconfV3
