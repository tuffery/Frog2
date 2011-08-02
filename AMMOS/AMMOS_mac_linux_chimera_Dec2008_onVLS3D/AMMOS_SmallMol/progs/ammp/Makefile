STATIC= -static
CC = gcc
CFLAGS =   -O2 -c -g -DSGI 
CLFLAGS =  -O2 -DSGI
LFLAGS = $(STATIC) -lm
#LFLAGS = -lm



ammp : ammp.o eval.o analyze.o anonbon.o  damp.o box.o cube.o  clone.o central.o n3.o ghost.o noel.o orthorb.o orbit.o dscf.o av5.o dill.o abc.o animate.o atoms.o bfgs.o bonds.o angle.o doubletime.o dipole.o grasp.o gtor.o gene.o dgeom.o gauss.o bell.o gdock.o abuild.o field.o gsdg.o hard.o hybrid.o image.o unonbon.o rectmm.o vnonbon.o debye.o screen.o tgroup.o tset.o tailor.o ttarget.o torsion.o periodic.o optimist.o variable.o kdock.o kohonen.o math.o mom.o monitor.o morse.o normal.o random.o react.o  richard.o swarm.o scf.o shadow.o simplex.o trunc.o trust.o  table.o tether.o restrain.o
	$(CC) $(CLFLAGS) -o ammp ammp.o eval.o analyze.o anonbon.o  damp.o box.o  cube.o clone.o central.o n3.o ghost.o noel.o orthorb.o orbit.o dscf.o av5.o dill.o abc.o animate.o atoms.o bfgs.o bonds.o angle.o doubletime.o dipole.o grasp.o gtor.o gene.o dgeom.o gauss.o bell.o gdock.o abuild.o field.o gsdg.o hard.o hybrid.o  image.o unonbon.o rectmm.o vnonbon.o debye.o screen.o tgroup.o tset.o tailor.o ttarget.o torsion.o periodic.o optimist.o variable.o kdock.o kohonen.o math.o mom.o monitor.o morse.o normal.o random.o react.o  richard.o swarm.o scf.o shadow.o simplex.o trunc.o trust.o table.o tether.o restrain.o $(LFLAGS) 

ammp.o : ammp.c ammp.h
	$(CC) $(CFLAGS)  ammp.c
eval.o : eval.c ammp.h orbit.h
	$(CC) $(CFLAGS)  eval.c
orthorb.o : orthorb.c ammp.h orbit.h
	$(CC)  $(CFLAGS) orthorb.c
orbit.o : orbit.c ammp.h orbit.h
	$(CC)  $(CFLAGS)  orbit.c
dscf.o : dscf.c ammp.h orbit.h
	$(CC)  $(CFLAGS)  dscf.c
animate.o : animate.c ammp.h
	$(CC) $(CFLAGS)  animate.c
analyze.o : analyze.c ammp.h
	$(CC) $(CFLAGS)  analyze.c
anonbon.o : anonbon.c ammp.h
	$(CC) $(CFLAGS)  anonbon.c
dill.o : dill.c ammp.h
	$(CC) $(CFLAGS)  dill.c
abc.o : abc.c ammp.h
	$(CC) $(CFLAGS)  abc.c
av5.o : av5.c ammp.h
	$(CC) $(CFLAGS)  av5.c
box.o : box.c ammp.h
	$(CC) $(CFLAGS)  box.c
damp.o : damp.c ammp.h
	$(CC) $(CFLAGS)  damp.c
cube.o : cube.c ammp.h
	$(CC) $(CFLAGS)  cube.c
clone.o : clone.c ammp.h
	$(CC) $(CFLAGS)  clone.c
central.o : central.c ammp.h
	$(CC) $(CFLAGS)  central.c
n3.o : n3.c ammp.h
	$(CC) $(CFLAGS)  n3.c
ghost.o : ghost.c ammp.h
	$(CC) $(CFLAGS)  ghost.c
noel.o : noel.c ammp.h
	$(CC) $(CFLAGS)  noel.c
atoms.o : atoms.c ammp.h
	$(CC) $(CFLAGS)  atoms.c
bfgs.o : bfgs.c ammp.h
	$(CC) $(CFLAGS)  bfgs.c
bonds.o : bonds.c ammp.h
	$(CC) $(CFLAGS)  bonds.c
angle.o : angle.c ammp.h
	$(CC) $(CFLAGS)  angle.c
image.o : image.c ammp.h
	$(CC) $(CFLAGS)  image.c
unonbon.o : unonbon.c ammp.h
	$(CC) $(CFLAGS)  unonbon.c
rectmm.o : rectmm.c ammp.h
	$(CC) $(CFLAGS)  rectmm.c
vnonbon.o : vnonbon.c ammp.h
	$(CC) $(CFLAGS)  vnonbon.c
debye.o : debye.c ammp.h
	$(CC) $(CFLAGS)  debye.c
screen.o : screen.c ammp.h
	$(CC) $(CFLAGS)  screen.c
tgroup.o : tgroup.c ammp.h
	$(CC) $(CFLAGS)  tgroup.c
tset.o : tset.c ammp.h
	$(CC) $(CFLAGS)  tset.c
tailor.o : tailor.c ammp.h
	$(CC) $(CFLAGS)  tailor.c
ttarget.o : ttarget.c ammp.h
	$(CC) $(CFLAGS)  ttarget.c
torsion.o : torsion.c ammp.h
	$(CC) $(CFLAGS)  torsion.c
doubletime.o : doubletime.c ammp.h
	$(CC) $(CFLAGS)  doubletime.c
dipole.o : dipole.c ammp.h
	$(CC) $(CFLAGS)  dipole.c
grasp.o : grasp.c ammp.h
	$(CC) $(CFLAGS)  grasp.c
gtor.o : gtor.c ammp.h
	$(CC) $(CFLAGS)  gtor.c
gene.o : gene.c ammp.h
	$(CC) $(CFLAGS)  gene.c
dgeom.o : dgeom.c ammp.h
	$(CC) $(CFLAGS)  dgeom.c
gauss.o : gauss.c ammp.h
	$(CC) $(CFLAGS)  gauss.c
bell.o : bell.c ammp.h
	$(CC) $(CFLAGS)  bell.c
gdock.o : gdock.c ammp.h
	$(CC) $(CFLAGS)  gdock.c
abuild.o : abuild.c ammp.h
	$(CC) $(CFLAGS)  abuild.c
field.o : field.c ammp.h
	$(CC) $(CFLAGS)  field.c
gsdg.o : gsdg.c ammp.h
	$(CC) $(CFLAGS)  gsdg.c
hard.o : hard.c ammp.h
	$(CC) $(CFLAGS)  hard.c
hybrid.o : hybrid.c ammp.h
	$(CC) $(CFLAGS)  hybrid.c
optimist.o : optimist.c ammp.h
	$(CC) $(CFLAGS)  optimist.c
variable.o : variable.c ammp.h
	$(CC) $(CFLAGS)  variable.c
kdock.o : kdock.c ammp.h
	$(CC) $(CFLAGS)  kdock.c
kohonen.o : kohonen.c ammp.h
	$(CC) $(CFLAGS)  kohonen.c
math.o : math.c ammp.h
	$(CC) $(CFLAGS)  math.c
mom.o : mom.c ammp.h
	$(CC) $(CFLAGS)  mom.c
monitor.o : monitor.c ammp.h
	$(CC) $(CFLAGS)  monitor.c
morse.o : morse.c ammp.h
	$(CC) $(CFLAGS)  morse.c
normal.o : normal.c ammp.h
	$(CC) $(CFLAGS)  normal.c
random.o : random.c ammp.h
	$(CC) $(CFLAGS)  random.c
react.o : react.c ammp.h
	$(CC) $(CFLAGS)  react.c
richard.o : richard.c ammp.h
	$(CC) $(CFLAGS) richard.c
swarm.o : swarm.c ammp.h 
	$(CC) $(CFLAGS)  swarm.c
scf.o : scf.c ammp.h orbit.h
	$(CC) $(CFLAGS)  scf.c
shadow.o : shadow.c ammp.h
	$(CC) $(CFLAGS)  shadow.c
simplex.o : simplex.c ammp.h
	$(CC) $(CFLAGS)  simplex.c
trunc.o : trunc.c ammp.h
	$(CC) $(CFLAGS)  trunc.c
trust.o : trust.c ammp.h
	$(CC) $(CFLAGS)  trust.c
table.o : table.c ammp.h
	$(CC) $(CFLAGS)  table.c
tether.o : tether.c ammp.h
	$(CC) $(CFLAGS)  tether.c
restrain.o : restrain.c ammp.h
	$(CC) $(CFLAGS)  restrain.c
periodic.o : periodic.c ammp.h
	$(CC) $(CFLAGS)  periodic.c
