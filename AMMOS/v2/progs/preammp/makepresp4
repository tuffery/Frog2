

CC = cc
CFLAGS =  -O0 -g
LFLAGS = -lm

presp4: presp4.o dobonds.o dohybridsp3.o dotorsion.o doangle.o uffsupport.o
	$(CC) $(CFLAGS) presp4.o doangle.o dobonds.o dohybridsp3.o dotorsion.o uffsupport.o $(LFLAGS) -o presp4

doangle.o : doangle.c  uff.h 
	$(CC) $(CFLAGS) -c doangle.c
dobonds.o : dobonds.c uff.h
	$(CC) $(CFLAGS) -c dobonds.c
dohybridsp3.o : dohybridsp3.c uff.h
	$(CC) $(CFLAGS) -c dohybridsp3.c
dotorsion.o : dotorsion.c uff.h
	$(CC) $(CFLAGS) -c dotorsion.c
presp4.o : presp4.c uff.h
	$(CC) $(CFLAGS) -c presp4.c 
uffsupport.o : uffsupport.c uff.h
	$(CC) $(CFLAGS) -c uffsupport.c 
