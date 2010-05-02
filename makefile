all: monte

CC = gcc
#Options for development
CFLAGS= -gdwarf-2 -g3 -Wall  
#Options for development
#CFLAGS= -O3 -Wall

monte: main.o conformation.o overlap.o random.o energy.o sample.o io.o
	$(CC) -o monte main.o conformation.o overlap.o random.o energy.o sample.o io.o  mesch/meschach.a -lm
main.o: main.c conformation.h overlap.h random.h energy.h global.h sample.h io.h
	$(CC) $(CFLAGS) -c main.c 
conformation.o: conformation.c conformation.h global.h random.h
	$(CC) $(CFLAGS) -c conformation.c 
overlap.o: overlap.c overlap.h global.h
	$(CC) $(CFLAGS) -c overlap.c
random.o: random.c random.h
	$(CC) $(CFLAGS) -c random.c
energy.o: energy.c energy.h global.h conformation.h
	$(CC) $(CFLAGS) -c energy.c
sample.o: sample.c sample.h global.h mesch/matrix.h
	$(CC) $(CFLAGS) -c sample.c
io.o: io.c io.h global.h
	$(CC) $(CFLAGS) -c io.c

clean:
	-rm *.o
