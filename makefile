all: monte

CC = gcc
#Options for development
CFLAGS= -g -Wall  
#Options for development
#CFLAGS= -O3 -Wall

monte: main.o conformation.o overlap.o random.o energy.o sample.o
	$(CC) -o monte main.o conformation.o overlap.o random.o energy.o sample.o mesch/meschach.a -lm
main.o: main.c conformation.h overlap.h random.h energy.h global.h sample.h
	$(CC) $(CFLAGS) -c main.c 
conformation.o: conformation.c conformation.h global.h random.h
	$(CC) $(CFLAGS) -c conformation.c 
overlap.o: overlap.c overlap.h global.h
	$(CC) $(CFLAGS) -c overlap.c
random.o: random.c random.h
	$(CC) $(CFLAGS) -c random.c
energy.o: energy.c energy.h global.h conformation.h
	$(CC) $(CFLAGS) -c energy.c
sample.o: sample.c sample.h global.h
	$(CC) $(CFLAGS) -c sample.c	

clean:
	-rm *.o
