all: monte

CC = gcc
#Options for development
CFLAGS= -g -Wall 
#Options for development
#CFLAGS= -O -Wall 
monte: main.o conformation.o overlap.o random.o
	$(CC) -o monte main.o conformation.o overlap.o -lm
main.o: main.c conformation.h global.h
	$(CC) $(CFLAGS) -c main.c 
conformation.o: conformation.c conformation.h global.h random.h
	$(CC) $(CFLAGS) -c conformation.c 
overlap.o: overlap.c overlap.h global.h
	$(CC) $(CFLAGS) -c overlap.c
random.o: random.c random.h
	$(CC) -c random.c

clean:
	-rm *.o
