/*conformation.h*/
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>

#define pi 3.14159265358979
#define kB 1.3806505E-23 //Boltzmann constant
#define H 100	//Hooken spring constant.
#define Q 100 //Angel constant.

#define delx 0.5	//max size of MC move
#define dely 0.5
#define delz 0.5
#define delalpha 1
#define delbeta 1

//define struct bead,x,y,z is cord,rad is radii.
typedef struct  {
	double x,y,z,r;
}bead ;

typedef struct  {
	/*neighborhood beads index.*/
	int start;
	int end;
	double len;
} spring;

typedef struct  {
	/*neighborhood springs index.*/
	int start;
	int end;
	double angle;
} ang;

typedef struct  {
	int nbd;
	int nsp;
	int nang;
	bead *beads;
	spring *springs;
	ang *angs;
} confor;

extern confor *confor_get(int nbd,int nsp,int nang);
extern void confor_free(confor *p);
extern void confor_copy(confor *from,confor *to);
#endif
