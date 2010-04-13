/*conformation.h*/
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <math.h>
#include <strings.h>

#define pi 3.14159265358979
#define kB 1.3806505E-23 //Boltzmann constant
#define H 100	//Hooken spring constant.
#define Q 100 //Angel constant.
#define nmax 200 //max number of beads a conformation could have.
#define delx 0.5	//max size of MC move
#define dely 0.5
#define delz 0.5

//define struct bead,x,y,z is cord,rad is radii.
struct bead {
	double x,y,z,r;
} ;
struct spring {
	double x0,x1;
	double y0,y1;
	double z0,z1;
	double len;
};

struct confor {
	struct bead beads[nmax];
	struct spring springs[nmax];
} ;

#endif
