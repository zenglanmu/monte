/*conformation.h*/
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <math.h>
#include <strings.h>

#define pi 3.14159265358979
#define kB 1.3806505E-23 //Boltzmann constant
#define nmax 2000 //max number of beads a conformation could have.
#define delx 1	//max size of MC move
#define dely 1
#define delz 1

//define struct bead,x,y,z is cord,rad is radii.
struct bead {
	double x,y,z,r;
} ;
struct spring {
	double len;
};

struct confor {
	struct bead beads[nmax];
	//struct spring springs[nmax-1];
} ;

#endif
