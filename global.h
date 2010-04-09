/*conformation.h*/
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <math.h>
#include <strings.h>

#define pi 3.14159265358979
#define nmax 2000	//max number of beads a conformation could have.

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
