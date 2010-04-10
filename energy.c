#include "global.h"
#include "energy.h"

extern int ntotal;

double LJ(struct confor p)
//Lennard-Jones potential
{
	int i,j;
	double r,E,epsilon,sigma,rmax;
	epsilon=0.3;sigma=1.8;rmax=5.2;
	E=0;
	for(i=0;i<ntotal;i++){
		for(j=i;j<ntotal-1;j++){
			r=sqrt(pow(((p.beads[j].x)-(p.beads[j+1].x)),2)+pow(((p.beads[j].y)-(p.beads[j+1].y)),2)+pow(((p.beads[j].z)-(p.beads[j+1].z)),2));
			if(r>rmax)
				E+=0;
			else
				E+=4*epsilon*(pow((sigma/r),12)-pow((sigma/r),6));
		}	
	}
	return E;
}
double Energy(struct confor p)
{
	double Etotal;
	Etotal=LJ(p);
	return Etotal;
}
