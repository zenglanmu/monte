#include "global.h"
#include "energy.h"
#include "conformation.h"

extern int ntotal;
extern int nspring,nang;
extern double temp;

double EBond(const confor *p)
//Bond potential.
{
	const double H = 10;	//Hooken spring constant.
	const double l0 = 0.72; //equilibrium spring lens
	int i;
	double E=0;
	
	for(i=0;i<nspring;i++)
		E += kB*temp*H*pow((p->springs[i].len-l0),2);
		
	return E;
}

double EAng(const confor *p)
//ang potential
{
	const double Q = 10; //Angel constant.
	const double ang0 = 1; //equilibrium angle
	int i;
	double E=0;
	
	for(i=0;i<nang;i++)
		E += kB*temp*Q*pow((p->angs[i].angle - ang0),2);
		
	return E;
}
		
double EVpair(const confor *p)
//Lennard-Jones potential
{
	int i,j;
	double r,E;
	double epsilon=0.3,sigma=3.2,rmax=7.2;
	double x,y,z;
	
	E = 0;
	
	for(i=0;i<ntotal-1;i++){
		for(j=i+1;j<ntotal;j++){
			x = p->beads[i].x - p->beads[j].x;
			y = p->beads[i].y - p->beads[j].y;
			z = p->beads[i].z - p->beads[j].z;
			r = sqrt(x*x+y*y+z*z);
			
			if(r>rmax)
				E +=  0;
			else
				E +=  4*epsilon*(pow((sigma/r),12)-pow((sigma/r),6));
		}	
	}
	return E;
}

double CHpair(const confor *p)
//Debye Huckel potential
{
	int i,j;

	const double l=0.072; //equilibrium lens.
	const double rD=3.07e-9;
	const double kappa=1/rD;
	const double nu=0.243;
	const double D=800;
	double x,y,z,r,E;
	
	E=0;
	
	for(i=0;i<ntotal-1;i++){
		for(j=i+1;j<ntotal;j++){
			x = p->beads[i].x - p->beads[j].x;
			y = p->beads[i].y - p->beads[j].y;
			z = p->beads[i].z - p->beads[j].z;
			r = sqrt(x*x+y*y+z*z);
			
			E +=  nu*nu*l*l/D*exp(-kappa*r)/r;
		}	
	}
	return E;
}

double Energy(const confor *p)
{
	double Etotal;
	Etotal = EBond(p)+EAng(p)+EVpair(p)+CHpair(p);
	return Etotal;
}
