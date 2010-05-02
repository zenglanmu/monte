#include "global.h"
#include "energy.h"
#include "conformation.h"

extern int ntotal;
extern int nspring,nang;
extern double temp;
extern double *Edist,*Eang;	//equilibrium spring lens and equilibrium angel.
extern confor *InitialConf;

double EBond(const confor *p)
//Bond potential.
{
	int i;
	double E=0;
	
	for(i=0;i<nspring;i++)
		E += kB*temp*H*pow((p->springs[i].len-Edist[i]),2);
		
	return E;
}

double EAng(const confor *p)
//ang potential
{
	int i;
	double E=0;
	
	for(i=0;i<nang;i++)
		E += kB*temp*Q*pow((p->angs[i].angle - Eang[i]),2);
		
	return E;
}
		
double EVpair(const confor *p)
//Lennard-Jones potential
{
	int i,j;
	double r,E;
	double epsilon=0.3,sigma=1.8,rmax=7.2;
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
	double r,E,nu,D,rD,kappa,l;
	double x,y,z;
	
	rD=3.07e-9;
	kappa=1/rD;
	nu=0.243;
	D=80;
	E=0;
	
	for(i=0;i<ntotal-1;i++){
		for(j=i+1;j<ntotal;j++){
			x = p->beads[i].x - p->beads[j].x;
			y = p->beads[i].y - p->beads[j].y;
			z = p->beads[i].z - p->beads[j].z;
			r = sqrt(x*x+y*y+z*z);
					
			x = InitialConf->beads[i].x - InitialConf->beads[j].x;
			y = InitialConf->beads[i].y - InitialConf->beads[j].y;
			z = InitialConf->beads[i].z - InitialConf->beads[j].z;	
			l = sqrt(x*x+y*y+z*z);
			
			E +=  nu*nu*l*l/D*exp(-kappa*r)/r;
		}	
	}
	return E;
}

double Energy(const confor *p)
{
	double Etotal;
	//Etotal = EBond(p)+EAng(p)+EVpair(p)+CHpair(p);
	Etotal = EVpair(p)+CHpair(p);
	return Etotal;
}
