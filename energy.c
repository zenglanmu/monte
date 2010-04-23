#include "global.h"
#include "energy.h"
#include "conformation.h"

extern int ntotal;
extern int nspring,nang;
extern double temp;
extern double Edist[nmax],Eang[nmax];	//equilibrium spring lens and equilibrium angel.
extern struct confor InitialConf;

double EBond(struct confor p)
//Bond potential.
{
	int i;
	double E=0;
	for(i=0;i<nspring;i++) 	E+=	kB*temp*H*pow((p.springs[i].len-Edist[i]),2);
	return E;
}
double EAng(struct confor p)
//ang potential
{
	int i;
	double E=0,ang[nmax];
	GetAngel(p,ang);
	for(i=0;i<nang;i++)	E+=	kB*temp*Q*pow((ang[i]-Eang[i]),2);
	return E;
}
		
double EVpair(struct confor p)
//Lennard-Jones potential
{
	int i,j;
	double r,E,epsilon,sigma,rmax;
	epsilon=0.3;sigma=1.8;rmax=7.2;
	E=0;
	for(i=0;i<ntotal-1;i++){
		for(j=i+1;j<ntotal;j++){
			r=sqrt(pow(((p.beads[i].x)-(p.beads[j].x)),2)+pow(((p.beads[i].y)-(p.beads[j].y)),2)+pow(((p.beads[i].z)-(p.beads[j].z)),2));
			if(r>rmax)
				E+=0;
			else
				E+=4*epsilon*(pow((sigma/r),12)-pow((sigma/r),6));
		}	
	}
	return E;
}

double CHpair(struct confor p)
//Debye Huckel potential
{
	int i,j;
	double r,E,nu,D,rD,kappa,l;
	rD=3.07e-9;
	kappa=1/rD;
	nu=0.243;
	D=80;
	E=0;
	for(i=0;i<ntotal-1;i++){
		for(j=i+1;j<ntotal;j++){
			r=sqrt(pow(((p.beads[i].x)-(p.beads[j].x)),2)+pow(((p.beads[i].y)-(p.beads[j].y)),2)+pow(((p.beads[i].z)-(p.beads[j].z)),2));
			l=sqrt(pow(((InitialConf.beads[i].x)-(InitialConf.beads[j].x)),2)+pow(((InitialConf.beads[i].y)-(InitialConf.beads[j].y)),2)+pow(((InitialConf.beads[i].z)-(InitialConf.beads[j].z)),2));
	
			E+=nu*nu*l*l/D*exp(-kappa*r)/r;
		}	
	}
	return E;
}

double Energy(struct confor p)
{
	double Etotal;
	Etotal=EBond(p)+EAng(p)+EVpair(p)+CHpair(p);
	return Etotal;
}
