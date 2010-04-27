//      multi.c
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      



#include "global.h"
#include "conformation.h"
#include "overlap.h"
#include "energy.h"
#include "random.h"
#include "sample.h"

//user specific code,var for run.
double temp;	//temperature,Kelvin.	
double eta0;	//Solvent viscosity,poise;
double rm;	//Molecular weight 
double vbar;	//Partial specific volume of solute,g/cm3.
double solden;	//Solution density, cm3/g. 
char *title;	//limit 20 chars 
char *filename;
int nstep;	//times of MC steps.
int nreject;

//user code,model data.
//This is for the DNA
int ntotal,nspring,nang;
double r,R,theta,h;
double Edist[nmax],Eang[nmax];	//equilibrium spring lens and equilibrium angel.
struct confor InitialConf;

void UserData()
{
	temp=293;
	eta0=0.0010;
	rm=110000;
	vbar=0.710;
	solden=1.0;
	title="The DNA"	; 
	filename="DNA.pdb";
	nstep=100;
	nreject=0;	
	ntotal=100;
	nspring=ntotal-2;
	nang=nspring-2;
	r=3;R=10;theta=36;h=3.4;
}

void RejectConfor(char *s)
{
	printf("Conformation rejected,reason:");
	printf("%s\n",s);
	nreject++;
}

void AcceptConfor(struct confor *newconf,struct confor *conf)
{
	sample(*conf);
	*conf=*newconf;		//Let old comformation = new comformation 

}
	
int main(int argc, char** argv)
{
	int i;
	double Eprev,Enew,u;
	struct confor conf,newconf;

	UserData();
	
	conf=InitialConf=InitialConfor();
	for(i=0;i<nspring;i++) Edist[i]=InitialConf.springs[i].len;
	GetAngel(InitialConf,Eang);
	
	SavePDBFile(conf,"DNAInitial.pdb");
	
	for(i=0;i<nstep;i++){
		printf("run nstep times %d\n",i);
		Eprev = Energy(conf);
		newconf = McMove(conf);
		if(ls_overlap(newconf,ntotal)){	/*ls_overlap(newconf,ntotal)*/
			RejectConfor("overlap");
			continue;
			}
		Enew = Energy(newconf);
		if(Enew<Eprev){	//accept conformation.
			AcceptConfor(&newconf,&conf);
		} else {
			u = rand();
			if(u<exp(-(Enew-Eprev)/(kB*temp))) /*u<exp(-(Enew-Eprev)/(kB*temp))*/
				AcceptConfor(&newconf,&conf);
			else{
				RejectConfor("Energy too high");
				continue;
			}
		}	
	}

	SavePDBFile(conf,"DNAFinal.pdb");
	printf("\naccept rate %f\n",(nstep-nreject)/(double)nstep);
	return 0;
}
