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
#include "io.h"

//user specific code,var for run.
double temp;	//temperature,Kelvin.	
double eta0;	//Solvent viscosity,poise;
double rm;	//Molecular weight 
double vbar;	//Partial specific volume of solute,g/cm3.
double solden;	//Solution density, cm3/g.
double Dt;	//translational diffusion coefficients.
double ft;	//translational friction coefficients.
double Dr;	//rotation diffusion coefficients.
double fr;	//rotation friction coefficients.
double tao[5] = {0,0,0,0,0};	//Relaxation time.
double taoh;	//Harmonic mean (correlation) time.
double Rg;	//Radius of gyration.
double eta;	//Intrinsic viscosity.

char *title;	
char *filename;
int nstep;	//times of MC steps.
int nreject; //reject times;
double acceptrate; //accept rate

//user code,model data.
int ntotal,nspring,nang;


void UserData()
{
	temp=293;
	eta0=0.010;
	rm=110000;
	vbar=0.710;
	solden=1.0;
	Dt=0;ft=0;
	Dr=0;fr=0;
	tao[5]=(0,0,0,0,0);
	taoh =0;
	Rg = 0;
	eta = 0;

	
	title="The DNA"	; 
	filename="DNA-results.txt";
	
	nstep=1000;
	nreject=0;
	acceptrate=0.5;
	
	ntotal=50;
	nspring=ntotal-2;
	nang=nspring-2;
	
}

confor *confor_get(int nbd,int nsp,int nang)
//malloc a confor which have nbd beads
//and nsp springs.
{
	confor *p;

	p = malloc(sizeof(confor));
	if (p == NULL) {
		printf("out of memory\n");
		exit(1);
	}
	p->nbd = nbd;
	p->nsp = nsp;
	p->nang = nang;
	p->beads = malloc(sizeof(bead)*nbd);
	p->springs = malloc(sizeof(spring)*nsp);
	p->angs = malloc(sizeof(ang)*nang);

	return p;
}

void confor_free(confor *p)
//free a confor.
{
	free(p->beads);
	free(p->springs);
	free(p->angs);
	free(p);
}

void confor_copy(confor *from,confor *to)
//copy a confor to another confor.
//both confors need to malloc first.
{
	int i;
	if((from->nbd != to->nbd) || (from->nsp != to->nsp) || (from->nang != to->nang)){
		printf("confors don't have equal size,can't copy!\n");
		exit(1);
	}
	for(i=0;i < from->nbd;i++){
		to->beads[i].x = from->beads[i].x;
		to->beads[i].y = from->beads[i].y;
		to->beads[i].z = from->beads[i].z;
		to->beads[i].r = from->beads[i].r;
	}
	for(i=0;i < from->nsp;i++){
		to->springs[i].start = from->springs[i].start;
		to->springs[i].end = from->springs[i].end;
		to->springs[i].len = from->springs[i].len;
	}
	for(i=0;i < from->nang;i++){
		to->angs[i].start = from->angs[i].start;
		to->angs[i].end = from->angs[i].end;
		to->angs[i].angle = from->angs[i].angle;
	}
	
}

void RejectConfor(char *s)
{
	printf("Conformation rejected,reason:");
	printf("%s\n",s);
	nreject++;
}

void AcceptConfor(confor *newconf,confor *conf)
{
	confor_copy(newconf,conf);	//Let old comformation = new comformation 
}
	
int main(int argc, char** argv)
{
	int i;
	double Eprev,Enew,u;
	confor *conf,*newconf;

	UserData();
	
	conf = confor_get(ntotal,nspring,nang);
	newconf = confor_get(ntotal,nspring,nang);
	
	InitialConfor(conf);
	confor_copy(conf,newconf);
	
	SavePDBFile(conf,"DNAInitial.pdb");
	
	for(i=0;i<=nstep;i++){
		printf("run nstep times %d\n",i);
		
		//sample every 100 run.
		//don't sample too often,or you would wait too long till stop.
		if(!(i%100)){
			sample(conf,i);
			results_output(stdout);
		}
		
		acceptrate = (i-nreject)/(double)i;	
		
		Eprev = Energy(conf);
		McMove(conf,newconf);
		
		if(ls_overlap(newconf)){	
			RejectConfor("overlap");
			continue;
		}
		
		Enew = Energy(newconf);
		if(Enew <= Eprev){	//accept conformation.
			AcceptConfor(newconf,conf);
		} else {
			u = rnd();
			if(u<exp(-(Enew-Eprev)/(kB*temp))) 
				AcceptConfor(newconf,conf);
			else{
				RejectConfor("Energy too high");
				continue;
			}
		}
		
		
	}

	SavePDBFile(conf,"DNAFinal.pdb");
	printf("\naccept rate %f\n",acceptrate);

	/*write results*/
	results_output(stdout);
	FILE *fp;
	fp=fopen(filename,"w");
	results_output(fp);
	fclose(fp);
	return 0;
}
