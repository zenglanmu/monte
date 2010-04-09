//      multi.c
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      



#include "conformation.h"
#include "overlap.h"
#include "global.h"

//user specific code,var for run.
//double temp,eta0,rm,vbar,solden;
double temp,eta0 ,rm ,vbar, solden;
char *title;	//limit 20 chars 
char *filename;
int nstep;	//times of MC steps.

//user code,model data.
//This is for the DNA
int ntotal;
double r,R,theta,h;

void UserData()
{
	temp=293;
	eta0=0.010;
	rm=110000;
	vbar=0.710;
	solden=1.0;
	title="The DNA"	; 
	filename="DNA.pdb";
	nstep=10000;	
	ntotal=200;
	r=3;R=10;theta=36;h=3.4;
}

int main(int argc, char** argv)
{
	int i,flag;
	UserData();
	struct confor conf;
	conf=InitialConfor();
	for(i=0;i<nstep;i++) conf = McMove(conf);
	SavePDBFile(conf,filename);
	flag=ls_overlap(conf,ntotal);
	if(flag) printf("overlap happens\n");	
	
	return 0;
}
