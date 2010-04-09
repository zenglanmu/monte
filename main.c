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
int nconf;	//max number of conformation.

//user code,model data.
//This is for the DNA
int n;
double r,R,theta,h;

void userdata()
{
	temp=293;
	eta0=0.010;
	rm=110000;
	vbar=0.710;
	solden=1.0;
	title="The DNA"	; 
	filename="confi";
	nconf=36;	//max number of conformation.
	n=200;
	r=3;R=10;theta=36;h=3.4;
}

int main(int argc, char** argv)
{
	int i,j;
	userdata();
	struct confor conf;
	conf=initconfor(r,R,theta,h,n);
	writepdb(conf,n,"dnaint.pdb");
		//Check overlap between spheres
		int flag;
		for(j=0;j<n;j++)
		{
			flag=ls_overlap(conf,j,n);
			if(flag)
			{
				printf("overlap happens at conf.beads[\%d]\n",j);
				break;
			}
		}
		
	
	return 0;
}
