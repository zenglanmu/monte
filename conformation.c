#include "conformation.h"
#include "global.h"
#include "random.h"

extern double r,R,theta,h;
extern int ntotal;

struct confor McMove(struct confor p)
//A MC routine
{
	int o;
	o = (int)(rand()*ntotal); //select one particles at random
	p.beads[o].x = p.beads[o].x+(rand()-0.5)*delx;
	p.beads[o].y = p.beads[o].y+(rand()-0.5)*dely;
	p.beads[o].z = p.beads[o].z+(rand()-0.5)*delz;
	return p; 
}

struct confor InitialConfor()
{
	struct confor p;
	int i,j;
	for(i=0;i<ntotal/2;i++)
	{
		p.beads[i].x=R*cos(theta*pi*i/180);
		p.beads[i].y=R*sin(theta*pi*i/180);
		p.beads[i].z=h*i;
		p.beads[i].r=r;
	}
	for(i=ntotal/2;i<ntotal;i++)
	{
		j=i-ntotal/2;
		p.beads[i].x=R*cos(theta*j*pi/180+180);
		p.beads[i].y=R*sin(theta*j*pi/180+180);
		p.beads[i].z=h*j;
		p.beads[i].r=r;
	}
	return p;
}

void SavePDBFile(struct confor p,char *filename)
{
	int i;
	char s[5];
	FILE *fp;
	fp=fopen(filename,"w");
	fprintf(fp,"REMARK   PDB file created by monte\n");
	for(i=0;i<ntotal;i++)
	{
		fprintf(fp,"ATOM  ");
		sprintf(s,"C%d",i+1);
		fprintf(fp,"%5d %-4.4s%c%3.3s %c%4d    ",i+1,"C"," "," ",'A',0);
		fprintf(fp,"%8.3f%8.3f%8.3f",p.beads[i].x,-p.beads[i].y,-p.beads[i].z);
		fprintf(fp,"  1.00%6.2f\n",0.0);
	}
	fprintf(fp,"END   \n");
	fclose(fp);
}


