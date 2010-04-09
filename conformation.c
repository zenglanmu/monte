#include "conformation.h"
#include "global.h"

//User code.
//Construction of each conformation
struct confor mcmove(struct confor p)
{
	return p; 
}

struct confor initconfor(double r,double R,double theta,double h,double N)
{
	struct confor p;
	int i,j;
	for(i=0;i<N/2;i++)
	{
		p.beads[i].x=R*cos(theta*pi*i/180);
		p.beads[i].y=R*sin(theta*pi*i/180);
		p.beads[i].z=h*i;
		p.beads[i].r=r;
	}
	for(i=N/2;i<N;i++)
	{
		j=i-N/2;
		p.beads[i].x=R*cos(theta*j*pi/180+180);
		p.beads[i].y=R*sin(theta*j*pi/180+180);
		p.beads[i].z=h*j;
		p.beads[i].r=r;
	}
	return p;
}

void writepdb(struct confor p,int N,char *filename)
{
	int i;
	char s[5];
	FILE *fp;
	fp=fopen(filename,"w");
	fprintf(fp,"REMARK   PDB file created by monte\n");
	for(i=0;i<N;i++)
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


