#include "global.h"
#include "io.h"

extern double Dt;	//translational diffusion coefficients.
extern double Dr;	//rotation diffusion coefficients.
extern double tao[5];	//Relaxation time.

void SavePDBFile(const confor *p,char *filename)
{
	int i;
	char s[5];
	FILE *fp;
	fp=fopen(filename,"w");
	fprintf(fp,"REMARK   PDB file created by monte\n");
	for(i=0;i< p->nbd ;i++)
	{
		fprintf(fp,"ATOM  ");
		sprintf(s,"C%d",i+1);
		fprintf(fp,"%5d %-4.4s%c%3.3s %c%4d    ",i+1,"C"," "," ",'A',0);
		fprintf(fp,"%8.3f%8.3f%8.3f",p->beads[i].x,-p->beads[i].y,-p->beads[i].z);
		fprintf(fp,"  1.00%6.2f\n",0.0);
	}
	fprintf(fp,"END   \n");
	fclose(fp);
}

void confor_foutput(FILE *fp, const confor *p)
{
	unsigned int i;

	if ( p == (confor *)NULL ){
		fprintf(fp,"confor: NULL\n");
		return;
	}
		
    fprintf(fp,"confor: beads: %d springs: %d angs: %d  \n"
					,p->nbd,p->nsp,p->nang);
    if ( p->beads == (bead *)NULL ){
		fprintf(fp,"NULL\n");
		return;
	}
     
	fprintf(fp,"given confor have %d beads:\n",p->nbd);
	for(i=0;i < p->nbd;i++)
		fprintf(fp,"beads[%d]:x = %-6.4f,y = %-6.4f,z = %-6.4f,r = %-6.4f\n",
			i,p->beads[i].x,p->beads[i].y,p->beads[i].z,p->beads[i].r);

	fprintf(fp,"given confor have %d springs\n",p->nsp);
	for(i=0;i < p->nsp ;i++)
		fprintf(fp,"springs[%d]:len = %-6.4f\n",i,p->springs[i].len);

	fprintf(fp,"given confor have %d angs\n",p->nang);
	for(i=0;i < p->nang;i++)
		fprintf(fp,"angs[%d]:ang = %-6.4f\n",i,p->angs[i].angle);		
	
}

void results_output(FILE *fp)
{
	char *intro="-------------------------------------------------\nMONTE by zenglanmu@126.com\n-------------------------------------------------\n";

	fprintf(fp,"%s",intro);			
	fprintf(fp,"\nTranslational diffusion coefficient: %-6.4e cm2/s\n",Dt);
	fprintf(fp,"\nRotational diffusion coefficient: %-6.4e s-1\n",Dr);

	int i;
	for(i = 0;i<5;i++)
		fprintf(fp,"Relaxation time (%d): %-6.4e s\n",i+1,tao[i]);
}
