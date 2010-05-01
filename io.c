#include "global.h"
#include "io.h"

void SavePDBFile(confor *p,char *filename)
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

void confor_output(confor *p)
{
	int i;
	printf("given confor have %d beads:\n",p->nbd);
	for(i=0;i < p->nbd;i++)
		printf("beads[%d]:x = %-6.4f,y = %-6.4f,z = %-6.4f,r = %-6.4f\n",
			i,p->beads[i].x,p->beads[i].y,p->beads[i].z,p->beads[i].r);

	printf("given confor have %d springs\n",p->nsp);
	for(i=0;i < p->nsp ;i++)
		printf("springs[%d]:len = %-6.4f\n",i,p->springs[i].len);

	printf("given confor have %d angs\n",p->nang);
	for(i=0;i < p->nang;i++)
		printf("angs[%d]:ang = %-6.4f\n",i,p->angs[i].angle);		
	
}
