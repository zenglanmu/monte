#include <string.h>
#include "global.h"
#include "io.h"
#include "libpdb/pdb.h"

extern double temp;	//temperature,Kelvin.
extern double Dt;	//translational diffusion coefficients.
extern double Dr;	//rotation diffusion coefficients.
extern double tao[5];	//Relaxation time.
extern double taoh;	//Harmonic mean (correlation) time.
extern double Rg;	//Radius of gyration.
extern double eta;	//Intrinsic viscosity.


void SavePDBFile(const confor *p,char *filename)
{
	int i;
	pdb_record r;
	pdb_residue s;
	FILE *fp;
	
	fp = fopen(filename,"w");
	r.record_type = PDB_ATOM;
	for(i=0;i< p->nbd ;i++)
	{
		r.pdb.atom.serial_num = i;
		strcpy(r.pdb.atom.name,"C");
		r.pdb.atom.alt_loc = 0;
		
		strcpy(s.name,"C");
		s.chain_id = 'A';
		s.seq_num = i;
		s.insert_code = 'A';
		r.pdb.atom.residue = s;
		
		r.pdb.atom.x = p->beads[i].x;
		r.pdb.atom.y = p->beads[i].y;
		r.pdb.atom.z = p->beads[i].z;
		r.pdb.atom.occupancy = p->beads[i].r;
		
		r.pdb.atom.temp_factor = temp;
		pdb_write_record(fp,&r,NULL,i);
	}
	fclose(fp);
}

void ReadPDBFILE(confor *p,char *filename)
{
	int i=0;
	pdb_record r;
	FILE *fp;
	fp = fopen(filename,"r");
	
	while(i<p->nbd){
		r = pdb_read_record(fp);
		if(r.record_type == PDB_ATOM){
			p->beads[i].x = r.pdb.atom.x;
			p->beads[i].y = r.pdb.atom.y;
			p->beads[i].z = r.pdb.atom.z;
			p->beads[i].r = r.pdb.atom.occupancy;
			i++;
		}
	}
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
	
	fprintf(fp,"Harmonic mean (correlation) time: %-6.4e s\n",taoh);	
	fprintf(fp,"Radius of gyration: %-6.4e cm\n",Rg);	
	fprintf(fp,"Intrinsic viscosity: %-6.4e cm^3/g\n",eta);
}
