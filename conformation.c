#include "conformation.h"
#include "global.h"
#include "random.h"

extern double r,R,theta,h;
extern int ntotal,nspring,nang;
extern double Edist[nmax],Eang[nmax];	//equilibrium spring lens and equilibrium angel.

struct confor GetSpring(struct confor p)
//calculate spring properties.
{
	int i,j;
	for(i=0;i<nspring/2;i++){
		p.springs[i].x0=p.beads[i].x;
		p.springs[i].y0=p.beads[i].y;
		p.springs[i].z0=p.beads[i].z;
		p.springs[i].x1=p.beads[i+1].x;
		p.springs[i].y1=p.beads[i+1].y;
		p.springs[i].z1=p.beads[i+1].z;
		p.springs[i].len=sqrt(pow(((p.springs[i].x0)-(p.springs[i].x1)),2)+pow(((p.springs[i].y0)-(p.springs[i].y1)),2)+pow(((p.springs[i].z0)-(p.springs[i].z1)),2));
	}
	for(i=nspring/2;i<nspring;i++){
		j=i+1;
		p.springs[i].x0=p.beads[j].x;
		p.springs[i].y0=p.beads[j].y;
		p.springs[i].z0=p.beads[j].z;
		p.springs[i].x1=p.beads[j+1].x;
		p.springs[i].y1=p.beads[j+1].y;
		p.springs[i].z1=p.beads[j+1].z;
		p.springs[i].len=sqrt(pow(((p.springs[i].x0)-(p.springs[i].x1)),2)+pow(((p.springs[i].y0)-(p.springs[i].y1)),2)+pow(((p.springs[i].z0)-(p.springs[i].z1)),2));
	}
	return p;
}

void GetAngel(struct confor p,double *ang)
{
	int i,j;
	double a,b,c;	//three side of the triangle.
	for(i=0;i<nang/2;i++){
		a=p.springs[i].len;
		b=p.springs[i+1].len;
		c=sqrt(pow(((p.springs[i].x0)-(p.springs[i+1].x1)),2)+pow(((p.springs[i].y0)-(p.springs[i+1].y1)),2)+pow(((p.springs[i].z0)-(p.springs[i+1].z1)),2));
		*ang=(a*a+b*b-c*c)/(2*a*b);
		ang++;
	}
	for(i=nang/2;i<nang;i++){
		j=i+1;
		a=p.springs[j].len;
		b=p.springs[j+1].len;
		c=sqrt(pow(((p.springs[j].x0)-(p.springs[j+1].x1)),2)+pow(((p.springs[j].y0)-(p.springs[j+1].y1)),2)+pow(((p.springs[j].z0)-(p.springs[j+1].z1)),2));
		*ang=(a*a+b*b-c*c)/(2*a*b);
		ang++;
	}
}

struct confor McMoveGlobal(struct confor p)
{
	int o,i,nhalf;
	double alpha,Dalpha,beta,Dbeta,x,y,z,l;
	nhalf=ntotal/2;
	Dalpha=rand()*delalpha;
	Dbeta=rand()*delbeta;
	o = (int)(rand()*(nhalf-2))+1;	//we don't want tail and nail beads.

	if((nhalf-o)>o){
		for(i=0;i<o;i++){
			x=p.beads[i].x-p.beads[o].x;
			y=p.beads[i].y-p.beads[o].y;
			z=p.beads[i].z-p.beads[o].z;
			l=sqrt(x*x+y*y+z*z);
			if(y<0) alpha=2*pi-acos(x/l);
			else alpha=acos(x/l);	
			beta=acos(z/l);

			alpha+=Dalpha;
			beta+=Dbeta;

			x=l*cos(alpha);
			y=l*sin(alpha);
			z=l*cos(beta);

			p.beads[i].x=p.beads[o].x+x;
			p.beads[i].y=p.beads[o].y+y;
			p.beads[i].z=p.beads[o].z+z;
		}
		for(i=nhalf;i<o+nhalf;i++){
			x=p.beads[i].x-p.beads[o+nhalf].x;
			y=p.beads[i].y-p.beads[o+nhalf].y;
			z=p.beads[i].z-p.beads[o+nhalf].z;
			l=sqrt(x*x+y*y+z*z);
			if(y<0) alpha=2*pi-acos(x/l);
			else alpha=acos(x/l);
			beta=acos(z/l);

			alpha+=Dalpha;
			beta+=Dbeta;

			x=l*cos(alpha);
			y=l*sin(alpha);
			z=l*cos(beta);

			p.beads[i].x=p.beads[o+nhalf].x+x;
			p.beads[i].y=p.beads[o+nhalf].y+y;
			p.beads[i].z=p.beads[o+nhalf].z+z;
		}
	} else {
		for(i=o+1;i<nhalf;i++){
			x=p.beads[i].x-p.beads[o].x;
			y=p.beads[i].y-p.beads[o].y;
			z=p.beads[i].z-p.beads[o].z;
			l=sqrt(x*x+y*y+z*z);
			if(y<0) alpha=2*pi-acos(x/l);
			else alpha=acos(x/l);
			beta=acos(z/l);

			alpha+=Dalpha;
			beta+=Dbeta;

			x=l*cos(alpha);
			y=l*sin(alpha);
			z=l*cos(beta);

			p.beads[i].x=p.beads[o].x+x;
			p.beads[i].y=p.beads[o].y+y;
			p.beads[i].z=p.beads[o].z+z;
		}
		for(i=o+nhalf;i<ntotal;i++){
			x=p.beads[i].x-p.beads[o+nhalf].x;
			y=p.beads[i].y-p.beads[o+nhalf].y;
			z=p.beads[i].z-p.beads[o+nhalf].z;
			l=sqrt(x*x+y*y+z*z);
			if(y<0) alpha=2*pi-acos(x/l);
			else alpha=acos(x/l);
			beta=acos(z/l);

			alpha+=Dalpha;
			beta+=Dbeta;

			x=l*cos(alpha);
			y=l*sin(alpha);
			z=l*cos(beta);

			p.beads[i].x=p.beads[o+nhalf].x+x;
			p.beads[i].y=p.beads[o+nhalf].y+y;
			p.beads[i].z=p.beads[o+nhalf].z+z;
		}
	}
	p = GetSpring(p);
	return p;
}

struct confor McMoveLocal(struct confor p)
{
	int o;
	o = (int)(rand()*ntotal); //select one particles at random
	p.beads[o].x = p.beads[o].x+(rand()-0.5)*delx;
	p.beads[o].y = p.beads[o].y+(rand()-0.5)*dely;
	p.beads[o].z = p.beads[o].z+(rand()-0.5)*delz;
	p = GetSpring(p);
	return p; 
}

struct confor McMove(struct confor p)
//A MC routine
{
	int i;
	i=rand()*2;
	if(i==0) p=McMoveLocal(p);
	else p=McMoveGlobal(p);
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
	
	p = GetSpring(p);

	for(i=0;i<nspring;i++) Edist[i]=p.springs[i].len;
	GetAngel(p,Eang);
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


