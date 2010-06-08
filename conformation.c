#include "global.h"
#include "conformation.h"
#include "random.h"

extern double r,R,theta,h;
extern int ntotal,nspring,nang;

extern double acceptrate; //accept rate

void GetSpring(confor *p)
//calculate spring properties.
//User supply code as your need.
{
	int i;
	for(i=0;i<(p->nsp/2);i++){
		p->springs[i].start = i;
		p->springs[i].end = i+1;
		p->springs[i].len=sqrt(
						pow((p->beads[i].x - p->beads[i+1].x),2)+
						pow((p->beads[i].y - p->beads[i+1].y),2)+
						pow((p->beads[i].z - p->beads[i+1].z),2));
	}
	for(i=nspring/2;i<nspring;i++){
		p->springs[i].start = i+1;
		p->springs[i].end = i+2;
		p->springs[i].len=sqrt(
						pow((p->beads[i+1].x - p->beads[i+2].x),2)+
						pow((p->beads[i+1].y - p->beads[i+2].y),2)+
						pow((p->beads[i+1].z - p->beads[i+2].z),2));
	}

}

void GetAngel(confor *p)
//calculate angs properties.
//User supply code as your need.
{
	int i;
	double a,b,c;	//three side of the triangle.
	for(i=0;i<nang/2;i++){
		p->angs[i].start = i;
		p->angs[i].end = i+1;
		a = p->springs[i].len;
		b = p->springs[i+1].len;
		c = sqrt(pow((p->beads[i].x-p->beads[i+2].x),2)+
				pow((p->beads[i].y-p->beads[i+2].y),2)+
				pow((p->beads[i].z-p->beads[i+2].z),2));
		p->angs[i].angle = acos((a*a+b*b-c*c)/(2*a*b));
	}
	for(i=nang/2;i<nang;i++){
		p->angs[i].start = i+1;
		p->angs[i].end = i+2;
		a = p->springs[i+1].len;
		b = p->springs[i+2].len;
		c = sqrt(pow((p->beads[i+2].x-p->beads[i+4].x),2)+
				pow((p->beads[i+2].y-p->beads[i+4].y),2)+
				pow((p->beads[i+2].z-p->beads[i+4].z),2));
		p->angs[i].angle = acos((a*a+b*b-c*c)/(2*a*b));
	}

}

void McMoveGlobal(const confor *old_p,confor *new_p)
{
	
	const double delalpha = 1;
	const double delbeta = 1;

	int o,i,nhalf;
	double alpha,Dalpha,beta,Dbeta,x,y,z,l;
	nhalf=ntotal/2;
	Dalpha=(rnd()-0.5)*delalpha;
	Dbeta=(rnd()-0.5)*delbeta;
	o = (int)(rnd()*(nhalf-2))+1;	//we don't want tail and nail beads.
	
	if((nhalf-o)>o){
		for(i=0;i<o;i++){
			x = old_p->beads[i].x - old_p->beads[o].x;
			y = old_p->beads[i].y - old_p->beads[o].y;
			z = old_p->beads[i].z - old_p->beads[o].z;
			l=sqrt(x*x+y*y+z*z);
			if(y<0) alpha=2*pi-acos(x/l);
			else alpha=acos(x/l);	
			if(z<0) beta=2*pi-acos(z/l);
			else beta=acos(z/l);

			alpha+=Dalpha;
			beta+=Dbeta;

			x=l*cos(alpha);
			y=l*sin(alpha);
			z=l*cos(beta);

			new_p->beads[i].x = old_p->beads[i].x+x;
			new_p->beads[i].y = old_p->beads[i].y+y;
			new_p->beads[i].z = old_p->beads[i].z+z;
		}
		for(i=nhalf;i<o+nhalf;i++){
			x = old_p->beads[i].x - old_p->beads[o+nhalf].x;
			y = old_p->beads[i].y - old_p->beads[o+nhalf].y;
			z = old_p->beads[i].z - old_p->beads[o+nhalf].z;
			l=sqrt(x*x+y*y+z*z);
			if(y<0) alpha=2*pi-acos(x/l);
			else alpha=acos(x/l);
			if(z<0) beta=2*pi-acos(z/l);
			else beta=acos(z/l);

			alpha+=Dalpha;
			beta+=Dbeta;

			x=l*cos(alpha);
			y=l*sin(alpha);
			z=l*cos(beta);

			new_p->beads[i].x = old_p->beads[o+nhalf].x+x;
			new_p->beads[i].y = old_p->beads[o+nhalf].y+y;
			new_p->beads[i].z = old_p->beads[o+nhalf].z+z;
		}
	} else {
		for(i=o+1;i<nhalf;i++){
			x = old_p->beads[i].x - old_p->beads[o].x;
			y = old_p->beads[i].y - old_p->beads[o].y;
			z = old_p->beads[i].z - old_p->beads[o].z;
			l=sqrt(x*x+y*y+z*z);
			if(y<0) alpha=2*pi-acos(x/l);
			else alpha=acos(x/l);
			if(z<0) beta=2*pi-acos(z/l);
			else beta=acos(z/l);

			alpha+=Dalpha;
			beta+=Dbeta;

			x=l*cos(alpha);
			y=l*sin(alpha);
			z=l*cos(beta);

			new_p->beads[i].x = old_p->beads[o].x+x;
			new_p->beads[i].y = old_p->beads[o].y+y;
			new_p->beads[i].z = old_p->beads[o].z+z;
		}
		for(i=o+nhalf;i<ntotal;i++){
			x = old_p->beads[i].x - old_p->beads[o+nhalf].x;
			y = old_p->beads[i].y - old_p->beads[o+nhalf].y;
			z = old_p->beads[i].z - old_p->beads[o+nhalf].z;
			l=sqrt(x*x+y*y+z*z);
			if(y<0) alpha=2*pi-acos(x/l);
			else alpha=acos(x/l);
			if(z<0) beta=2*pi-acos(z/l);
			else beta=acos(z/l);

			alpha+=Dalpha;
			beta+=Dbeta;

			x=l*cos(alpha);
			y=l*sin(alpha);
			z=l*cos(beta);

			new_p->beads[i].x = old_p->beads[o+nhalf].x+x;
			new_p->beads[i].y = old_p->beads[o+nhalf].y+y;
			new_p->beads[i].z = old_p->beads[o+nhalf].z+z;
		}
	}
	
	GetSpring(new_p);
	GetAngel(new_p);
	
}

void McMoveLocal(const confor *old_p,confor *new_p)
/*local single bead move*/
{
	static double delx = 0.5;	//max size of MC move
	static double dely = 0.5;
	static double delz = 0.5;
	int o;
	o = (int)(rnd()*ntotal); //select one particles at random
	
	//dynamic adjust monte carlo step length,according to accept rate.
	if(acceptrate>0.55){
		delx = delx*1.2;
		dely = dely*1.2;
		delz = delz*1.2;
	}
	if(acceptrate<0.3){
		delx = delx*0.8;
		dely = dely*0.8;
		delz = delz*0.8;
	}
	
	new_p->beads[o].x = old_p->beads[o].x + (rnd()-0.5)*delx;
	new_p->beads[o].y = old_p->beads[o].y + (rnd()-0.5)*dely;
	new_p->beads[o].z = old_p->beads[o].z + (rnd()-0.5)*delz;

	GetSpring(new_p);
	GetAngel(new_p);
	
}

void McMove(const confor *old_p,confor *new_p)
//A MC routine
{
	int i;
	i=rnd()*2;
	if(i==0) McMoveLocal(old_p,new_p);
	else McMoveLocal(old_p,new_p);
}

void InitialConfor(confor *p)
{
	//This is for the DNA
	double r,R,theta,h;
	r=3;R=10;theta=36;h=3.4;

	int i,j;
	
	for(i=0;i<ntotal/2;i++)
	{
		p->beads[i].x=R*cos(theta*pi*i/180);
		p->beads[i].y=R*sin(theta*pi*i/180);
		p->beads[i].z=h*i;
		p->beads[i].r=r;
	}
	for(i=ntotal/2;i<ntotal;i++)
	{
		j=i-ntotal/2;
		p->beads[i].x=R*cos(theta*j*pi/180+180);
		p->beads[i].y=R*sin(theta*j*pi/180+180);
		p->beads[i].z=h*j;
		p->beads[i].r=r;
	}
	
	GetSpring(p);
	GetAngel(p);

}




