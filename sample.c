#include "global.h"
#include "sample.h"

#include	"mesch/matrix.h"

extern int ntotal;
extern double temp;	//temperature,Kelvin.	
extern double eta0;	//Solvent viscosity,poise;
extern double rm;	//Molecular weight 
extern double vbar;	//Partial specific volume of solute,g/cm3.
extern double solden;	//Solution density, cm3/g.
extern double Dt;	//translational diffusion coefficients.
extern double ft;	//translational friction coefficients.
extern double Dr;	//rotation diffusion coefficients.
extern double fr;	//rotation friction coefficients.

static MAT *I;
void GetI()
{
	static int is_malloc;	//we just need malloc one time.
	
	if(!is_malloc){
		I = m_get(3,3);
		int i,j;
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				if(i==j) I->me[i][j] = 1;
				else I->me[i][j] = 0;
			}
		}
		is_malloc = 1;
	}
}

double m_trace(MAT *D)
{
	double trace=0;
	int i,n;

	n = D->m;
	if(n != D->n) {
		printf("illegal MAT\n");
		exit(1);
	}

	for(i=0;i<n;i++) 	trace += D->me[i][i];

	return trace;
}

MAT *tensorT(struct confor *p,int i,int j)
/*Rotne-Pager-Yamakawa tensor*/
{
	static MAT *T,*RR,*I3;
	static int is_malloc;
	
	if(!is_malloc){
		T = m_get(3,3);
		RR = m_get(3,3);
		I3 = m_get(3,3); 
		is_malloc = 1;
	}
	
	
	double a,b,c,r;
	a = p->beads[j].x-p->beads[i].x;
	b = p->beads[j].y-p->beads[i].y;
	c = p->beads[j].z-p->beads[i].z;
	r = sqrt(a*a+b*b+c*c);

	RR->me[0][0] = a*a/r;
	RR->me[0][1] = a*b/r;
	RR->me[0][2] = a*c/r;
	RR->me[1][0] = b*a/r;
	RR->me[1][1] = b*b/r;
	RR->me[1][2] = b*c/r;
	RR->me[2][0] = c*a/r;
	RR->me[2][1] = c*b/r;
	RR->me[2][2] = c*c/r;

	sm_mlt(1/3.0,I,I3);
	m_sub(I3,RR,T);
	sm_mlt((p->beads[i].r*p->beads[i].r+p->beads[j].r*p->beads[j].r)/(r*r),T,T);
	m_add(RR,T,T);
	m_add(I,T,T);
	T = sm_mlt(1/(8*pi*eta0*r),T,T);
	 
	return T;
}

MAT *GetBigB(struct confor *p)
{
	int i,j,k,l;
	static MAT *BigB,*B;

	static int is_malloc;
	
	if(!is_malloc){
		B = m_get(3,3);
		BigB = m_get(3*ntotal,3*ntotal);
		is_malloc = 1;
	}

	for(i=0;i<ntotal;i++){
		for(j=0;j<ntotal;j++){
			if(i==j) sm_mlt((1/(6*pi*eta0*p->beads[i].r)),I,B);
			else B = tensorT(p,i,j);
			for(k=0;k<3;k++){
				for(l=0;l<3;l++){
					BigB->me[i*3+k][j*3+l] = B->me[k][l];
				}
			}
		}
	}
	
	return BigB;
}

MAT *GetBigE(struct confor *p)
{
	int i,j,k,l;
	static MAT *C,*Ett,*Etr,*Err,*BigB,*BigC,*BigE,*Ui,*Uj;

	static MAT *TEMP,*TEMP1;

	static int is_malloc;
	
	if(!is_malloc){
		TEMP = m_get(3,3);
		TEMP1 = m_get(3,3);
		C = m_get(3,3);
		Ett = m_get(3,3);
		Etr = m_get(3,3);
		Err = m_get(3,3);
		Ui = m_get(3,3);
		Uj = m_get(3,3);
		BigC = m_get(3*ntotal,3*ntotal);
		BigE = m_get(6,6);
		is_malloc = 1;
	}
	
	
	m_zero(Ett);m_zero(Etr);m_zero(Err);
	m_zero(Ui);m_zero(Uj);
	
	BigB = GetBigB(p);
	
	m_inverse(BigB,BigC);
	for(i=0;i<ntotal;i++){

		Ui->me[0][1] = -p->beads[i].z;
		Ui->me[0][2] = p->beads[i].y;
		Ui->me[1][0] = p->beads[i].z;
		Ui->me[1][2] = -p->beads[i].x;
		Ui->me[2][0] = -p->beads[i].y;
		Ui->me[2][1] = p->beads[i].x;
		
		for(j=0;j<ntotal;j++){
			for(k=0;k<3;k++){
				for(l=0;l<3;l++){
					C->me[k][l] = BigC->me[i*3+k][j*3+l];
				}
			}
			Uj->me[0][1] = -p->beads[j].z;
			Uj->me[0][2] = p->beads[j].y;
			Uj->me[1][0] = p->beads[j].z;
			Uj->me[1][2] = -p->beads[j].x;
			Uj->me[2][0] = -p->beads[j].y;
			Uj->me[2][1] = p->beads[j].x;
			
			m_add(C,Ett,Ett);

			m_mlt(Ui,C,TEMP);
			m_add(TEMP,Etr,Etr);

			m_mlt(Ui,C,TEMP);
			m_mlt(TEMP,Uj,TEMP1);
			m_add(TEMP1,Err,Err);
		}
	}

	for(k=0;k<3;k++){
		for(l=0;l<3;l++){
			BigE->me[k][l] = Ett->me[k][l];
		}
	}

	m_transp(Etr,TEMP);
	for(k=0;k<3;k++){
		for(l=0;l<3;l++){
			BigE->me[k][l+3] = TEMP->me[k][l];
		}
	}
	
	for(k=0;k<3;k++){
		for(l=0;l<3;l++){
			BigE->me[k+3][l] = Etr->me[k][l];
		}
	}
	for(k=0;k<3;k++){
		for(l=0;l<3;l++){
			BigE->me[k+3][l+3] = Err->me[k][l];
		}
	}

	//m_output(BigE);
	return BigE;
}

MAT *GetBigD(struct confor *p)
{
	static MAT *BigD;
	static int is_malloc;

	if(!is_malloc){
		BigD = m_get(3*ntotal,3*ntotal);
		is_malloc = 1;
	}
	
	sm_mlt(kB*temp,GetBigE(p),BigD);
	return BigD;
}

double add_average(double a,double sum,int n)
{
	return (1/((double)(n+1))*(n*sum+a));
}

void sample(struct confor *p,int n)
//n is MC steps when sampling
{

	GetI();	//Initial I
	
	static MAT *BigD;
	BigD = GetBigD(p);
	
	static MAT *Dtt,*Dtr,*Drr;
	static int is_malloc;

	if(!is_malloc){
		Dtt = m_get(3,3);
		Dtr = m_get(3,3);
		Drr = m_get(3,3);
		is_malloc = 1;
	}
	
	int k,l;
	for(k=0;k<3;k++){
		for(l=0;l<3;l++){
			Dtt->me[k][l] = BigD->me[k][l]; 
		}
	}
		
	for(k=0;k<3;k++){
		for(l=0;l<3;l++){
			Dtt->me[k][l] = BigD->me[k+3][l];
		}
	}
	for(k=0;k<3;k++){
		for(l=0;l<3;l++){
			Dtt->me[k][l] = BigD->me[k+3][l+3];
		}
	}

	double currDt = 1/3.0*m_trace(Dtt);	
	Dt = add_average(currDt,Dt,n);
	double currft = kB*temp/(currDt);
	ft = add_average(currft,ft,n);

	double currDr = 1/3.0*m_trace(Drr);	
	Dr = add_average(currDr,Dr,n);
	double currfr = kB*temp/(currDr);
	fr = add_average(currfr,fr,n);
		
}
