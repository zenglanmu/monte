#include "global.h"
#include "sample.h"

#include	"mesch/matrix.h"

extern int ntotal;
extern double temp;	//temperature,Kelvin.	
extern double eta0;	//Solvent viscosity,poise;
extern double rm;	//Molecular weight 
extern double vbar;	//Partial specific volume of solute,g/cm3.
extern double solden;	//Solution density, cm3/g.
extern double Dt;
extern double ft;
extern double Dr;
extern double fr;

static MAT *I;
void GetI()
{
	I = m_get(3,3);
	int i,j;
   	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			if(i==j) I->me[i][j] = 1;
			else I->me[i][j] = 0;
		}
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

MAT *tensorT(int i,int j)
{
	static MAT *T;
	T = m_get(3,3);
	T = m_copy(I,T);
	return T;
}

MAT *GetBigB(struct confor *p)
{
	int i,j,k,l;
	static MAT *BigB,*B;
	
	B = m_get(3,3);
	BigB = m_get(3*ntotal,3*ntotal);
	
	for(i=0;i<ntotal;i++){
		for(j=0;j<ntotal;j++){
			if(i==j) sm_mlt((1/(6*pi*eta0*p->beads[i].r)),I,B);
			else B = tensorT(i,j);
			for(k=0;k<3;k++){
				for(l=0;l<3;l++){
					BigB->me[i+k][j+l] = B->me[k][l];
				}
			}
		}
	}
	
	M_FREE(B);
	return BigB;
}

MAT *GetBigE(struct confor *p)
{
	int i,j,k,l;
	static MAT *C,*Ett,*Etr,*Err,*BigB,*BigC,*BigE,*Ui,*Uj;

	static MAT *TEMP;
	TEMP = m_get(3,3);
	
	C = m_get(3,3);
	Ett = m_get(3,3);Etr = m_get(3,3);Err = m_get(3,3);
	Ui = m_get(3,3);Uj = m_get(3,3);
	m_zero(Ett);m_zero(Etr);m_zero(Err);
	m_zero(Ui);m_zero(Uj);
	
	BigC = m_get(3*ntotal,3*ntotal);
	BigE = m_get(6,6);
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
					C->me[k][l] = BigC->me[i+k][j+l];
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
			m_mlt(TEMP,Uj,TEMP);
			m_add(TEMP,Err,Err);
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

	m_output(BigE);
	return BigE;
}

MAT *GetBigD(struct confor *p)
{
	static MAT *BigD;
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
	Dtt = m_get(3,3);Dtr = m_get(3,3);Drr = m_get(3,3);

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
