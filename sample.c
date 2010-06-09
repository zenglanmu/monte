#include "mesch/matrix.h"
#include "mesch/matrix2.h"
#include "global.h"
#include "sample.h"



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
extern double tao[5];	//Relaxation time.
extern double Rg;	//Radius of gyration.
extern double taoh;	//Harmonic mean (correlation) time.
extern double eta;	//Intrinsic viscosity.

MAT *BigC; //make global here to caculate Intrinsic viscosity.

static MAT *I;
void GetI()
{
	static int is_malloc = 0;	//we just need malloc one time.
	
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

static MAT *etaD;
void GetetaD()
{
	static int is_malloc = 0;	//we just need malloc one time.
	
	if(!is_malloc){
		etaD = m_get(3,3);
		int i,j;
		m_zero(etaD);
		etaD->me[0][1] = 1;
		is_malloc = 1;
	}
}

static MAT *etaE;
void GetetaE()
{
	static int is_malloc = 0;	//we just need malloc one time.
	
	if(!is_malloc){
		etaE = m_get(3,3);
		int i,j;
		m_zero(etaE);
		etaE->me[0][1] = 1;
		etaE->me[1][0] = 1;
		is_malloc = 1;
	}
}

double m_trace(const MAT *D)
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

MAT *tensorT(const confor *p,int i,int j)
/*Rotne-Pager-Yamakawa tensor*/
{
	static MAT *T,*RR,*I3;
	static int is_malloc = 0;
	
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

MAT *GetBigB(const confor *p)
{
	int i,j,k,l;
	static MAT *BigB,*B;

	static int is_malloc = 0;
	
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

MAT *GetBigE(const confor *p)
{
	int i,j,k,l;
	static MAT *C,*Ett,*Etr,*Err,*BigB,*BigE,*Ui,*Uj;

	static MAT *TEMP,*TEMP1;

	static int is_malloc = 0;
	
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
	
	
	m_zero(Ett);
	m_zero(Etr);
	m_zero(Err);
	m_zero(Ui);
	m_zero(Uj);
	
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

	/*we need some unit converts here*/
	/*so we have SI units*/
	sm_mlt(10E-10,Ett,Ett);
	sm_mlt(10E-19,Etr,Etr);
	sm_mlt(10E-28,Err,Err);
	
	//Volume correction 
	double Volume;
	Volume = 0;
	for(i=0;i<p->nbd;i++)
		Volume += 4/3.0*pi*pow(p->beads[i].r,3);
		
	sm_mlt(6*eta0*Volume,I,TEMP);
	sm_mlt(10E-28,TEMP,TEMP);
	m_add(TEMP,Err,Err);
	
	
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

	static MAT *Dtt;
	static int is_malloc2 = 0;

	if(!is_malloc2){
		Dtt = m_get(3,3);
		is_malloc2 = 1;
	}
	m_inverse(Ett,Dtt);
	sm_mlt(kB*temp,Dtt,Dtt);
	
	return BigE;
}

MAT *GetBigD(const confor *p)
{
	static MAT *BigD,*TEMP;
	static int is_malloc = 0;

	if(!is_malloc){
		BigD = m_get(6,6);
		TEMP = m_get(6,6);
		is_malloc = 1;
	}

	m_inverse(GetBigE(p),TEMP); 
	sm_mlt(kB*temp,TEMP,BigD);
	return BigD;
}

double add_average(double a,double sum,int n)
{
	return (1/((double)(n+1))*(n*sum+a));
}

void sample(const confor *p,int n)
//n is MC steps when sampling
{

	GetI();	//Initial I
	GetetaD();	//Initial etaD
	GetetaE();	//Initial etaE
	
	static MAT *BigD;
	BigD = GetBigD(p);

	printf("Matrix BigD:\n");
	m_output(BigD);
	
	static MAT *Dtt,*Dtr,*Drr;
	static int is_malloc = 0;

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
			Dtr->me[k][l] = BigD->me[k+3][l];
		}
	}
	for(k=0;k<3;k++){
		for(l=0;l<3;l++){
			Drr->me[k][l] = BigD->me[k+3][l+3];
		}
	}

	/*we need some unit converts here*/
	/*we convert m to cm*/
	sm_mlt(10E4,Dtt,Dtt);
	sm_mlt(10E2,Dtr,Dtr);
	
	double currDt = 1/3.0*m_trace(Dtt);
	currDt = fabs(currDt);
	Dt = add_average(currDt,Dt,n);
	double currft = kB*temp/(currDt);
	ft = add_average(currft,ft,n);

	double currDr = 1/3.0*m_trace(Drr);
	currDr = fabs(currDr);
	Dr = add_average(currDr,Dr,n);
	double currfr = kB*temp/(currDr);
	fr = add_average(currfr,fr,n);

	/*rotation relaxation time */
	/*we need Drr's EIGENVALUES at first*/
	static MAT *TT, *QQ;
	static VEC *evals_re, *evals_im;
	static int is_malloc1 = 0;

	if(!is_malloc1){
		QQ = m_get(Drr->m,Drr->n);
		TT = m_get(Drr->m,Drr->n);
		evals_re = v_get(Drr->m);
		evals_im = v_get(Drr->m);
		is_malloc = 1;
	}

	TT = m_copy(Drr,TT);

	/* compute Schur form: Drr = QQ*TT*QQ^TT */
	schur(TT,QQ);
 
	/* extract eigenvalues */
	schur_evals(TT,evals_re,evals_im);
	/*We have gotten Drr's EIGENVALUES in evals_re*/

  	double lamda1,lamda2,lamda3,delta;
  	double currtao[5];
  	lamda1 = evals_re->ve[0];
  	lamda2 = evals_re->ve[1];
  	lamda3 = evals_re->ve[2];
  	delta = sqrt(lamda1*lamda1+lamda2*lamda2+lamda3*lamda3-
		lamda1*lamda2-lamda1*lamda3-lamda2*lamda3);
	currtao[0] = 1/(6*Dr - 2*delta);
	currtao[0] = fabs(currtao[0]);
	tao[0] = add_average(currtao[0],tao[0],n);
	
	currtao[1] = 1/(3*(Dr + lamda1));
	currtao[1] = fabs(currtao[1]);
	tao[1] = add_average(currtao[1],tao[1],n);

	currtao[2] = 1/(3*(Dr + lamda2));
	currtao[2] = fabs(currtao[2]);
	tao[2] = add_average(currtao[2],tao[2],n);

	currtao[3] = 1/(3*(Dr + lamda3));
	currtao[3] = fabs(currtao[3]);
	tao[3] = add_average(currtao[3],tao[3],n);

	currtao[4] = 1/(6*Dr + 2*delta);
	currtao[4] = fabs(currtao[4]);
	tao[4] = add_average(currtao[4],tao[4],n);
	
	double currtaoh = 0; 
	int i,j;
	for(i = 0;i<5;i++) currtaoh += 1/currtao[i];
	currtaoh = 1/currtaoh;
	taoh = add_average(currtaoh,taoh,n);
		
	//Radius of gyration:
	double currRg = 0;
	double fi,fj;
	double a,b,c,r; 
	double RadiusSum=0;
	for(i=0;i<ntotal;i++) RadiusSum += pow(p->beads[i].r,3);
	for(i=0;i<ntotal;i++){
		fi = pow(p->beads[i].r,3)/RadiusSum;
		for(j=0;j<ntotal;j++){
			fj = pow(p->beads[j].r,3)/RadiusSum;
			a = p->beads[j].x-p->beads[i].x;
			b = p->beads[j].y-p->beads[i].y;
			c = p->beads[j].z-p->beads[i].z;
			r = sqrt(a*a+b*b+c*c);
			
			currRg += fi*fj*r*r;
		}
	}	
	
	currRg = currRg/2.0;
	
	//Volume correction for radius of gyration.
	double Si;
	for(i=0;i<ntotal;i++){
		fi = pow(p->beads[i].r,3)/RadiusSum;
		Si = 3*p->beads[i].x*p->beads[i].x/5.0;
		currRg += fi*Si; 
	}
	
	currRg = sqrt(currRg);
	/*we need some unit converts here*/
	/*we convert nm to cm*/
	currRg = 1e-7*currRg;
	Rg = add_average(currRg,Rg,n);
	
	//Intrinsic viscosity.Cij is needed when caculate it.
	double curreta = 0;
	static MAT *riT, *rj,*tempvector,*tempvector1,*tempx,*C;
	static int is_malloc2 = 0;
	if(!is_malloc2){
		riT = m_get(1,3);
		tempvector = m_get(1,3);
		tempvector1 = m_get(1,3);
		rj = m_get(3,1);
		tempx = m_get(1,1);
		C = m_get(3,3);
		is_malloc2 = 1;
	}
	
	for(i=0;i<ntotal;i++){
		//may have problem here.rotation center O?
		riT->me[0][0] = p->beads[i].x;
		riT->me[0][1] = p->beads[i].y;
		riT->me[0][2] = p->beads[i].z;
		
		for(j=0;j<ntotal;j++){
			rj->me[0][0] = p->beads[j].x;
			rj->me[1][0] = p->beads[j].y;
			rj->me[2][0] = p->beads[j].z;
			m_mlt(riT,etaD,tempvector);

			for(k=0;k<3;k++){
				for(l=0;l<3;l++){
					C->me[k][l] = BigC->me[i*3+k][j*3+l];
				}
			}

			m_mlt(tempvector,C,tempvector1);
			m_mlt(tempvector1,etaE,tempvector);
			m_mlt(tempvector,rj,tempx);
			
			curreta += tempx->me[0][0];
		}
	}
	
	curreta = Avogadro/(2.0*rm)*curreta;
	
	//Volume correction
	double Volume = 0;
	for(i=0;i<p->nbd;i++)
		Volume += 4/3.0*pi*pow(p->beads[i].r,3);
	curreta += 5*Avogadro*Volume/(2.0*rm);
	
	/*we need some unit converts here*/
	/*convert nm to cm*/
	curreta = curreta*1E-21;
	
	//Intrinsic viscosity.Monte Carlo average.
	eta = add_average(curreta,eta,n);
}
