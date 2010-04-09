/*conformation.h*/
#ifndef CONFORMATION_H
#define CONFORMATION_H
extern struct confor mcmove(struct confor p);
extern struct confor initconfor(double r,double R,double theta,double h,double N);
extern void SavePDBFile(struct confor p,int N,char *);
#endif
