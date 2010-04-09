/*conformation.h*/
#ifndef CONFORMATION_H
#define CONFORMATION_H
extern struct confor McMove(struct confor p);
extern struct confor InitialConfor();
extern void SavePDBFile(struct confor p,char *);
#endif
