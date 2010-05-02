/*io.h*/
#ifndef CONFORMATION_H
#define CONFORMATION_H
extern void SavePDBFile(const confor *p,char *filename);
extern void confor_foutput(FILE *fp, const confor *p);
#define confor_output(conf)	confor_foutput(stdout,conf)
#endif
