/*io.h*/
#ifndef IO_H
#define IO_H
extern void SavePDBFile(const confor *p,char *filename);
extern void ReadPDBFILE(confor *p,char *filename);
extern void confor_foutput(FILE *fp, const confor *p);
#define confor_output(conf)	confor_foutput(stdout,conf)
extern void results_output(FILE *fp);
#endif
