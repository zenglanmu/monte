/*conformation.h*/
#ifndef CONFORMATION_H
#define CONFORMATION_H
extern void GetSpring(confor *p);
extern void GetAngel(confor *p);
extern void McMove(confor *old_p,confor *new_p);
extern void InitialConfor(confor *p);
extern void SavePDBFile(confor *p,char *filename);
#endif
