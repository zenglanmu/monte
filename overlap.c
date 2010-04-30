#include "global.h"
#include "overlap.h"

int ls_overlap(confor *p)
//Check overlap of bead i with all the other beads.
//return 0 if there is no overlap
//Otherwise, return 1
{
	int i;
	int ii;
	double dist,distmin;
	
	for(i=0;i < p->nbd;i++){
		for(ii=0;ii < p->nbd;ii++){
			if(ii != i){
				dist=sqrt(pow((p->beads[i].x-p->beads[ii].x),2)+
						pow((p->beads[i].y-p->beads[ii].y),2)+
						pow((p->beads[i].z-p->beads[ii].z),2));
				distmin=(p->beads[ii].r+p->beads[i].r)*0.99;
				if(dist<distmin) return 1;
			}
		}
	}
	return 0;
}
