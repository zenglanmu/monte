#include "global.h"
#include "overlap.h"

int ls_overlap(struct confor p,int n)
//Check overlap of bead i with all the other beads.
//n,total number of beads.
//return 0 if there is no overlap
//Otherwise, return 1
{
	int i;
	int ii;
	double dist,distmin;
	for(i=0;i<n;i++){
		for(ii=0;ii<n;ii++){
			if(ii != i){
				dist=sqrt(pow(((p.beads[i].x)-(p.beads[ii].x)),2)+pow(((p.beads[i].y)-(p.beads[ii].y)),2)+pow(((p.beads[i].z)-(p.beads[ii].z)),2));
				distmin=(p.beads[ii].r+p.beads[i].r)*0.99;
				if(dist<distmin) return 1;
			}
		}
	}
	return 0;
}
