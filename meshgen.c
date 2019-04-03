#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "core/meshgen.h"

int nquad=0,npts=0,nelem=0,nvert=0,ncide=0,nfld=1;

block blocks[max_block];
quad elems[max_elem];
point verts[4*max_elem];
edge cides[4*max_elem];

point *points;

#include "core/spaces.c"
#include "core/triangles.c"
#include "core/ptops.c"
#include "core/conops.c"
#include "core/output.c"
#include "core/read_inp.c"

int main(int argc, char *argv[]){

  char fname[64];
  char bcs[4][2][4];

  point corners[4];

  double rr[6]={-1.3,1.6,-1.3,1.6,-1.602,2.002};
  int i,j,k,n;
  int box[6][4] = { {1, 9, 14, 13}, {10, 2, 13, 14}, {11, 3, 15, 16}, {4, 12, 16, 15}, {8, 5, 17, 18},{6, 7, 18, 17}};
  int nx[6] = {15,15,15,15,25,25};
  int ny = 3;

  sprintf(fname,"newmesh.inp");
  read_inp(fname);

  for(i=0;i<npts;i++){
    printf("%d: %15.8e, %15.8e\n",i,points[i].x,points[i].y);
  } 
  printf("\n");

  sprintf(fname,"newmesh.rea");

  for(i=1;i<argc;i++){
    if(strncmp(argv[i],"-f",2)==0){
      i++;
      strcpy(fname,argv[i]);
      strcat(fname,".rea");
    }
  }

  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[1][0],"W  ");
  sprintf(bcs[2][0],"W  ");
  sprintf(bcs[3][0],"W  ");

  for(i=0;i<6;i++){
    corners[0]=points[box[i][0]-1];
    corners[1]=points[box[i][1]-1];
    corners[2]=points[box[i][2]-1];
    corners[3]=points[box[i][3]-1];
    for(j=0;j<4;j++){
      printf("%d: %15.8e, %15.8e\n",j,corners[j].x,corners[j].y);
    }
    printf("\n");

    make_cquad_space(ny,nx[i],rr[i],0.0,0.0,corners,bcs);
  }

  write_rea(fname);
    
//sprintf(fname,"pts.dat");
//output_pts(points,npts,fname);
  sprintf(fname,"vts.dat");
  if(nvert>0) output_pts(verts,nvert,fname);

return 0;
}
