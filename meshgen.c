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
  point axis;
  int i,j,k,n;
//  double x[12] = {1.2999996, 1.5999997, 1.2559620, 1.5457994, 1.5884249, 1.9850029, 1.9847418, 1.5881638, 1.3, 1.6, 1.2557036, 1.5454813};
//  double y[12] = {0.001, 0.001, 0.33549873, 0.41292152, 0.20811144, 0.26032196, 0.26230485, 0.21009436, 0.0, 0.0, 0.33646476, 0.41411047};
  int box[3][4] = { {9, 10, 2, 1}, {3, 4, 12, 11}, {5, 6, 7, 8}};
  int nx[3] = {15,15,25};
  int ny = 3;

  sprintf(fname,"newmesh.inp");
  read_inp(fname);

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

  for(i=0;i<3;i++){
    corners[0]=points[box[i][0]-1];
    corners[1]=points[box[i][1]-1];
    corners[2]=points[box[i][2]-1];
    corners[3]=points[box[i][3]-1];

    make_quad_space(nx[i],ny,corners,bcs);
  }

  write_rea(fname);
    
//sprintf(fname,"pts.dat");
//output_pts(points,npts,fname);
  sprintf(fname,"vts.dat");
  if(nvert>0) output_pts(verts,nvert,fname);

return 0;
}
