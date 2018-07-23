#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "core/meshgen.h"

int nelem=0,nvert=0,ncide=0,nfld=1,npts=0;

mesh meshes[max_mesh];
quad elems[max_elem];
point verts[4*max_elem];
edge cides[4*max_elem];

#include "core/spaces.c"
#include "core/triangles.c"
#include "core/ptops.c"
#include "core/conops.c"
#include "core/output.c"

int main(int argc, char *argv[]){

  char fname[64];
  char bcs[4][2][4];

  point corners[4];
  point axis;
  int i,j,k,n;

  sprintf(fname,"newmesh.rea");

  for(i=1;i<argc;i++){
    if(strncmp(argv[i],"-f",2)==0){
      i++;
      strcpy(fname,argv[i]);
      strcat(fname,".rea");
    }
  }

  point points[24];
  connector cons[4];
  
  points[0].x=0.0;
  points[0].y=0.0;
  points[1].x=1.0;
  points[1].y=0.0;
  points[2].x=1.0;
  points[2].y=1.0;
  points[3].x=0.0;
  points[3].y=1.0;
  n=5;

  npts=4;

  cons[0]=make_lin_side(points[0],points[1],n);
  cons[1]=make_lin_side(points[1],points[2],n);
  cons[2]=make_lin_side(points[2],points[3],n);
  cons[3]=make_lin_side(points[3],points[0],n);

  sprintf(fname,"pts.dat");
  output_pts(points,npts,fname);
  sprintf(fname,"vts.dat");
  if(nvert>0) output_pts(verts,nvert,fname);
  sprintf(fname,"cns.dat");
  output_cons(cons,4,fname);

return 0;
}
