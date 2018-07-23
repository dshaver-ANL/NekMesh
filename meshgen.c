#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "core/meshgen.h"

int nelem=0,nvert=0,ncide=0,nfld=1,npts=0;

block blocks[max_block];
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

#include "VTR/VTR.c"

/*  
  corners[0].x=0.0;corners[0].y=0.0;
  corners[1].x=0.8;corners[1].y=-0.1;
  corners[2].x=1.0;corners[2].y=1.0;
  corners[3].x=-0.1;corners[3].y=0.86;
  sprintf(bcs[0],"W  ");
  sprintf(bcs[1],"O  ");
  sprintf(bcs[2],"SYM");
  sprintf(bcs[3],"E  ");

  make_cquad_space(5,6,-1.1,0.03,-0.04,corners,bcs);

  corners[1]=corners[0];
  corners[2]=corners[3];
  corners[0].x=-0.8;corners[0].y=0.1;
  corners[3].x=-0.6;corners[3].y=0.7;
  sprintf(bcs[1],"E  ");
  sprintf(bcs[3],"v  ");

  make_cquad_space(5,6,1.1,0.03,0.00,corners,bcs);
  write_rea(fname);
 */ 
  sprintf(fname,"pts.dat");
  output_pts(points,npts,fname);
  sprintf(fname,"vts.dat");
  if(nvert>0) output_pts(verts,nvert,fname);

return 0;
}
