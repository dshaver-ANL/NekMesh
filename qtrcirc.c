#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "core/meshgen.h"

int nquad=0,npts=0,nelem=0,nvert=0,ncide=0,nfld=2;

block blocks[max_block];
quad elems[max_elem];
point verts[4*max_elem];
edge cides[4*max_elem];

point *points;
point origin={0.0, 0.0};

#include "core/spaces.c"
#include "core/triangles.c"
#include "core/ptops.c"
#include "core/conops.c"
#include "core/output.c"
#include "core/read_inp.c"

int main(int argc, char *argv[]){

  char inname[64],reaname[64];
  char bcs[4][2][4];

  point corners[4];

  int i,j,k,n;

  double delta=0.040, R=2.5;
  int a=2,b=3,c=3,d=3,e=3,f=2,g=3;

  for(i=1;i<argc;i++){
    if(strncmp(argv[i],"-fi",3)==0){
      i++;
      strcpy(inname,argv[i]);
    }else if(strncmp(argv[i],"-rea",4)==0){
      i++;
      strcpy(reaname,argv[i]);
    }
  }

//read the points from the input file
  read_inp(inname);

  corners[0]=points[0];
  corners[1]=points[1];
  corners[2]=points[9];
  corners[3]=points[6];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[3][0],"SYM");
  sprintf(bcs[0][1],"I  ");
  sprintf(bcs[3][1],"I  ");
  make_gquad_space(2*a,b,delta,corners,bcs);

  corners[0]=points[9];
  corners[1]=points[4];
  corners[2]=points[3];
  corners[3]=points[6];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[2][0],"SYM");
  sprintf(bcs[2][1],"I  ");
  make_quad_space(g,2*a,corners,bcs);

  corners[0]=points[2];
  corners[1]=points[7];
  corners[2]=points[9];
  corners[3]=points[1];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[3][0],"W  ");
  sprintf(bcs[0][1],"I  ");
  sprintf(bcs[3][1],"I  ");
  make_cquad_space(b,c,R,delta,delta,corners,bcs);

  corners[0]=points[7];
  corners[1]=points[8];
  corners[2]=points[4];
  corners[3]=points[9];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[0][1],"I  ");
  make_cquad_space(g,c,R,delta,0.0,corners,bcs);

  corners[0]=points[8];
  corners[1]=points[5];
  corners[2]=points[3];
  corners[3]=points[4];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[1][0],"SYM");
  sprintf(bcs[0][1],"I  ");
  sprintf(bcs[1][1],"I  ");
  make_cquad_space(2*a,c,R,delta,0.0,corners,bcs);
  
  write_rea(reaname);
    
  sprintf(reaname,"pts.dat");
  output_pts(points,npts,reaname);
  sprintf(reaname,"vts.dat");
  if(nvert>0) output_pts(verts,nvert,reaname);

return 0;
}
