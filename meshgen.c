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
point origin={0.0, 0.0};

#include "core/spaces.c"
#include "core/triangles.c"
#include "core/ptops.c"
#include "core/conops.c"
#include "core/output.c"
#include "core/read_inp.c"

int main(int argc, char *argv[]){

#define Nbox 7

  char inname[64],reaname[64];
  char bcs[4][2][4];

  point corners[4];

  double rr[Nbox]={-1.3,1.6,-1.602,2.002};
  int i,j,k,n;
  int box[Nbox][4] = { {1, 2, 16, 15}, {3, 4, 15, 16}, {10, 9, 17, 18}, {12, 11, 18, 17} };
  int nx[Nbox] = {15,15,25,25};
  int ny = 3;
  char bndry[Nbox][4][2][4];


  sprintf(bndry[0][0][0],"W  ");
  sprintf(bndry[0][1][0],"SYM");
  sprintf(bndry[0][2][0],"E  ");
  sprintf(bndry[0][3][0],"W  ");

  sprintf(bndry[1][0][0],"E  ");
  sprintf(bndry[1][1][0],"W  ");
  sprintf(bndry[1][2][0],"E  ");
  sprintf(bndry[1][3][0],"SYM");

  sprintf(bndry[2][0][0],"E  ");
  sprintf(bndry[2][1][0],"W  ");
  sprintf(bndry[2][2][0],"E  ");
  sprintf(bndry[2][3][0],"SYM");

  sprintf(bndry[3][0][0],"E  ");
  sprintf(bndry[3][1][0],"SYM");
  sprintf(bndry[3][2][0],"E  ");
  sprintf(bndry[3][3][0],"W  ");

  sprintf(inname,"newmesh.inp");
  sprintf(reaname,"newmesh.rea");

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
//read_inp(inname);

  corners[1].x=1.0;corners[1].y=0.0;
  corners[2].x=1.3;corners[1].y=0.0;
  corners[3]=rotate_point(corners[2],15.0*M_PI/180.0,origin);
  corners[0]=rotate_point(corners[1],15.0*M_PI/180.0,origin);

  sprintf(bndry[0][0][0],"W  ");
  sprintf(bndry[0][1][0],"W  ");
  sprintf(bndry[0][2][0],"W  ");
  sprintf(bndry[0][3][0],"W  ");

  make_arc_space(5,4,-1.0,corners,bndry[0]);

/*for(i=0;i<Nbox;i++){
    corners[0]=points[box[i][0]-1];
    corners[1]=points[box[i][1]-1];
    corners[2]=points[box[i][2]-1];
    corners[3]=points[box[i][3]-1];

    make_cquad_space(ny,nx[i],rr[i],0.0,0.0,corners,bndry[i]);
  }
*/
  write_rea(reaname);
    
  sprintf(reaname,"pts.dat");
  output_pts(points,npts,reaname);
  sprintf(reaname,"vts.dat");
  if(nvert>0) output_pts(verts,nvert,reaname);

return 0;
}
