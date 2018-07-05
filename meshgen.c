#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "meshgen.h"

int nelem=0,nvert=0,ncide=0;

block blocks[max_block];
quad elems[max_elem];
point verts[4*max_elem];
edge cides[4*max_elem];

int main(int argc, char *argv[]){

  char fname[64];
  char bcs[4][4];

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

return 0;
}

int write_rea(char *fname){

  char line[256],bcout[4];

  FILE *reain,*reaout;

  int i,j,ivert;
  double dx;
  double xmin=1.0e30,xmax=-1.0e30,ymin=1.0e30,ymax=-1.0e30;

  for(i=0;i<nvert;i++){
    xmax=fmax(xmax,(verts+i)->x);
    xmin=fmin(xmin,(verts+i)->x);
    ymax=fmax(ymax,(verts+i)->y);
    ymin=fmin(ymin,(verts+i)->y);
  }

  dx=fmax((xmax-xmin),(ymax-ymin));

  reaout=fopen(fname,"w");

  printf("\n\twriting %d elements to file \"%s\"\n\n",nelem,fname);

  reain=fopen("reatop.txt","r");
  while(fgets(line,256,reain)!=NULL) fprintf(reaout,"%s",line);
  fclose(reain);

  fprintf(reaout,"%10.6f%10.6f%10.6f%10.6f\n",1.2*dx,1.2*dx,xmin-0.1*dx,ymin-0.1*dx);
  fprintf(reaout,"**MESH DATA**\n");

  fprintf(reaout,"%5d 2 %5d     NEL, NDIM, NELV\n",nelem,nelem);
  for(i=0;i<nelem;i++){
    fprintf(reaout,"            ELEMENT%5d [    1a]    GROUP     0\n",i+1);
    for(j=0;j<4;j++){
      ivert=(elems+i)->vid[j];
      fprintf(reaout,"%15.7e",(verts+ivert)->x);
    }
    fprintf(reaout,"\n");
    for(j=0;j<4;j++){
      ivert=(elems+i)->vid[j];
      fprintf(reaout,"%15.7e",(verts+ivert)->y);
    }
    fprintf(reaout,"\n");
  }

  fprintf(reaout,"  ***** CURVED SIDE DATA *****\n");
  if(ncide==0){
    fprintf(reaout,"     0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n");
  }else{
    fprintf(reaout,"%6d Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n",ncide);
    if(nelem<1000){for(i=0;i<ncide;i++)
      fprintf(reaout,"%3d%3d %9.6f %13.5f %13.5f %13.5f %13.5f     %c\n",
        (cides+i)->esid,(cides+i)->elid+1,(cides+i)->curve,0.0,0.0,0.0,0.0,(cides+i)->ccurve);
    }else{ for(i=0;i<ncide;i++)
      fprintf(reaout,"%2d%6d %9.6f %13.5f %13.5f %13.5f %13.5f     %c\n",
        (cides+i)->esid,(cides+i)->elid+1,(cides+i)->curve,0.0,0.0,0.0,0.0,(cides+i)->ccurve);
    }
  }

  fprintf(reaout,"  ***** BOUNDARY CONDITIONS *****\n");
  fprintf(reaout,"  ***** FLUID   BOUNDARY CONDITIONS *****\n");
  for(i=0;i<nelem;i++){
    for(j=0;j<4;j++){
       strcpy(bcout,(elems+i)->BC[j]);
       fprintf(reaout," %s%3d  %1d   %6f       %6f       %6f       %6f       %6f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
    }
  }

  reain=fopen("reabot.txt","r");
  while(fgets(line,256,reain)!=NULL) fprintf(reaout,"%s",line);
  fclose(reain);

  fclose(reaout);

return 0;
}

int make_quad_space(int nr,int ns,point *p,char bcs[4][4]){

  int i,j,nvert0=nvert;
  double *r,*s,dr,ds,r0,s0;

  printf("Generating %d elements in a quadrilateral space\n",ns*nr);

  r=(double *)malloc(sizeof(double)*2*(ns+1));
  s=(double *)malloc(sizeof(double)*2*(ns+1));

//      3---s3---2
//      |        |
// j,s  s4      s2
//  ᴧ   |        |
//  |   0---s1---1
//  |
//  *---->i,r

//side 2 (1-2)
  dr=(p[2].x-p[1].x)/(double)ns;
  ds=(p[2].y-p[1].y)/(double)ns;
  for(j=0;j<=ns;j++){
    *(r+j)=p[1].x+(double)j*dr;
    *(s+j)=p[1].y+(double)j*ds;
  }

//side 4 (0-3)
  dr=(p[3].x-p[0].x)/(double)ns;
  ds=(p[3].y-p[0].y)/(double)ns;
  for(j=0;j<=ns;j++){
    *(r+j+ns+1)=p[0].x+(double)j*dr;
    *(s+j+ns+1)=p[0].y+(double)j*ds;
  }

  for(j=0;j<=ns;j++){
    r0=*(r+j+ns+1);
    s0=*(s+j+ns+1);
    dr=(*(r+j)-r0)/(double)nr;
    ds=(*(s+j)-s0)/(double)nr;
    for(i=0;i<=nr;i++){
      (verts+nvert0+i+j*(nr+1))->x=r0+(double)i*dr;
      (verts+nvert0+i+j*(nr+1))->y=s0+(double)i*ds;
      nvert++;
    }
  }

  for(j=0;j<ns;j++){
    for(i=0;i<nr;i++){
      (elems+nelem)->vid[0]=(i+0)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[1]=(i+1)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[2]=(i+1)+(j+1)*(nr+1)+nvert0;
      (elems+nelem)->vid[3]=(i+0)+(j+1)*(nr+1)+nvert0;
      if(i==0){
        strcpy((elems+nelem)->BC[3],bcs[3]);
        sprintf((elems+nelem)->BC[1],"E  ");
      }else if(i==nr-1){
        strcpy((elems+nelem)->BC[1],bcs[1]);
        sprintf((elems+nelem)->BC[3],"E  ");
      }else{
        sprintf((elems+nelem)->BC[1],"E  ");
        sprintf((elems+nelem)->BC[3],"E  ");
      }
      if(j==0){
        strcpy((elems+nelem)->BC[0],bcs[0]);
        sprintf((elems+nelem)->BC[2],"E  ");
      }else if(j==ns-1){
        strcpy((elems+nelem)->BC[2],bcs[2]);
        sprintf((elems+nelem)->BC[0],"E  ");
      }else{
        sprintf((elems+nelem)->BC[0],"E  ");
        sprintf((elems+nelem)->BC[2],"E  ");
      }
      nelem++;
    }
  }

  free(r);free(s);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_gquad_space(int nr,int ns,double dn0,point *p,char bcs[4][4]){

  int i,j,nvert0=nvert;
  double dr,ds;
  connector r2,r4;

  printf("Generating %d elements in a quadrilateral space with geometric growth\n",ns*nr);

  r2.n=ns+1;
  r2.p=(point *)malloc(sizeof(point)*r2.n);
  r4.n=ns+1;
  r4.p=(point *)malloc(sizeof(point)*r4.n);

//side 2 (1-2)
  get_g_side(p[1],p[2],dn0,&r2);
//side 4 (0-3)
  get_g_side(p[0],p[3],dn0,&r4);

  for(j=0;j<=ns;j++){
    dr=(r2.p[j].x-r4.p[j].x)/nr;
    ds=(r2.p[j].y-r4.p[j].y)/nr;
    for(i=0;i<=nr;i++){
      verts[nvert0+i+j*(nr+1)].x=r4.p[j].x+(double)i*dr;
      verts[nvert0+i+j*(nr+1)].y=r4.p[j].y+(double)i*ds;
      nvert++;
    }
  }

  for(j=0;j<ns;j++){
    for(i=0;i<nr;i++){
      (elems+nelem)->vid[0]=(i+0)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[1]=(i+1)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[2]=(i+1)+(j+1)*(nr+1)+nvert0;
      (elems+nelem)->vid[3]=(i+0)+(j+1)*(nr+1)+nvert0;
      if(i==0){
        strcpy((elems+nelem)->BC[3],bcs[3]);
        sprintf((elems+nelem)->BC[1],"E  ");
      }else if(i==nr-1){
        strcpy((elems+nelem)->BC[1],bcs[1]);
        sprintf((elems+nelem)->BC[3],"E  ");
      }else{
        sprintf((elems+nelem)->BC[1],"E  ");
        sprintf((elems+nelem)->BC[3],"E  ");
      }
      if(j==0){
        strcpy((elems+nelem)->BC[0],bcs[0]);
        sprintf((elems+nelem)->BC[2],"E  ");
      }else if(j==ns-1){
        strcpy((elems+nelem)->BC[2],bcs[2]);
        sprintf((elems+nelem)->BC[0],"E  ");
      }else{
        sprintf((elems+nelem)->BC[0],"E  ");
        sprintf((elems+nelem)->BC[2],"E  ");
      }
      nelem++;
    }
  }

  free(r2.p);free(r4.p);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_g2quad_space(int nr,int ns,double dr0,double ds0,point *p,char bcs[4][4]){

  int i,j,nvert0=nvert;
  point r0,r1;
  connector r,r2,r4;

  printf("Generating %d elements with geometric growth in two directions\n",ns*nr);

  r2.n=ns+1;
  r2.p=(point *)malloc(sizeof(point)*r2.n);
  r4.n=ns+1;
  r4.p=(point *)malloc(sizeof(point)*r4.n);

  r.n=nr+1;
  r.p=(point *)malloc(sizeof(point)*r.n);

//side 2 (1-2)
  get_g_side(p[1],p[2],ds0,&r2);
//side 4 (0-3)
  get_g_side(p[0],p[3],ds0,&r4);

  for(j=0;j<=ns;j++){
    r0=r4.p[j];
    r1=r2.p[j];
    get_g_side(r0,r1,dr0,&r);
    for(i=0;i<=nr;i++) verts[nvert0+i+j*(nr+1)]=r.p[i];
    nvert+=r.n;
  }

  for(j=0;j<ns;j++){
    for(i=0;i<nr;i++){
      (elems+nelem)->vid[0]=(i+0)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[1]=(i+1)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[2]=(i+1)+(j+1)*(nr+1)+nvert0;
      (elems+nelem)->vid[3]=(i+0)+(j+1)*(nr+1)+nvert0;
      if(i==0){
        strcpy((elems+nelem)->BC[3],bcs[3]);
        sprintf((elems+nelem)->BC[1],"E  ");
      }else if(i==nr-1){
        strcpy((elems+nelem)->BC[1],bcs[1]);
        sprintf((elems+nelem)->BC[3],"E  ");
      }else{
        sprintf((elems+nelem)->BC[1],"E  ");
        sprintf((elems+nelem)->BC[3],"E  ");
      }
      if(j==0){
        strcpy((elems+nelem)->BC[0],bcs[0]);
        sprintf((elems+nelem)->BC[2],"E  ");
      }else if(j==ns-1){
        strcpy((elems+nelem)->BC[2],bcs[2]);
        sprintf((elems+nelem)->BC[0],"E  ");
      }else{
        sprintf((elems+nelem)->BC[0],"E  ");
        sprintf((elems+nelem)->BC[2],"E  ");
      }
      nelem++;
    }
  }

  free(r2.p);free(r4.p),free(r.p);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_cquad_space(int nt,int nr,double R,double dr0,double dt0,point *p,char bcs[4][4]){
  
  int i,j,nvert0=nvert;
  double theta,thetb,thet0,delta,rr;
  double *ta,*tb;
  point pc,r0,r1;
  connector r,r2,r3,r4;

  printf("Generating %d elements in a curved quadrilateral space with geometric growth with BCs:\n",nt*nr);
  for(i=0;i<4;i++)printf("\t\"%s\"\n",bcs[i]);
  
  pc=circle_center_2pR(p[0],p[1],R);
  rr = R-dr0;

  r.n=nr+1;
  r.p=(point *)malloc(sizeof(point)*r.n);
  r2.n=nr+1;
  r2.p=(point *)malloc(sizeof(point)*r2.n);
  r4.n=nr+1;
  r4.p=(point *)malloc(sizeof(point)*r4.n);
  r3.n=nt+1;
  r3.p=(point *)malloc(sizeof(point)*r3.n);
  ta=(double *)malloc(sizeof(double)*(nt+1));
  tb=(double *)malloc(sizeof(double)*(nt+1));

//      3-----s3-----2
//      |            |
//      |            |
//      s4          s2
// j,r  |            |
//  ᴧ   | /---s1---\ |
//  |   0/          \1
//  |
//  *---->i,t  (r,theta)


//side 2 (1-2)
  r0=p[1];
  r1=line_circle_intercept(p[1],p[2],pc,rr);
  delta=distance(r0,r1);
  get_g_side(p[1],p[2],delta,&r2);

//side 4 (0-3)
  r0=p[0];
  r1=line_circle_intercept(p[0],p[3],pc,rr);
  delta=distance(r0,r1);
  get_g_side(p[0],p[3],delta,&r4);

  delta=distance(r4.p[0],r2.p[0]);
  theta=2.0*asin(delta/(2.0*fabs(R)));
  
  delta=distance(r4.p[1],r2.p[1]);
  thetb=2.0*asin(delta/(2.0*fabs(rr)));

//side 3 (3-2), theta1a, theta1b
  if(dt0>0.0){ //geometric growth from side 4
    get_g_side(p[3],p[2],dt0,&r3);
    thet0=get_theta_0(r4.p[0],p[3],0,dt0,R,pc);
    get_g1D(0.0,-theta,nt,thet0,ta);
    thet0=get_theta_0(r4.p[1],p[3],0,dt0,rr,pc);
    get_g1D(0.0,-thetb,nt,thet0,tb);
  }else if(dt0<0.0){ //geometric growth from side 2
    get_g_side(p[2],p[3],-dt0,&r3);
    thet0=get_theta_0(r2.p[0],p[2],1,-dt0,R,pc);
    get_g1D(0.0,-theta,nt,thet0,ta);
    //re-orient from side 4
    invert_connector(&r3);
    for(i=0;i<=nt;i++) *(ta+i)=-(theta+*(ta+i));
    invert(ta,nt+1);
  }else{ //linear distribution
    for(i=0;i<=nt;i++){
      r3.p[i].x=p[3].x+(p[2].x-p[3].x)*((double)i/(double)(nt));
      r3.p[i].y=p[3].y+(p[2].y-p[3].y)*((double)i/(double)(nt));
      *(ta+i)=-((double)i)/((double)nt)*theta;
    }
  }

  for(i=0;i<=nt;i++){
    //j=0
    theta=-R/fabs(R)* *(ta+i);
    verts[nvert0+i]=rotate_point(r4.p[0],theta,pc);
    //j=1
    verts[nvert0+i+nt+1]=line_circle_intercept(verts[nvert0+i],r3.p[i],pc,rr);
   
    for(j=2;j<=nr;j++){
      if(j==2) get_g_side(verts[nvert0+i],r3.p[i],dr0,&r);
      verts[nvert0+i+j*(nt+1)]=r.p[j];
    }
  }
  nvert+=(nr+1)*(nt+1);

  for(i=0;i<nr;i++){
    for(j=0;j<nt;j++){
      (elems+nelem)->vid[0]=(j+0)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[1]=(j+1)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[2]=(j+1)+(i+1)*(nt+1)+nvert0;
      (elems+nelem)->vid[3]=(j+0)+(i+1)*(nt+1)+nvert0;
      if(i==0){
        strcpy((elems+nelem)->BC[0],bcs[0]);
        sprintf((elems+nelem)->BC[2],"E  ");
      }else if(i==nr-1){
        strcpy((elems+nelem)->BC[2],bcs[2]);
        sprintf((elems+nelem)->BC[0],"E  ");
      }else{
        sprintf((elems+nelem)->BC[0],"E  ");
        sprintf((elems+nelem)->BC[2],"E  ");
      }
      if(j==0){ //j is build from side 4 [3] to 2 [1]
        strcpy((elems+nelem)->BC[3],bcs[3]);
        sprintf((elems+nelem)->BC[1],"E  ");
      }else if(j==nt-1){
        strcpy((elems+nelem)->BC[1],bcs[1]);
        sprintf((elems+nelem)->BC[3],"E  ");
      }else{
        sprintf((elems+nelem)->BC[1],"E  ");
        sprintf((elems+nelem)->BC[3],"E  ");
      }
      if(i<3){
        if(i==0){
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=1;
          (cides+ncide)->curve=R;
          (cides+ncide)->ccurve='C';
          ncide++;
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=3;
          (cides+ncide)->curve=-(R-dr0);
          (cides+ncide)->ccurve='C';
          ncide++;
        }else if(i==1){
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=1;
          (cides+ncide)->curve=R-dr0;
          (cides+ncide)->ccurve='C';
          ncide++;
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=3;
          (cides+ncide)->curve=-2*(R-dr0);
          (cides+ncide)->ccurve='C';
          ncide++;
        }else{
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=1;
          (cides+ncide)->curve=2*(R-dr0);
          (cides+ncide)->ccurve='C';
          ncide++;
        }
      }
      nelem++;
    }
  }

  free(r.p);
  free(r2.p);
  free(r3.p);
  free(r4.p);
  free(ta);free(tb);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_tri_space(int ns,double *x,double *y){

  int i,j;

  double mx[3],my[3],gx=0.0,gy=0.0;
  double xq[4],yq[4];

  for(i=0;i<3;i++){
    gx+=x[i];
    gy+=y[i];
  }
  gx/=3.0;
  gy/=3.0;

  mx[0]=0.5*(x[0]+x[1]);
  my[0]=0.5*(y[0]+y[1]);
  mx[1]=0.5*(x[1]+x[2]);
  my[1]=0.5*(y[1]+y[2]);
  mx[2]=0.5*(x[2]+x[0]);
  my[2]=0.5*(y[2]+y[0]);

  for(i=0;i<3;i++){
    j=(i+2)%3;
    xq[0]=x[i];  yq[0]=y[i];
    xq[1]=mx[i]; yq[1]=my[i];
    xq[2]=gx;    yq[2]=gy;
    xq[3]=mx[j]; yq[3]=my[j];

//  make_quad_space(ns,ns,xq,yq);
  }

return 0;
}

int make_vtri_space(int nb,int nh,double *x,double *y){

  int i,j,ns,nr;

  double mx[3],my[3],gx=0.0,gy=0.0;
  double xq[4],yq[4];

/*
  for(i=0;i<3;i++){
    j=(i+1)%3;
    mx[i]=0.5*(x[i]+x[j]);
    my[i]=0.5*(y[i]+y[j]);
    gx+=x[i];
    gy+=y[i];
  }
  gx/=3.0;
  gy/=3.0;
*/
  mx[0]=0.5*(x[0]+x[1]);
  my[0]=0.5*(y[0]+y[1]);

  mx[1]=(nb*x[1]+nh*x[2])/(nh+nb);
  my[1]=(nb*y[1]+nh*y[2])/(nh+nb);

  mx[2]=(nh*x[2]+nb*x[0])/(nh+nb);
  my[2]=(nh*y[2]+nb*y[0])/(nh+nb);

  gx=(nh*x[2]+nb*nb*mx[0])/(nh+nb*nb);
  gy=(nh*y[2]+nb*nb*my[0])/(nh+nb*nb);

  for(i=0;i<3;i++){
    j=(i+2)%3;
    xq[0]=x[i];  yq[0]=y[i];
    xq[1]=mx[i]; yq[1]=my[i];
    xq[2]=gx;    yq[2]=gy;
    xq[3]=mx[j]; yq[3]=my[j];

    ns=nb; nr=nb;
    if(i==0) nr=nh;
    if(i==1) ns=nh;

//  make_quad_space(ns,nr,xq,yq);
  }

return 0;
}

int make_ctri_space(int nb,int nh,double R,double dr0,double *x,double *y){

  int i,j,ns,nr;

  double mx[3],my[3],gx=0.0,gy=0.0;
  double xq[4],yq[4],xc,yc,aa,w[3];
  point p1,p2,pc;

  w[0]=2.0*(double)nb;
  w[1]=w[0];
  w[2]=(double)nh;
  for(i=0;i<3;i++){
    j=(i+1)%3;
    mx[i]=(x[i]*w[i]+x[j]*w[j])/(w[i]+w[j]);
    my[i]=(y[i]*w[i]+y[j]*w[j])/(w[i]+w[j]);
  }

  p1.x=x[0];
  p1.y=y[0];
  p2.x=x[1];
  p2.y=y[1];
  pc=circle_center_2pR(p1,p2,R);
  xc=pc.x;
  yc=pc.y;

  aa=sqrt((mx[0]-xc)*(mx[0]-xc)+(my[0]-yc)*(my[0]-yc));
  mx[0]=(mx[0]-xc)*fabs(R)/aa+xc;
  my[0]=(my[0]-yc)*fabs(R)/aa+yc;

  gx=(3.0*nb*mx[0]+nh*x[2])/(nh+3.0*nb);
  gy=(3.0*nb*my[0]+nh*y[2])/(nh+3.0*nb);

  xq[0]=x[0];  yq[0]=y[0];
  xq[1]=mx[0]; yq[1]=my[0];
  xq[2]=gx;    yq[2]=gy;
  xq[3]=mx[2]; yq[3]=my[2];
//make_cquad_space(nb,nh,R,dr0,0.0,xq,yq);

  xq[0]=mx[0]; yq[0]=my[0];
  xq[1]=x[1];  yq[1]=y[1];
  xq[2]=mx[1]; yq[2]=my[1];
  xq[3]=gx;    yq[3]=gy;
//make_cquad_space(nb,nh,R,dr0,0.0,xq,yq);

  xq[0]=x[2];  yq[0]=y[2];
  xq[1]=mx[2]; yq[1]=my[2];
  xq[2]=gx;    yq[2]=gy;
  xq[3]=mx[1]; yq[3]=my[1];
//make_quad_space(nb,nb,xq,yq);

return 0;
}

double growth_ratio(int nx,double delta,double del0){

  double eps,f;
  int istp;
  double N=(double)nx;
  double del1=delta-del0;

  f=del1/(del0*(nx-1.0));
  if(f==1.0) return f;
  for(istp=0,eps=1.0;fabs(eps)>1.0e-12&&istp<1000;istp++){
    double fN=pow(f,N);
    double fNm1=pow(f,N-1.0);
    double onemf=1.0-f;
    double dedf=del0*(-N*fNm1/onemf+(1.-fN)/(onemf*onemf));
    eps=del0*(1.0-fN)/onemf-delta;
//  printf("istp = %d, f = %f\n",istp,f);
    f-=eps/dedf;
  }

return f;
}

int get_g_side(point p1,point p2,double d0,connector *r){

  int j;
  double dr=(p2.x-p1.x);
  double ds=(p2.y-p1.y);
  double delta=sqrt(dr*dr+ds*ds);
  double ratio=growth_ratio((r->n)-1,delta,d0);
  double dr0=d0*dr/delta;
  double ds0=d0*ds/delta;
  (r->p[0])=p1;
  for(j=1;j<(r->n);j++){
    (r->p[j]).x=(r->p[j-1]).x+dr0*pow(ratio,(double)(j-1));
    (r->p[j]).y=(r->p[j-1]).y+ds0*pow(ratio,(double)(j-1));
  }

return j;
}

int get_g1D(double x1,double x2,int n,double d0,double *r){

  int j;
  double dr=(x2-x1);
  double ratio=growth_ratio(n,fabs(dr),d0);
  *r=x1;
   for(j=1;j<=n;j++) *(r+j)=*(r+j-1)+d0*dr/fabs(dr)*pow(ratio,(double)(j-1));

return n;
}

point line_circle_intercept(point p1,point p2,point pc,double rr){

  point pi;
  double dx2,dy2;
  double dx1=p2.x-p1.x;
  double dy1=p2.y-p1.y;
  if(fabs(dx1)>1.0e-12){
    double m=dy1/dx1;
    double c=-m*p1.x+p1.y;
    double A=(m*m+1.0);
    double B=2.0*(m*c-m*pc.y-pc.x);
    double C=pc.y*pc.y-rr*rr+pc.x*pc.x-2.0*c*pc.y+c*c;
    double des=B*B-4.0*A*C;
    if(des<0.0){
      printf("error: undefined intercept\n");
      pi.x=0.0;
      pi.y=0.0;
      return pi;
    }
    //Always want the intercept closest to x1,y1
    dx1=(-B-sqrt(des))/(2.0*A);
    dx2=(-B+sqrt(des))/(2.0*A);
    dy1=m*dx1+c-p1.y;
    dy2=m*dx2+c-p1.y;
    dx1-=p1.x;
    dx2-=p1.x;
    dx1=sqrt(dx1*dx1+dy1*dy1);
    dx2=sqrt(dx2*dx2+dy2*dy2);
    if(dx1<dx2){
      pi.x=(-B-sqrt(des))/(2.0*A);
    }else{
      pi.x=(-B+sqrt(des))/(2.0*A);
    }
    pi.y=m*pi.x+c;
  }else{
    pi.x=p1.x;
    dy1= sqrt(rr*rr-(p1.x-pc.x)*(p1.x-pc.x))+pc.y;
    dy2=-sqrt(rr*rr-(p1.x-pc.x)*(p1.x-pc.x))+pc.y;
    if(fabs(dy1-p1.y)<fabs(dy2-p1.y)){
     pi.y=dy1;
    }else{
     pi.y=dy2;
    }
  }

  return pi;
}

point line_line_intercept(point p1,point p2,point p3, point p4){
  point p0={0.0,0.0};
  double denom=(p1.x-p2.x)*(p3.y-p4.y)-(p1.y-p2.y)*(p3.x-p4.x);
  if(fabs(denom)<1.0e-12){
    printf("Error in line_line_intercept: lines are parallel!");
    return p0;
  }
  p0.x=(p1.x*p2.y-p1.y*p2.x)*(p3.x-p4.x)-(p1.x-p2.x)*(p3.x*p4.y-p3.y*p4.x);
  p0.x/=denom;
  p0.y=(p1.x*p2.y-p1.y*p2.x)*(p3.y-p4.y)-(p1.y-p2.y)*(p3.x*p4.y-p3.y*p4.x);
  p0.y/=denom;

  return p0;
}

point circle_center_2pR(point p0,point p1,double R){

  double delta = sqrt((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y));
  double aa = delta/2.0;
  double bb = sqrt(R*R-aa*aa);
  double xa = 0.5*(p1.x-p0.x); double ya = 0.5*(p1.y-p0.y);
  point pc;
  if(R>0.0){ 
    pc.x = p0.x+xa-bb*ya/aa; 
    pc.y = p0.y+ya+bb*xa/aa;
  }else if(R<0.0){
    pc.x = p0.x+xa+bb*ya/aa;
    pc.y = p0.y+ya-bb*xa/aa;
  }

  return pc;
}

double get_theta_0(point p0,point p1,int icor,double dt0,double R,point pc){
  double cp,del,fact;
  point p2,p3,p4;

  double m=(p1.y-p0.y)/(p1.x-p0.x);
  double B1=p0.y-m*p0.x;
  double B2=B1+dt0*sqrt(1+m*m);
  if(fabs(p1.x-p0.x)<1.0e-10){
    p3.x=p0.x-dt0;
    p3.y=p0.y;
    p4.x=p1.x-dt0;
    p4.y=p0.y;
  }else if(fabs(m)>1){ 
    p3.y=p0.y;
    p3.x=(p3.y-B2)/m;
    p4.y=p1.y;
    p4.x=(p4.y-B2)/m;
  }else{
    p3.x=p0.x;
    p3.y=m*p3.x+B2;
    p4.x=p1.x;
    p4.y=m*p4.x+B2;
  }
  p2=line_circle_intercept(p3,p4,pc,R);

  cp=(p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
  if(((cp>0.0)&&(icor==1))||((cp<0.0)&&(icor==0))){
    del=sqrt((p2.x-p0.x)*(p2.x-p0.x)+(p2.y-p0.y)*(p2.y-p0.y));
    return 2.0*asin(del/(2.0*fabs(R)));
  }
  B2=B1-dt0*sqrt(1+m*m);
  if(fabs(p1.x-p0.x)<1.0e-10){
    p3.x=p0.x+dt0;
    p3.y=p0.y;
    p4.x=p1.x+dt0;
    p4.y=p0.y;
  }else if(fabs(m)>1){ 
    p3.y=p0.y;
    p3.x=(p3.y-B2)/m;
    p4.y=p1.y;
    p4.x=(p4.y-B2)/m;
  }else{
    p3.x=p0.x;
    p3.y=m*p3.x+B2;
    p4.x=p1.x;
    p4.y=m*p4.x+B2;
  }
  p2=line_circle_intercept(p3,p4,pc,R);
  del=sqrt((p2.x-p0.x)*(p2.x-p0.x)+(p2.y-p0.y)*(p2.y-p0.y));
  return 2.0*asin(del/(2.0*fabs(R)));
}

int invert(double *x,int n){

  int i,j;
  double *y;

  y=malloc(sizeof(double)*n);

  for(i=n-1,j=0;i>=0;i--,j++) y[j]=x[i];
  for(i=0;i<n;i++) x[i]=y[i];

  free(y);

return i;
}

int invert_connector(connector *r){

  int i,j;
  point *q;

  q=(point *)malloc(sizeof(point)*(r->n));
  for(i=(r->n-1),j=0;i>=0;i--,j++) q[j]=r->p[i];
  for(i=0;i<(r->n);i++) r->p[i]=q[i];

  free(q);

return i;
}

point rotate_point(point pin,double theta,point zero){
  point pout;

  pout.x = (pin.x-zero.x)*cos(theta)-(pin.y-zero.y)*sin(theta)+zero.x;
  pout.y = (pin.x-zero.x)*sin(theta)+(pin.y-zero.y)*cos(theta)+zero.y;

 return pout;
}

double distance(point p1,point p2){
  double dist=sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
  return dist;
}

point reflect_point(point pin,point p0,point p1){
  point pout;
  double a=p1.y-p0.y;
  double b=p0.x-p1.x;
  double c=(p1.x-p0.x)*p0.y-(p1.y-p0.y)*p0.x;
  double RHS=-2.0*(a*pin.x+b*pin.y+c)/(a*a+b*b);
  pout.x=RHS*a+pin.x;
  pout.y=RHS*b+pin.y;
  return pout;
}

point translate_point(point pin,point del){
  point pout;
  pout.x=pin.x+del.x;
  pout.y=pin.y+del.y;
  pout.z=pin.z+del.z;
  return pout;
}

int reset(void){
  printf("\n\treseting mesh generator\n\n");
  nelem = 0;
  nvert = 0;
  ncide = 0;
}
