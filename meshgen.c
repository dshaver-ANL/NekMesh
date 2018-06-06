#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "meshgen.h"

int nelem=0,nvert=0,ncide=0;

int main(int argc, char *argv[]){

  quad elems[max_elem];
  point verts[4*max_elem];
  edge cides[4*max_elem];

  char fname[64];

  double x[4],y[4],x0[4],y0[4],A,B,C,R,delta,RR;
  double xc,yc,phi,dphi,theta,psi;

  int i,j,newvert,newelem,newside;

  sprintf(fname,"newmesh.rea");

  for(i=1;i<argc;i++){
    if(strncmp(argv[i],"-f",2)==0){
      i++;
      strcpy(fname,argv[i]);
      strcat(fname,".rea");
    }
  }

  R=0.5;
  delta=0.005;
  x[2]=0.0; y[2]=0.0;
  x[0]=R;   y[0]=0.0;
  x[1]=0.0 ;y[1]=R;

  make_ctri_space(elems,verts,cides,3,8,R,delta,x,y,&newelem,&newvert);

  write_rea(elems,verts,cides,fname);

return 0;
}

int write_rea(quad *elems,point *verts,edge *cides,char *fname){

  char line[256];

  FILE *reain,*reaout;

  int i,j,ivert;

  reaout=fopen(fname,"w");

  reain=fopen("reatop.txt","r");
  while(fgets(line,256,reain)!=NULL) fprintf(reaout,"%s",line);
  fclose(reain);

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
  fprintf(reaout,"  ***** NO FLUID   BOUNDARY CONDITIONS *****\n");

  reain=fopen("reabot.txt","r");
  while(fgets(line,256,reain)!=NULL) fprintf(reaout,"%s",line);
  fclose(reain);

  fclose(reaout);

return 0;
}

int make_quad_space(quad *elem,point *vert,int nr,int ns,double *x,double *y,int *ielem,int *ivert){

  int i,j;
  double *r,*s,dr,ds,r0,s0;
  double tx,ty;

  printf("Generating %d elements in a quadrilateral space\n",ns*nr);

  r=(double *)malloc(sizeof(double)*2*(ns+1));
  s=(double *)malloc(sizeof(double)*2*(ns+1));

  *ivert=0;
  *ielem=0;

//      3---s3---2
//      |        |
//  j   s4      s2
//  ᴧ   |        |
//  |   0---s1---1
//  |
//  *---->i

//side 2 (1-2)
  dr=(x[2]-x[1])/(double)ns;
  ds=(y[2]-y[1])/(double)ns;
  for(j=0;j<=ns;j++){
    *(r+j)=x[1]+(double)j*dr;
    *(s+j)=y[1]+(double)j*ds;
  }

//side 4 (0-3)
  dr=(x[3]-x[0])/(double)ns;
  ds=(y[3]-y[0])/(double)ns;
  for(j=0;j<=ns;j++){
    *(r+j+ns+1)=x[0]+(double)j*dr;
    *(s+j+ns+1)=y[0]+(double)j*ds;
  }

  for(j=0;j<=ns;j++){
    r0=*(r+j+ns+1);
    s0=*(s+j+ns+1);
    dr=(*(r+j)-r0)/(double)nr;
    ds=(*(s+j)-s0)/(double)nr;
    for(i=0;i<=nr;i++){
      (vert+i+j*(nr+1))->x=r0+(double)i*dr;
      (vert+i+j*(nr+1))->y=s0+(double)i*ds;
      (*ivert)++;
    }
  }

  for(j=0;j<ns;j++){
    for(i=0;i<nr;i++){
      (elem+*ielem)->vid[0]=(i+0)+(j+0)*(nr+1)+nvert;
      (elem+*ielem)->vid[1]=(i+1)+(j+0)*(nr+1)+nvert;
      (elem+*ielem)->vid[2]=(i+1)+(j+1)*(nr+1)+nvert;
      (elem+*ielem)->vid[3]=(i+0)+(j+1)*(nr+1)+nvert;
      (*ielem)++;
    }
  }

  free(r);free(s);
  nvert+=*ivert;
  nelem+=*ielem;
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_gquad_space(quad *elem,point *vert,int nr,int ns,double dn0,double *x,double *y,int *ielem,int *ivert){

  int i,j;
  double *r,*s,dr,ds,r0,s0,dr0,ds0;
  double tx,ty,delta,ratio;

  printf("Generating %d elements in a quadrilateral space with geometric growth\n",ns*nr);

  r=(double *)malloc(sizeof(double)*2*(ns+1));
  s=(double *)malloc(sizeof(double)*2*(ns+1));

  *ivert=0;
  *ielem=0;

//side 2 (1-2)
  get_g_side(x[1],x[2],y[1],y[2],ns,dn0,r,s);
//side 4 (0-3)
  get_g_side(x[0],x[3],y[0],y[3],ns,dn0,r+ns+1,s+ns+1);

  for(j=0;j<=ns;j++){
    r0=*(r+j+ns+1);
    s0=*(s+j+ns+1);
    dr=(*(r+j)-r0)/(double)nr;
    ds=(*(s+j)-s0)/(double)nr;
    for(i=0;i<=nr;i++){
      (vert+i+j*(nr+1))->x=r0+(double)i*dr;
      (vert+i+j*(nr+1))->y=s0+(double)i*ds;
      (*ivert)++;
    }
  }

  for(j=0;j<ns;j++){
    for(i=0;i<nr;i++){
      (elem+*ielem)->vid[0]=(i+0)+(j+0)*(nr+1)+nvert;
      (elem+*ielem)->vid[1]=(i+1)+(j+0)*(nr+1)+nvert;
      (elem+*ielem)->vid[2]=(i+1)+(j+1)*(nr+1)+nvert;
      (elem+*ielem)->vid[3]=(i+0)+(j+1)*(nr+1)+nvert;
      (*ielem)++;
    }
  }

  free(r);free(s);
  nvert+=*ivert;
  nelem+=*ielem;
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_cquad_space(quad *elem,point *vert,edge *cide,int nt,int nr,double R,double dr0,double *x,double *y,int *ielem,int *ivert,int *icide){
  
  int i,j,k,istp;
  double theta,dtheta,delta;
  double *r,*s,r0[2],r1[2],s0[2],s1[2];
  double rr,xc,yc,xa,ya,xb,yb,xd,yd;

  printf("Generating %d elements in a curved quadrilateral space with geometric growth\n",nt*nr);
  
  *ivert=0;
  *ielem=0;
  *icide=0;

  circle_center_2pR(x[0],y[0],x[1],y[1],R,&xc,&yc);
  rr = R-dr0;

  r=(double *)malloc(sizeof(double)*(nr+1));
  s=(double *)malloc(sizeof(double)*(nr+1));

//side 2 (1-2)
  r1[0]=x[1];s1[0]=y[1];
  line_circle_intercept(x[1],y[1],x[2],y[2],xc,yc,rr,r1+1,s1+1);

//side 4 (0-3)
  r0[0]=x[0];s0[0]=y[0];
  line_circle_intercept(x[0],y[0],x[3],y[3],xc,yc,rr,r0+1,s0+1);

  for(i=0;i<=nt;i++){
    xb=x[3]+(x[2]-x[3])*((double)i/(double)(nt));
    yb=y[3]+(y[2]-y[3])*((double)i/(double)(nt));
    for(j=0;j<=nr;j++){
      if(j<2){
        rr=sqrt((r0[j]-xc)*(r0[j]-xc)+(s0[j]-yc)*(s0[j]-yc));
        delta=sqrt((r0[j]-r1[j])*(r0[j]-r1[j])+(s0[j]-s1[j])*(s0[j]-s1[j]));
        theta=R/fabs(R)*2.0*asin(delta/(2.0*rr))/((double)nt)*(double)i;
        xa=(r0[j]-xc)*cos(theta)-(s0[j]-yc)*sin(theta)+xc;
        ya=(r0[j]-xc)*sin(theta)+(s0[j]-yc)*cos(theta)+yc;
        (vert+i+j*(nt+1))->x=xa;
        (vert+i+j*(nt+1))->y=ya;
      }else{
        if(j==2){
          xd=(vert+i)->x;
          yd=(vert+i)->y;
          dr0=sqrt((xa-xd)*(xa-xd)+(ya-yd)*(ya-yd));
          get_g_side(xd,xb,yd,yb,nr,dr0,r,s);
        }
        (vert+i+j*(nt+1))->x=*(r+j);
        (vert+i+j*(nt+1))->y=*(s+j);
      }
    }
  }
  (*ivert)+=(nr+1)*(nt+1);

  for(i=0;i<nr;i++){
    for(j=0;j<nt;j++){
      (elem+*ielem)->vid[0]=(j+0)+(i+0)*(nt+1)+nvert;
      (elem+*ielem)->vid[1]=(j+1)+(i+0)*(nt+1)+nvert;
      (elem+*ielem)->vid[2]=(j+1)+(i+1)*(nt+1)+nvert;
      (elem+*ielem)->vid[3]=(j+0)+(i+1)*(nt+1)+nvert;
      if(i<2){
        if(i==0){
          (cide+*icide)->elid=*ielem+nelem;
          (cide+*icide)->esid=1;
          (cide+*icide)->curve=R;
          (cide+*icide)->ccurve='C';
          (*icide)++;
          (cide+*icide)->elid=*ielem+nelem;
          (cide+*icide)->esid=3;
          (cide+*icide)->curve=-(R-dr0);
          (cide+*icide)->ccurve='C';
          (*icide)++;
        }else{
          (cide+*icide)->elid=*ielem+nelem;
          (cide+*icide)->esid=1;
          (cide+*icide)->curve=R-dr0;
          (cide+*icide)->ccurve='C';
          (*icide)++;
        }
      }
      (*ielem)++;
    }
  }

  free(r),free(s);
  nvert+=*ivert;
  nelem+=*ielem;
  ncide+=*icide;
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return *ielem;
}

int make_tri_space(quad *elem,point *vert,int ns,double *x,double *y,int *ielem,int *ivert){

  int i,j,qelem=0,qvert=0;

  double mx[3],my[3],gx=0.0,gy=0.0;
  double xq[4],yq[4];

  *ielem=0;
  *ivert=0;

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

    make_quad_space(elem+nelem,vert+nvert,ns,ns,xq,yq,&qelem,&qvert);
  }

return *ielem;
}

int make_vtri_space(quad *elem,point *vert,int nb,int nh,double *x,double *y,int *ielem,int *ivert){

  int i,j,qelem=0,qvert=0,ns,nr;

  double mx[3],my[3],gx=0.0,gy=0.0;
  double xq[4],yq[4];

  *ielem=0;
  *ivert=0;

  for(i=0;i<3;i++){
    j=(i+1)%3;
    mx[i]=0.5*(x[i]+x[j]);
    my[i]=0.5*(y[i]+y[j]);
    gx+=x[i];
    gy+=y[i];
  }
  gx/=3.0;
  gy/=3.0;

  for(i=0;i<3;i++){
    j=(i+2)%3;
    xq[0]=x[i];  yq[0]=y[i];
    xq[1]=mx[i]; yq[1]=my[i];
    xq[2]=gx;    yq[2]=gy;
    xq[3]=mx[j]; yq[3]=my[j];

    ns=nb; nr=nb;
    if(i==0) nr=nh;
    if(i==1) ns=nh;

    make_quad_space(elem+nelem,vert+nvert,ns,nr,xq,yq,&qelem,&qvert);
  }

return *ielem;
}

int make_ctri_space(quad *elem,point *vert,edge *cide,int nb,int nh,double R,double dr0,double *x,double *y,int *ielem,int *ivert){

  int i,j,qelem=0,qvert=0,qcide=0,ns,nr;

  double mx[3],my[3],gx=0.0,gy=0.0;
  double xq[4],yq[4],xc,yc,aa,w[3];

  *ielem=0;
  *ivert=0;

  w[0]=2.0*(double)nb;
  w[1]=w[0];
  w[2]=(double)nh;
  for(i=0;i<3;i++){
    j=(i+1)%3;
    mx[i]=(x[i]*w[i]+x[j]*w[j])/(w[i]+w[j]);
    my[i]=(y[i]*w[i]+y[j]*w[j])/(w[i]+w[j]);
  }

  circle_center_2pR(x[0],y[0],x[1],y[1],R,&xc,&yc);

  aa=sqrt((mx[0]-xc)*(mx[0]-xc)+(my[0]-yc)*(my[0]-yc));
  mx[0]=(mx[0]-xc)*fabs(R)/aa+xc;
  my[0]=(my[0]-yc)*fabs(R)/aa+yc;

  gx=(3.0*nb*mx[0]+nh*x[2])/(nh+3.0*nb);
  gy=(3.0*nb*my[0]+nh*y[2])/(nh+3.0*nb);

  xq[0]=x[0];  yq[0]=y[0];
  xq[1]=mx[0]; yq[1]=my[0];
  xq[2]=gx;    yq[2]=gy;
  xq[3]=mx[2]; yq[3]=my[2];
  make_cquad_space(elem+nelem,vert+nvert,cide+ncide,nb,nh,R,dr0,xq,yq,&qelem,&qvert,&qcide);

  xq[0]=mx[0]; yq[0]=my[0];
  xq[1]=x[1];  yq[1]=y[1];
  xq[2]=mx[1]; yq[2]=my[1];
  xq[3]=gx;    yq[3]=gy;
  make_cquad_space(elem+nelem,vert+nvert,cide+ncide,nb,nh,R,dr0,xq,yq,&qelem,&qvert,&qcide);

  xq[0]=x[2];  yq[0]=y[2];
  xq[1]=mx[2]; yq[1]=my[2];
  xq[2]=gx;    yq[2]=gy;
  xq[3]=mx[1]; yq[3]=my[1];
  make_quad_space(elem+nelem,vert+nvert,nb,nb,xq,yq,&qelem,&qvert);

return *ielem;
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

int get_g_side(double x1,double x2,double y1, double y2,int n,double d0,double *r, double *s){

  int j;

  double dr=(x2-x1);
  double ds=(y2-y1);
  double delta=sqrt(dr*dr+ds*ds);
  double ratio=growth_ratio(n,delta,d0);
  double dr0=d0*dr/delta;
  double ds0=d0*ds/delta;
  *r=x1;
  *s=y1;
  for(j=1;j<=n;j++){
    *(r+j)=*(r+j-1)+dr0*pow(ratio,(double)(j-1));
    *(s+j)=*(s+j-1)+ds0*pow(ratio,(double)(j-1));
  }

return n;
}

int line_circle_intercept(double x1,double y1,double x2,double y2,double xc,double yc,double rr,double *xi,double *yi){

  double dx2,dy2;
  double dx1=x2-x1;
  double dy1=y2-y1;
  if(fabs(dx1)>1.0e-16){
    double m=dy1/dx1;
    double c=-m*x1+y1;
    double A=(m*m+1.0);
    double B=2.0*(m*c-m*yc-xc);
    double C=yc*yc-rr*rr+xc*xc-2.0*c*yc+c*c;
    double des=B*B-4.0*A*C;
    if(des<0.0){
      printf("error: undefined intercept\n");
      *xi=0.0;
      *yi=0.0;
      return -1;
    }
    //Always want the intercept closest to x1,y1
    dx1=(-B-sqrt(des))/(2.0*A);
    dx2=(-B+sqrt(des))/(2.0*A);
    dy1=m*dx1+c-y1;
    dy2=m*dx2+c-y1;
    dx1-=x1;
    dx2-=x1;
    dx1=sqrt(dx1*dx1+dy1*dy1);
    dx2=sqrt(dx2*dx2+dy2*dy2);
    if(dx1<dx2){
      *xi=(-B-sqrt(des))/(2.0*A);
    }else{
      *xi=(-B+sqrt(des))/(2.0*A);
    }
    *yi=m* *xi+c;
  }else{
    *xi=x1;
    dy1= sqrt(rr*rr-(x1-xc)*(x1-xc))+yc;
    dy2=-sqrt(rr*rr-(x1-xc)*(x1-xc))+yc;
    if(fabs(dy1-y1)<fabs(dy2-y1)){
     *yi=dy1;
    }else{
     *yi=dy2;
    }
  }

return 0;
}

int circle_center_2pR(double x0,double y0,double x1,double y1,double R,double *xc,double *yc){

  double delta = sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
  double aa = delta/2.0;
  double bb = sqrt(R*R-aa*aa);
  double xa = 0.5*(x1-x0); double ya = 0.5*(y1-y0);
  if(R>0.0){ 
    *xc = x0+xa-bb*ya/aa; 
    *yc = y0+ya+bb*xa/aa;
  }else if(R<0.0){
    *xc = x0+xa+bb*ya/aa;
    *yc = y0+ya-bb*xa/aa;
  }

return 0;
}





