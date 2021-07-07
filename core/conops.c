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

int initialize_con(connector *r,int ne){

  r->n=ne+1;
  r->p=malloc(sizeof(point)*(r->n));

return 0;
}

int dealloc_con(connector *r){free(r->p); return 0;}

int set_lin_side(point p1,point p2,connector *r){

  int i;
  double dx=(p2.x-p1.x)/(double)(r->n-1);
  double dy=(p2.y-p1.y)/(double)(r->n-1);

  r->p[0]=p1;
  for(i=1;i<r->n;i++){
    r->p[i].x=r->p[i-1].x+dx;
    r->p[i].y=r->p[i-1].y+dy;
  }

return 0;
}

int set_geo_side(point p1,point p2,double d0,connector *r){

  int i;
  double dx=(p2.x-p1.x),dy=(p2.y-p1.y);
  double delta=sqrt(dx*dx+dy*dy),gr=growth_ratio((r->n)-1,delta,d0),ratio=1.0;
  double dx0=d0*dx/delta,dy0=d0*dy/delta;

  r->p[0]=p1;
  for(i=1;i<r->n;i++){
    r->p[i].x=r->p[i-1].x+dx0*ratio;
    r->p[i].y=r->p[i-1].y+dy0*ratio;
    ratio*=gr;
  }

return 0;
}

int set_arc_side(point p1,point p2,double R,connector *r){

  int i;
  point pc=circle_center_2pR(p1,p2,R);
  double dtheta=angle(pc,p1,p2)/(double)(r->n-1)*R/fabs(R);

  for(i=0;i<r->n;i++) r->p[i]=rotate_point(p1,(double)i*dtheta,pc);

return 0;
}

int set_garc_side(point p1,point p2,point pc,double s0,connector *r){

  int i;
  double *th,sgn;
  double R=distance(p1,pc);
  double theta=angle(pc,p1,p2);

  sgn=crossprod2d(pc,p1,p2);
  sgn=sgn/fabs(sgn);

//printf("%f, %f, %f\n",sgn,theta*180./M_PI,R);

  th=malloc(sizeof(double)* r->n);
  get_g1D(0.0,theta,r->n-1,s0,th);
//for(i=0;i<r->n;i++)printf("th %f\n",th[i]*180./M_PI);
  for(i=0;i<r->n;i++)r->p[i]=rotate_point(p1,th[i]*sgn,pc);
  free(th);

return 0;
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

