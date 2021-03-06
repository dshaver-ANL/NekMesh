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

