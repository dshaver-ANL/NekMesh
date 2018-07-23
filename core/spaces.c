int make_quad_space(int nr,int ns,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
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
      for(k=0;k<nfld;k++){
        if(i==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
        if(j==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(j==ns-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
      }
      nelem++;
    }
  }

  free(r);free(s);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_gquad_space(int nr,int ns,double dn0,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
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
      for(k=0;k<nfld;k++){
        if(i==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
        if(j==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(j==ns-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
      }
      nelem++;
    }
  }

  free(r2.p);free(r4.p);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_g2quad_space(int nr,int ns,double dr0,double ds0,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
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
      for(k=0;k<nfld;k++){
        if(i==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
        if(j==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(j==ns-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
      }
      nelem++;
    }
  }

  free(r2.p);free(r4.p),free(r.p);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_cquad_space(int nt,int nr,double R,double dr0,double dt0,point *p,char bcs[4][2][4]){
  
  int i,j,k,nvert0=nvert;
  double theta,thetb,thet0,delta,rr;
  double *ta,*tb;
  point pc,r0,r1;
  connector r,r2,r3,r4;

  printf("Generating %d elements in a curved quadrilateral space with geometric growth\n",nt*nr);
  
//for(k=0;k<nfld;k++) for(i=0;i<4;i++) printf("\t\"%s\"\n",bcs[i][k]);
  
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
      for(k=0;k<nfld;k++){
        if(i==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
        if(j==0){ //j is build from side 4 [3][k] to 2 [1][k]
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(j==nt-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
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

