int write_rea(char *fname){

  char line[256],bcout[4];

  FILE *reain,*reaout;

  int i,j,k,ivert;
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

  reain=fopen("core/reatop.txt","r");
  while(fgets(line,256,reain)!=NULL) fprintf(reaout,"%s",line);
  fclose(reain);

  fprintf(reaout,"           13  LOGICAL SWITCHES FOLLOW\n");
  fprintf(reaout,"  T     IFFLOW\n");
  if(nfld==1){ 
    fprintf(reaout,"  F     IFHEAT\n");
  }else{
    fprintf(reaout,"  T     IFHEAT\n");
  }
  fprintf(reaout,"  T     IFTRAN\n");
  fprintf(reaout,"  T T F F F F F F F F F IFNAV & IFADVC (convection in P.S. fields)\n");
  fprintf(reaout,"  F F T T T T T T T T T T IFTMSH (IF mesh for this field is T mesh)\n");
  fprintf(reaout,"  F     IFAXIS\n");
  fprintf(reaout,"  F     IFSTRS\n");
  fprintf(reaout,"  F     IFSPLIT\n");
  fprintf(reaout,"  F     IFMGRID\n");
  fprintf(reaout,"  F     IFMODEL\n");
  fprintf(reaout,"  F     IFKEPS\n");
  fprintf(reaout,"  F     IFMVBD\n");
  fprintf(reaout,"  F     IFCHAR\n");

  fprintf(reaout,"%10.6f %10.6f %10.6f %10.6f   XFAC,YFAC,XZERO,YZERO\n",1.2*dx,1.2*dx,xmin-0.1*dx,ymin-0.1*dx);
  fprintf(reaout,"**MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y.\n");

  fprintf(reaout,"    %5d         2      %5d NEL,NDIM,NELV\n",nelem,nelem);
  for(i=0;i<nelem;i++){
    fprintf(reaout,"            ELEMENT%5d [    1a]    GROUP     0\n",i+1);
    for(j=0;j<4;j++){
      ivert=(elems+i)->vid[j];
      fprintf(reaout,"%14.6e",(verts+ivert)->x);
    }
    fprintf(reaout,"\n");
    for(j=0;j<4;j++){
      ivert=(elems+i)->vid[j];
      fprintf(reaout,"%14.6e",(verts+ivert)->y);
    }
    fprintf(reaout,"\n");
  }

  fprintf(reaout,"  ***** CURVED SIDE DATA *****\n");
  if(ncide==0){
    fprintf(reaout,"0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n");
  }else{
    fprintf(reaout,"%6d Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n",ncide);
    if(nelem<1000){for(i=0;i<ncide;i++)
      fprintf(reaout,"%3d%3d%14.6f%14.6g%14.6g%14.6g%14.6g %c\n",
        (cides+i)->esid,(cides+i)->elid+1,(cides+i)->curve,0.0,0.0,0.0,0.0,(cides+i)->ccurve);
    }else if(nelem<1000000){ for(i=0;i<ncide;i++)
      fprintf(reaout,"%2d%6d%14.6f%14.6g%14.6g%14.6g%14.6g %c\n",
        (cides+i)->esid,(cides+i)->elid+1,(cides+i)->curve,0.0,0.0,0.0,0.0,(cides+i)->ccurve);
    }else printf("Error in output: format for nel > 1e6 not supported!\n");
  }

  fprintf(reaout,"  ***** BOUNDARY CONDITIONS *****\n");
  fprintf(reaout,"  ***** FLUID   BOUNDARY CONDITIONS *****\n");
  k=0;
  for(i=0;i<nelem;i++){
    for(j=0;j<4;j++){
      strcpy(bcout,(elems+i)->BC[j][k]);
      if(nelem<999){
        fprintf(reaout," %s%3d  %1d %13.5f %13.5f %13.5f %13.5f %13.5f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
      }else{
        fprintf(reaout," %s%5d%1d %13.5f %13.5f %13.5f %13.5f %13.5f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
      }
    }
  }
  if(nfld==1){
    fprintf(reaout,"  ***** NO THERMAL BOUNDARY CONDITIONS *****\n");
  }else{
    fprintf(reaout,"  ***** THERMAL BOUNDARY CONDITIONS *****\n");
    k=1;
    for(i=0;i<nelem;i++){
      for(j=0;j<4;j++){
        strcpy(bcout,(elems+i)->BC[j][k]);
        if(nelem<999){
          fprintf(reaout," %s%3d  %1d %13.5f %13.5f %13.5f %13.5f %13.5f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
        }else{
          fprintf(reaout," %s%5d%1d %13.5f %13.5f %13.5f %13.5f %13.5f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
        }
      }
    }
  }

  reain=fopen("core/reabot.txt","r");
  while(fgets(line,256,reain)!=NULL) fprintf(reaout,"%s",line);
  fclose(reain);

  fclose(reaout);

return 0;
}
//-------------------------------------------------------------------
int reset(void){
  printf("\n\treseting mesh generator\n\n");
  nelem = 0;
  nvert = 0;
  ncide = 0;
}
//-------------------------------------------------------------------
int output_pts(point *pts,int npts,char *fnm){

  int i;
  FILE *ptout;

  ptout=fopen(fnm,"w");
  for(i=0;i<npts;i++) fprintf(ptout,"%13.7f %13.7f \"%d\"\n",(pts+i)->x,(pts+i)->y,i);
  fclose(ptout);

  return i;
}
//-------------------------------------------------------------------

