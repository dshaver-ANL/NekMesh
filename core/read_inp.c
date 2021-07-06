#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define max(x,y) ((x) > (y) ? (x) : (y))

#define MXLS 512
#define MXTOKL 16
/*----------------------------------------------------*/
int read_inp(char *fname){

  FILE *inp;
  var *vars=NULL;
  int nvars=0;

  inp=fopen(fname,"r");

  if(inp==NULL){
    printf("Unable to open file \"%s\"!\n",fname);
    return -1;
  }

  read_vars(inp,&vars,&nvars);
  rewind(inp);
  npts=read_points(inp,vars,nvars);
//read_spaces(inp);

  fclose(inp);

  printf("Done reading input\n");

return 0;
}
/*----------------------------------------------------*/
int read_vars(FILE *fp,var **vars,int *nvars){

  char line[MXLS],*tok[MXTOKL],*cpt;
  int i,j,novar;

// count the number of defined constants
  for(i=0;i<MXTOKL;i++) tok[i]=NULL;
  printf("Counting defined variables...\n"); 
  while((fgets(line,MXLS,fp)!=NULL)&&(novar=strncmp(line,"[variables]",11)!=0));
  if(novar!=0){
    printf("No defined variables found\n");
    return 0;
  }
  while((fgets(line,MXLS,fp)!=NULL)&&(strcmp(line,"\n")!=0)){
    novar=1;
    cpt=strtok(line," ");
    i=0;
    while((cpt!=NULL)&&(i<MXTOKL)){
      tok[i++]=cpt; 
      cpt=strtok(NULL," ");
    }
    if(isalpha(*tok[0])&&(isdigit(*tok[1])||(*tok[1]=='-'))){
      novar=0;
      (*nvars)++;
    }
    if(novar) printf("Invalid variable name: \"%s\", ignoring\n",tok[0]);
  }
  printf("...Done!\n\n");

  if(*nvars>0){
    *vars=malloc(sizeof(var)* *nvars);
    rewind(fp);

    printf("Reading %d variables...\n",*nvars);
    while((fgets(line,MXLS,fp)!=NULL)&&(novar=strncmp(line,"[variables]",11)!=0));
    for(j=0;j<*nvars;j++){
      fgets(line,MXLS,fp);
      cpt=strtok(line," ");
      i=0;
      while((cpt!=NULL)&&(i<MXTOKL)){
        tok[i++]=cpt; 
        cpt=strtok(NULL," ");
      }
      strcpy((*vars)[j].name,tok[0]);
      (*vars)[j].val=atof(tok[1]);
      printf("\"%s\" = %g\n",(*vars)[j].name,(*vars)[j].val);
    }
    printf("\n");
  }

return 0;
}
/*----------------------------------------------------*/
int read_points(FILE *fp,var *vars,int nvars){

#define LIST 0  //read list
#define MID 1   //midpoint
#define LIN 2   //linear distance between two points
#define ROT 3   //rotation
#define LCI 4   //line-circle intercept
#define LLI 5   //line-line intercept
#define REF 6   //reflection
#define TRN 7   //translation
#define APT 8   //add point
#define RPT 9   //add point r-theta

#define VOPS 10

  char line[MXLS],*tok[MXTOKL],*cpt;
  char ops[VOPS][MXTOKL]={"list\0","mid\0","linear\0","rotate\0","intercept-lc\0","intercept-ll\0","reflect\0","translate\0","point\0","point-rt\0"};
  char *lfile[256];

  double r1,r2;
  
  int i,j,k,noop=1,ipt=0,pt1,pt2,pt3,pt4;
  int nops[VOPS],iop[VOPS],tops=0;
  int *optype,tlen,olen;

  point p1;

  list *lists;
  mid *mids;
  lin *lins;
  rot *rots;
  lci *lcis;
  lli *llis;
  ref *refs;
  trn *trns;
  apt *apts;
  rpt *rpts;

  for(i=0;i<VOPS;i++){
    nops[i]=0;
    iop[i]=0;
  }
  for(i=0;i<MXTOKL;i++) tok[i]=NULL;
// count the valid operations
  printf("Counting point operations...\n");
  while((fgets(line,MXLS,fp)!=NULL)&&(noop=strncmp(line,"[points]",8)!=0));
  if(noop!=0){
    printf("No point operations found!\n");
    return 0;
  }
  npts=1; //origin is always pt 0
  while((fgets(line,MXLS,fp)!=NULL)&&(strcmp(line,"\n")!=0)){
//  printf("line = \"%s\"\n",line);
    noop=1;
    cpt=strtok(line," ");
    i=0;
    while(cpt != NULL){
      tok[i++]=cpt;
      cpt=strtok(NULL," "); 
    }
    tlen=strlen(tok[0]);
    for(i=0;i<VOPS;i++){
      olen=strlen(ops[i]);
      if(strncmp(tok[0],ops[i],max(tlen,olen))==0){
        tops++;
        nops[i]++;
        noop=0;
        if(i==LIST) {
          npts+=atoi(tok[1]);
        }else npts++;
        i=VOPS;
      }
    }
    if(noop) printf("Invalid point operation: \"%s\", ignoring\n",tok[0]);
  }
  printf("...Done!\n\n");

// allocate memory
  optype=malloc(tops*sizeof(int));
  points=malloc(npts*sizeof(point));
  lists =malloc(nops[LIST]*sizeof(list));
  mids  =malloc(nops[MID]*sizeof(mid));
  lins  =malloc(nops[LIN]*sizeof(lin));
  rots  =malloc(nops[ROT]*sizeof(rot));
  lcis  =malloc(nops[LCI]*sizeof(lci));
  llis  =malloc(nops[LLI]*sizeof(lli));
  refs  =malloc(nops[REF]*sizeof(ref));
  trns  =malloc(nops[TRN]*sizeof(trn));
  apts  =malloc(nops[APT]*sizeof(apt));
  rpts  =malloc(nops[RPT]*sizeof(rpt));

// read the operations
  rewind(fp);j=0;
  while((fgets(line,MXLS,fp)!=NULL)&&(strncmp(line,"[points]",8)!=0));
  while((fgets(line,MXLS,fp)!=NULL)&&(strcmp(line,"\n")!=0)){
    sanitize_string(line,MXLS);
    cpt=strtok(line," ");
    i=0;
    while(cpt != NULL){
      tok[i++]=cpt;
      cpt=strtok(NULL," "); 
    }
    tlen=strlen(tok[0]);
    for(i=0;i<VOPS;i++){
      olen=strlen(ops[i]);
      if(strncmp(tok[0],ops[i],max(tlen,olen))==0){
        optype[j++]=i;
        if(i==LIST) {
          (lists+iop[i])->npts=atoi(tok[1]);
          k=strlen(tok[2])-1;
          if(tok[2][k]=='\n') tok[2][k]='\0';
          strcpy((lists+iop[i])->fname,tok[2]);
        }else if(i==MID){
          (mids+iop[i])->pt1=atoi(tok[1]);
          (mids+iop[i])->pt2=atoi(tok[2]);
        }else if(i==LIN){
          (lins+iop[i])->pt1=atoi(tok[1]);
          (lins+iop[i])->pt2=atoi(tok[2]);
          (lins+iop[i])->fra=get_var(tok[3],vars,nvars);
        }else if(i==ROT){
          (rots+iop[i])->pt1=atoi(tok[1]);
          (rots+iop[i])->ang=get_var(tok[2],vars,nvars);
          (rots+iop[i])->org=atoi(tok[3]);
        }else if(i==LCI){
          (lcis+iop[i])->pt1=atoi(tok[1]);
          (lcis+iop[i])->pt2=atoi(tok[2]);
          (lcis+iop[i])->org=atoi(tok[3]);
          (lcis+iop[i])->rad=get_var(tok[4],vars,nvars);
        }else if(i==LLI){
          (llis+iop[i])->pt1=atoi(tok[1]);
          (llis+iop[i])->pt2=atoi(tok[2]);
          (llis+iop[i])->pt3=atoi(tok[3]);
          (llis+iop[i])->pt4=atoi(tok[4]);
        }else if(i==REF){
          (refs+iop[i])->pt1=atoi(tok[1]);
          (refs+iop[i])->pt2=atoi(tok[2]);
          (refs+iop[i])->pt3=atoi(tok[3]);
        }else if(i==TRN){
          (trns+iop[i])->pt1=atoi(tok[1]);
          (trns+iop[i])->dxt=get_var(tok[2],vars,nvars);
          (trns+iop[i])->dyt=get_var(tok[3],vars,nvars);
        }else if(i==APT){
          (apts+iop[i])->xc=get_var(tok[1],vars,nvars);
          (apts+iop[i])->yc=get_var(tok[2],vars,nvars);
        }else if(i==RPT){
          (rpts+iop[i])->rc=get_var(tok[1],vars,nvars);
          (rpts+iop[i])->tc=get_var(tok[2],vars,nvars);
        }
        iop[i]++;
      }
    }   
  } 

// do the operations IN ORDER!
  printf("Creating point data...\n");
  points[0]=origin;ipt++; //origin is always pt 0
  for(i=0;i<VOPS;i++) iop[i]=0;
  for(i=0;i<tops;i++){ j=iop[optype[i]]; switch(optype[i]){
    case LIST:
      read_list(lists[iop[LIST]],&ipt);
      break;
    case MID:
      pt1=(mids+j)->pt1;
      pt2=(mids+j)->pt2;
      printf("point %d: Creating mid point between points %d and %d\n",ipt,pt1,pt2);
      points[ipt]=midpoint(points[pt1],points[pt2]);
      ipt++;
      break;
    case LIN:
      pt1=(lins+j)->pt1;
      pt2=(lins+j)->pt2;
      r1=(lins+j)->fra;
      printf("point %d: Creating point %.2f%% between points %d and %d\n",ipt,r1*100.0,pt1,pt2);
      points[ipt]=linpoint(r1,points[pt1],points[pt2]);
      ipt++;
      break;
    case ROT:
      pt1=(rots+j)->pt1;
      r1=(rots+j)->ang;
      pt2=(rots+j)->org;
      printf("point %d: Rotating point %d %.2f degrees about point %d\n",ipt,pt1,r1,pt2);
      r1*=M_PI/180.0;
      points[ipt]=rotate_point(points[pt1],r1,points[pt2]);
      ipt++;
      break;
    case LCI:
      pt1=(lcis+j)->pt1;
      pt2=(lcis+j)->pt2;
      pt3=(lcis+j)->org;
      r1=(lcis+j)->rad;
      printf("point %d: Creating point at intercept of the line formed from points %d and %d with a circle centered at point %d with radius %.2f\n",ipt,pt1,pt2,pt3,r1);
      points[ipt]=line_circle_intercept(points[pt1],points[pt2],points[pt3],r1);
      ipt++;
      break;
    case LLI:
      pt1=(llis+j)->pt1;
      pt2=(llis+j)->pt2;
      pt3=(llis+j)->pt3;
      pt4=(llis+j)->pt4;
      printf("point %d: Creating point at intercept of the lines formed from points %d and %d with points %d and %d\n",ipt,pt1,pt2,pt3,pt4);
      points[ipt]=line_line_intercept(points[pt1],points[pt2],points[pt3],points[pt4]);
      ipt++;
      break;
    case REF:
      pt1=(refs+j)->pt1;
      pt2=(refs+j)->pt2;
      pt3=(refs+j)->pt3; 
      printf("point %d: Reflecting point %d over the line formed from points %d and %d\n",ipt,pt1,pt2,pt3);
      points[ipt]=reflect_point(points[pt1],points[pt2],points[pt3]);
      ipt++;
      break;
    case TRN:
      pt1=(trns+j)->pt1;
      r1=(trns+j)->dxt;
      r2=(trns+j)->dyt;
      printf("point %d: Translating point %d by %.2f, %.2f\n",ipt,pt1,r1,r2);
      p1.x=r1;p1.y=r2;
      points[ipt]=translate_point(points[pt1],p1);
      ipt++;
      break;
    case APT:
      r1=(apts+j)->xc;
      r2=(apts+j)->yc;
      printf("point %d: Creating point at %.2f, %.2f (x ,y)\n",ipt,r1,r2);
      points[ipt].x=r1;
      points[ipt].y=r2;
      ipt++;
      break;
    case RPT:
      r1=(rpts+j)->rc;
      r2=(rpts+j)->tc;
      printf("point %d: Creating point at %.2f, %.2f (r, theta)\n",ipt,r1,r2);
      points[ipt].x=r1*cos(r2*M_PI/180.);
      points[ipt].y=r1*sin(r2*M_PI/180.);
      ipt++;
      break;

  } iop[optype[i]]++; }
  printf("...Done!\n\n");

return npts;
}
/*----------------------------------------------------*/
int read_list(list l0,int *ipt){

  char line[256];
  double x,y;
  FILE *fp;
  int i;

  printf("Reading list of %d points from file \"%s\"\n",l0.npts,l0.fname);
  fp=fopen(l0.fname,"r");

  if(fp==NULL){
    printf("\tError opening file \"%s\"!\n",l0.fname);
    return -1;
  }

  for(i=0;i<l0.npts;i++){
    if(fgets(line,256,fp)==NULL){
      printf("\tWarning: only %d points in file \"%s\"!\n",i,l0.fname);
      return -1;
    }
    sscanf(line,"%lf %lf",&x,&y);
    (points+*ipt)->x=x;(points+*ipt)->y=y;
    (*ipt)++;
  }

return 0;
}
/*----------------------------------------------------*/
double get_var(char *tok,var *vars,int nvars){

  int i,tlen,vlen;
  double val;

  if(*tok=='-'){
    return -1.0*get_var(&tok[1],vars,nvars);
  }else{
    if(isdigit(*tok)){ 
      return atof(tok);
    }else{
      for(i=0;i<nvars;i++){
        if(strcmp(tok,vars[i].name)==0){
          return vars[i].val;
        }
      }
    }
  } 

return 0.0;
}
/*----------------------------------------------------*/
int sanitize_string(char *str,int len){

  int i;

  for(i=0;i<len;i++){
    if(str[i]=='\n'){
      str[i]=' ';
      return i;
    }
  }

return 0;
}
