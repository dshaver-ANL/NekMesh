#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MXLS 512
/*----------------------------------------------------*/
int read_inp(char *fname){

  FILE *inp;

  inp=fopen(fname,"r");

  if(inp==NULL){
    printf("Unable to open file \"%s\"!\n",fname);
    return -1;
  }

  npts=read_points(inp);
//read_spaces(inp);

  fclose(inp);

return 0;
}
/*----------------------------------------------------*/
int read_points(FILE *fp){

#define LIST 0  //read list
#define MID 1   //midpoint
#define LIN 2   //linear distance between two points
#define ROT 3   //rotation
#define LCI 4   //line-circle intercept
#define LLI 5   //line-line intercept
#define REF 6   //reflection
#define TRN 7   //translation
#define APT 8   //add point

#define VOPS 9
#define MXTOKL 16

  char line[MXLS];
  char ops[VOPS][MXTOKL]={"list\0","mid\0","linear\0","rotate\0","intercept-lc\0","intercept-ll\0","reflect\0","translate\0","point\0"};
//char op[MXTOKL];
  char *lfile[256];
  char *tok[MXTOKL];
  char *cpt;

  double r1,r2;
  
  int i,j,k,noop,ipt=0,pt1,pt2,pt3,pt4;
  int nops[VOPS],iop[VOPS],tops=0;
  int *optype;

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

  for(i=0;i<VOPS;i++){
    nops[i]=0;
    iop[i]=0;
  }
// count the valid operations
  printf("Counting point operations...\n");
  while((fgets(line,MXLS,fp)!=NULL)&&(noop=strncmp(line,"[points]",8)!=0));
  if(noop!=0){
    printf("No point operations found!\n");
    return 0;
  }
  npts=1; //origin is always pt 0
  while((fgets(line,MXLS,fp)!=NULL)&&(strcmp(line,"\n")!=0)){
    noop=1;
    cpt=strtok(line," ");
    i=0;
    while(cpt != NULL){
      tok[i++]=cpt;
      cpt=strtok(NULL," "); 
    }
    for(i=0;i<VOPS;i++){
      if(strncmp(tok[0],ops[i],strlen(ops[i]))==0){
        tops++;
        nops[i]++;
        noop=0;
        if(i==LIST) {
          npts+=atoi(tok[1]);
        }else npts++;
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

// read the operations
  rewind(fp);j=0;
  while((fgets(line,MXLS,fp)!=NULL)&&(strncmp(line,"[points]",8)!=0));
  while((fgets(line,MXLS,fp)!=NULL)&&(strcmp(line,"\n")!=0)){
    cpt=strtok(line," ");
    i=0;
    while(cpt != NULL){
      tok[i++]=cpt;
      cpt=strtok(NULL," "); 
    }
    for(i=0;i<VOPS;i++){
      if(strncmp(tok[0],ops[i],strlen(ops[i]))==0){
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
          (lins+iop[i])->fra=atof(tok[3]);
        }else if(i==ROT){
          (rots+iop[i])->pt1=atoi(tok[1]);
          (rots+iop[i])->ang=atof(tok[2]);
          (rots+iop[i])->org=atoi(tok[3]);
        }else if(i==LCI){
          (lcis+iop[i])->pt1=atoi(tok[1]);
          (lcis+iop[i])->pt2=atoi(tok[2]);
          (lcis+iop[i])->org=atoi(tok[3]);
          (lcis+iop[i])->rad=atof(tok[4]);
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
          (trns+iop[i])->dxt=atof(tok[2]);
          (trns+iop[i])->dyt=atof(tok[3]);
        }else if(i==APT){
          (apts+iop[i])->xc=atof(tok[1]);
          (apts+iop[i])->yc=atof(tok[2]);
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
      printf("point %d: Creating point %.2f\% between points %d and %d\n",ipt,r1*100.0,pt1,pt2);
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
      printf("point %d: Creating point at %.2f, %.2f\n",ipt,r1,r2);
      points[ipt].x=r1;
      points[ipt].y=r2;
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
