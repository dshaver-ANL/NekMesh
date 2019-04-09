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

#define LIST 0
#define MID 1
#define LIN 2
#define ROT 3

#define VOPS 4
#define MXTOKL 16

  char line[MXLS];
  char ops[VOPS][MXTOKL]={"list\0","mid\0","linear\0","rotate\0"};
//char op[MXTOKL];
  char *lfile[256];
  char *tok[MXTOKL];
  char *cpt;

  double r1;
  
  int i,j,noop,ipt=0,pt1,pt2,pt3;
  int nops[VOPS],iop[VOPS],tops=0;
  int *optype;

  list *lists;
  mid *mids;
  lin *lins;
  rot *rots;

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
          strcpy((lists+iop[i])->fname,tok[2]);
          iop[i]++;
        }else if(i==MID){
          (mids+iop[i])->pt1=atoi(tok[1]);
          (mids+iop[i])->pt2=atoi(tok[2]);
          iop[i]++;
        }else if(i==LIN){
          (lins+iop[i])->pt1=atoi(tok[1]);
          (lins+iop[i])->pt2=atoi(tok[2]);
          (lins+iop[i])->fra=atof(tok[3]);
          iop[i]++;
        }else if(i==ROT){
          (rots+iop[i])->pt1=atoi(tok[1]);
          (rots+iop[i])->ang=atof(tok[2]);
          (rots+iop[i])->org=atoi(tok[3]);
          iop[i]++;
        }
      }
    }   
  } 

// do the operations IN ORDER!
  printf("Generating point data...\n");
  for(i=0;i<VOPS;i++) iop[i]=0;
  for(i=0;i<tops;i++){ j=iop[optype[i]]; switch(optype[i]){
    case LIST:
      read_list(lists[iop[LIST]],&ipt);
      break;
    case MID:
      pt1=(mids+j)->pt1;
      pt2=(mids+j)->pt2;
      printf("point %d: Generating mid point between points %d and %d\n",ipt+1,pt1,pt2);
      points[ipt]=midpoint(points[pt1-1],points[pt2-1]);
      ipt++;
      break;
    case LIN:
      pt1=(lins+j)->pt1;
      pt2=(lins+j)->pt2;
      r1=(lins+j)->fra;
      printf("point %d: Generating point %.2f\% between points %d and %d\n",ipt+1,r1*100.0,pt1,pt2);
      points[ipt]=linpoint(r1,points[pt1-1],points[pt2-1]);
      ipt++;
      break;
    case ROT:
      pt1=(rots+j)->pt1;
      r1=(rots+j)->ang;
      pt2=(rots+j)->org;
      if(pt2==0){
        printf("point %d: Rotating point %d %.2f degrees about the origin\n",ipt+1,pt1,r1);
        r1*=M_PI/180.0;
        points[ipt]=rotate_point(points[pt1-1],r1,origin);
      }else{
        printf("point %d: Rotating point %d %.2f degrees about point %d\n",ipt+1,pt1,r1,pt2);
        r1*=M_PI/180.0;
        points[ipt]=rotate_point(points[pt1-1],r1,points[pt2-1]);
      }
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
      printf("\tError: only %d points in file \"%s\"!\n",i,l0.fname);
      return -1;
    }
    sscanf(line,"%lf %lf",&x,&y);
    (points+*ipt)->x=x;(points+*ipt)->y=y;
    (*ipt)++;
  }

return 0;
}
/*----------------------------------------------------*/
