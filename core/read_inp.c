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

#define VOPS 2
#define MXTOKL 16

  char line[MXLS];
  char ops[VOPS][MXTOKL]={"list\0","mid\0"};
//char op[MXTOKL];
  char *lfile[256];
  char *tok[MXTOKL];
  char *cpt;
  
  int i,j,noop,ipt=0;
  int nops[VOPS]={0,0},iop[VOPS]={0,0};

  list *lists;
  mid *mids;

// count the valid operations
  printf("Counting point operations...\n");
  while((fgets(line,MXLS,fp)!=NULL)&&(strcmp(line,"[points]\n")!=0));
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
  points=malloc(npts*sizeof(point));
  lists =malloc(nops[0]*sizeof(list));
  mids  =malloc(nops[1]*sizeof(mid));

// read the operations
  rewind(fp);
  while((fgets(line,MXLS,fp)!=NULL)&&(strcmp(line,"[points]\n")!=0));
  while((fgets(line,MXLS,fp)!=NULL)&&(strcmp(line,"\n")!=0)){
    cpt=strtok(line," ");
    i=0;
    while(cpt != NULL){
      tok[i++]=cpt;
      cpt=strtok(NULL," "); 
    }
    for(i=0;i<VOPS;i++){
      if(strncmp(tok[0],ops[i],strlen(ops[i]))==0){
        if(i==LIST) {
          (lists+iop[i])->npts=atoi(tok[1]);
          strcpy((lists+iop[i])->fname,tok[2]);
          iop[i]++;
        }else if(i==MID){
          (mids+iop[i])->pt1=atoi(tok[1]);
          (mids+iop[i])->pt2=atoi(tok[2]);
          iop[i]++;
        }
      }
    }   
  } 

// do the operations
  printf("Generating point data...\n");
  for(i=0;i<VOPS;i++){ for(j=0;j<nops[i];j++){
    if(i==LIST) read_list(lists[j],&ipt);
    else if(i==MID){
      printf("Generating mid point between points %d and %d\n",(mids+j)->pt1,(mids+j)->pt2);
      points[ipt]=midpoint(points[(mids+j)->pt1],points[(mids+j)->pt2]);
      ipt++;
    }
  }}
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
