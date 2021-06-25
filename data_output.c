#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "constants.h"
#include "my.h"
#include "parameters.h"
#include "data_output.h"
#include "mersenne.h"


// //for now it is really just dumb printing, later maybe something better
// void PrintACell(int Time,CELL *ic)
// {
//   double R=*(ic->Rp) + *(ic->Rr)+ *(ic->R)*2. + *(ic->Crmt)*2. + *(ic->Crme)*2. + *(ic->Crmq)*2. + *(ic->Crmrp)*2.;
//   double Q=*(ic->Q);
//   double T=*(ic->T);// + *(ic->Cst) + *(ic->Cct);
//   double E=*(ic->E);// + *(ic->Cce) + *(ic->Cae);
//   double cellvol=ic->cellvol;
//   
//   fprintf(stderr,"%d %f %f %f %f %f %f %f %f %f %f %f %f\n", Time, cellvol, R,Q,T,E,
// 				    *(ic->mRp),*(ic->mQ),*(ic->mT),*(ic->mE), 
// 				    *(ic->S), *(ic->C),*(ic->A) );
// }

//

void OutputBackup(int Time, CELL *world, CELL *prev_world)
{
  char fname[MAXLEN];
  char command[MAXLEN];
  
  int printTime=Time-par_time_anc;
  
  sprintf(command,"%s %s","mkdir -p",dirname_backup);
  if(system(command)==-1){
    fprintf(stderr,"OutputBackup: Failed to mkdir %s. Save here\n",dirname_backup);
    sprintf(fname,"%s_t%d",file_backup,printTime);
  }
  else
    sprintf(fname,"%s/%s_t%d",dirname_backup,file_backup,printTime);
  
  SaveBackup(printTime,world,fname,prev_world);
}

//This writes backup files, lots of info is written here
// but this does also anc trace
void SaveBackup(int printTime, CELL *world, char *fname, CELL *prev_world)
{
  int i,j;
  CELL *ic;
  FILE *fp;
  
  int ancpos[par_pop];
  for(i=0;i<par_pop;i++) ancpos[i]=0;	//initialise ancpos table
  
  fp=fopen(fname,"w");
  if(fp==NULL) fprintf(stderr,"SaveBackup(): Warning. File not opened.\n");
  fprintf(fp,"Time %d par_pop %d kt alfat betat ktrscrT ktrscrQ ktrscrRr\n",printTime, par_pop);
  
  for(i=0;i<par_pop;i++){
    //go through current world and retrieve the position of the ancestors...
    int this_ancpos= world[i].anc;
    ancpos[this_ancpos]=1; //the position this_ancpos in the array ancpos is now 1
  }
  
  //Now ancpos contains 1 if that position is ancestral, go take that from prev_world and write it. 
  for(i=0;i<par_pop;i++){
    if(printTime%par_time_backup!=0){ 
      if(ancpos[i]==0) continue;
      //else save only ancestors
    }
    //else save everybody regardless of ancestry
    
    ic = &prev_world[i];
    
    fprintf(fp,"%d %d",i, ic->anc); 
    for(j=0;j<2;j++ ) fprintf(fp," %f",ic->ktrscrT[j]);
    for(j=0;j<2;j++ ) fprintf(fp," %f",ic->ktrscrQ[j]);
    for(j=0;j<2;j++ ) fprintf(fp," %f",ic->ktrscrRr[j]);
    for(j=0;j<2;j++ ) fprintf(fp," %f",ic->ktrscrRp[j]);
    
    //for(j=0;j<3;j++) fprintf(fp," %f",ic->mut_genes[j]);
    for(j=0;j<3;j++) fprintf(fp," %f",ic->max_ktrmut[j]);   //in this model!
    
    fprintf(fp," %f %f %f",ic->health, *(ic->S) , *(ic->A) ); // for now this
    fprintf(fp," %f %f",ic->av_tr_mut, ic->lifespan);
    
    fprintf(fp," %s",ic->genome);
    fprintf(fp,"\n");
  }
  
   
  fclose(fp);
}


//EXactly like SaveDataRelative3, but saves most abundant genotype in the pop
void SaveDataRelative4(double Time,CELL *world)
{
  int i,j,loopflag=0; 
  int bestguy; //at the beginning everyone is the same, so it doesn't matter
  //char pop_genomes[par_pop][MAXGENOMESIZE]; // this segfaults because it's too big for the stack 
  int ngenomes[par_pop];
  
  char **pop_genomes;
  pop_genomes= malloc(par_pop*sizeof(char*));
  for(i=0;i<par_pop;i++) pop_genomes[i]=malloc(MAXGENOMESIZE*sizeof(char));
  
  //pop_genomes contains unique genomes in the first position where they are found
  //ngenomes says, at that position, how many times that genome is found in the pop
  for(i=0;i<par_pop;i++){
    ngenomes[i]=0; //initialise ngenomes array
    pop_genomes[i][0]='\0';
  }
  strcpy(pop_genomes[0],world[0].genome); //we copy the first genome of the population
  ngenomes[0]=1;
  int ndiffgenomes=1;
  //printf("\n\n\n\n");
  
  for(i=1;i<par_pop;i++){
    //printf("i: %d,   genome: %s\n",i,world[i].genome);
    for(j=0;j<i;j++){
      //printf("  j: %d, genome: %s\n",j,pop_genomes[j] ,ndiffgenomes);
      
      if(world[i].health>0 && strcmp( world[i].genome, pop_genomes[j])==0 ){
	ngenomes[j]++;
	loopflag=1;
	break; // we found it, let's move on
      }
    }
    if(loopflag==1){
      loopflag=0;
      continue;
    }
    //if you didn't break before, you found a new genome
    if(pop_genomes[i][0]!='\0') fprintf(stderr,"Warning, overwriting some genome\n");
    //else fprintf(stderr,"New genome found\n");
    
    strcpy(pop_genomes[ndiffgenomes],world[i].genome);
    
    ngenomes[i]++;
    ndiffgenomes++;  
  }
  
  //in this second loop we retrieve the first guy whose genome is most abundant
  int pos=-1;
  int abundance=-1;
  double best_health=-1.;
  int this=0;
  for(i=0;i<par_pop;i++){
    if(world[i].health>best_health && ngenomes[i]>abundance){
      pos=i;
      abundance=ngenomes[i];
      this=1;
    }
  }
  if(this) bestguy=pos;
  else bestguy=par_pop*genrand_real2();
  
  bestguy=0;
  int ancestors=0;
  int anc_array[par_pop];
  
  
  for(i=0;i<par_pop;i++) anc_array[i]=0; //initialise anc array
  
  for(i=0;i<par_pop;i++){
    //here we are going to have a look at the diversity in the population
    anc_array[ world[i].id ]++;
    
    //if( world[i].health > 1.01*world[bestguy].health){
    //  bestguy=i;
      //printf("bestguy health = %f\n",besthealth );
    //}
  }
    
  for(i=0;i<par_pop;i++) if(anc_array[i]!=0) ancestors++;
  
  if(ancestors==1){
    for(i=0;i<par_pop;i++) world[i].id=i;
  }
  
  
  
  //printf("Hello save 0.1\n");
  //bestguy=par_pop*genrand_real2();
  //printf("printing pos: %d\n",bestguy);
  CELL *ic=&world[bestguy];
  FILE *fp;
  
  double cellvol=ic->cellvol;
  double tarvol=ic->tarvol;
  double health=ic->health;
  //double R=(*(ic->Rp) + *(ic->R) + *(ic->Crmt) + *(ic->Crme) + *(ic->Crmq) + *(ic->Crmrp))/cellvol;
  //double R1=(*(ic->Rr) + *(ic->Crmt) + *(ic->Crmq) )/cellvol;
  double Rp=(*(ic->Rp) )/cellvol;
  double Rr=(*(ic->Rr) )/cellvol;
  double Q=*(ic->Q)/cellvol;
  double T=*(ic->T)/cellvol;// + *(ic->Cst) + *(ic->Cct);
  
  fp=fopen(file_out_var,"a");
  if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
  
  fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d %f\n", 
				    Time, health, cellvol 
				    , Rp, Rr, Q , T , *(ic->mQ) , *(ic->mT) , *(ic->mRp)
				    ,*(ic->S) , *(ic->A)
				    , ic->growthrate, tarvol ,ic->lifespan
				    , ic->id , bestguy , ancestors , abundance , ndiffgenomes
				    , ic->av_tr_mut
	 );
  
  fclose(fp);	
  
  fp=fopen(file_out_par,"a");
  if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
  
//   fprintf(fp,"%f %f %f %f %f %f %#include <float.h>f %f %f %f %f\n",
// 	            Time, ic->ks1,ic->ks2,ic->ks3,ic->ksm1,ic->ksm3,
// 	            ic->kc1,ic->kc2,ic->kc3,ic->kcm1,ic->kcm3);
  
  //If you really want to print stuff you do not mutate, 
  // then you shoould also print sM,sR,sA
  
  fprintf(fp,"%f %f %f %f",
 	            Time, ic->kt,ic->alfat,ic->betat
	//	          ,ic->ke,ic->alfae,ic->betae
	 );
  for(i=0;i<2;i++ ) fprintf(fp," %f",ic->ktrscrT[i]);
  //for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrE[i]);
  for(i=0;i<2;i++ ) fprintf(fp," %f",ic->ktrscrQ[i]);
  for(i=0;i<2;i++ ) fprintf(fp," %f",ic->ktrscrRr[i]);
  for(i=0;i<2;i++ ) fprintf(fp," %f",ic->ktrscrRp[i]);
  
  for(i=0;i<N_GENE_TYPES;i++ ) fprintf(fp," %d",ic->genesum[i]);
  for(i=0;i<N_GENE_TYPES;i++ ) fprintf(fp," %d",ic->delet_genesum[i]);
  //for(i=0;i<3;i++) fprintf(fp," %f",ic->mut_genes[i]);
  for(j=0;j<3;j++) fprintf(fp," %f",ic->max_ktrmut[j]);   //in this model!
  
  fprintf(fp," %s",ic->genome);
  fprintf(fp,"\n");
  fclose(fp);
  
  for(i=0;i<par_pop;i++) free(pop_genomes[i]);
  free(pop_genomes);
}




//EXactly like SaveDataRelative2, but saves fittest guy
// (notice that fittest guy may not be ancestor)
//save a cell, but al vars are relative to volume
// and volume is scaled by target volume
// version 2 splits Rrna contribution from Rproteins
/*
void SaveDataRelative3(int Time,CELL *world)
{
  int i; 
  static int bestguy=0; //at the beginning everyone is the same, so it doesn't matter
  //double besthealth=-1.;
  int anc_array[par_pop];
  int ancestors=0;
  
  for(i=0;i<par_pop;i++) anc_array[i]=0; //initialise anc array
  
  for(i=0;i<par_pop;i++){
    //here we are going to have a look at the diversity in the population
    anc_array[ world[i].id ]++;
    if( world[i].health > 1.01*world[bestguy].health){
      bestguy=i;
      //printf("bestguy health = %f\n",besthealth );
    }
  }
    
  for(i=0;i<par_pop;i++) if(anc_array[i]!=0) ancestors++;
  
  if(ancestors==1){
    for(i=0;i<par_pop;i++) world[i].id=i;
  }
  
  //printf("Hello save 0.1\n");
  
  CELL *ic=&world[bestguy];
  FILE *fp;
  
  double cellvol=ic->cellvol;
  double tarvol=ic->tarvol;
  double health=ic->health;
  double R=(*(ic->Rp) + *(ic->R) + *(ic->Crmt) + *(ic->Crme) + *(ic->Crmq) + *(ic->Crmrp))/cellvol;
  double R1=(*(ic->Rr)+ *(ic->R) + *(ic->Crmt) + *(ic->Crme) + *(ic->Crmq) + *(ic->Crmrp))/cellvol;
  double Q=*(ic->Q)/cellvol;
  double T=*(ic->T)/cellvol;// + *(ic->Cst) + *(ic->Cct);
  double E=*(ic->E)/cellvol;// + *(ic->Cce) + *(ic->Cae);
  
  fp=fopen(file_out_var,"a");
  if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
  
  fprintf(fp,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d\n", Time, 
				    health, cellvol, 
				    R, R1, Q, T , E,
				    *(ic->mRp)/cellvol,*(ic->mQ)/cellvol,*(ic->mT)/cellvol,*(ic->mE)/cellvol 
				    ,*(ic->S), *(ic->C),*(ic->A)
				    , ic->growthrate, tarvol,ic->previous_lifespan, ic->id,bestguy,ancestors
	 );
  
  fclose(fp);	
  
  fp=fopen(file_out_par,"a");
  if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
  
//   fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f\n",
// 	            Time, ic->ks1,ic->ks2,ic->ks3,ic->ksm1,ic->ksm3,
// 	            ic->kc1,ic->kc2,ic->kc3,ic->kcm1,ic->kcm3);
  
  fprintf(fp,"%d %f %f %f %f %f %f",
 	            Time, ic->kt,ic->alfat,ic->betat,
		          ic->ke,ic->alfae,ic->betae);
  for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrT[i]);
  for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrE[i]);
  for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrQ[i]);
  for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrRp[i]);
  for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrRr[i]);
  
  for(i=0;i<N_GENE_TYPES;i++ ) fprintf(fp," %d",ic->genesum[i]);
  for(i=0;i<N_GENE_TYPES;i++ ) fprintf(fp," %d",ic->delet_genesum[i]);
  for(i=0;i<3;i++) fprintf(fp," %f",ic->mut_genes[i]);
  fprintf(fp," %s",ic->genome);
  fprintf(fp,"\n");
  fclose(fp);
}
*/




//for now we save only a cell
// void SaveData(double Time,CELL *world)
// {
//   
//   CELL *ic=&world[0];
//   FILE *fp;
//   
//   double R=*(ic->Rp) + *(ic->Rr)+ *(ic->R)*2. + *(ic->Crmt)*2. + *(ic->Crme)*2. + *(ic->Crmq)*2. + *(ic->Crmrp)*2.;
//   double Q=*(ic->Q);
//   double T=*(ic->T);// + *(ic->Cst) + *(ic->Cct);
//   double E=*(ic->E);// + *(ic->Cce) + *(ic->Cae);
//   double cellvol=ic->cellvol;
//   
//   fp=fopen(file_out_var,"a");
//   if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
//   
//   fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f %f %f\n", Time, cellvol, R,Q,T,E,
// 				    *(ic->mRp),*(ic->mQ),*(ic->mT),*(ic->mE), 
// 				    *(ic->S), *(ic->C),*(ic->A) );
//   
//   fclose(fp);	
//   
//   fp=fopen(file_out_par,"a");
//   if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
//   
// //   fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f\n",
// // 	            Time, ic->ks1,ic->ks2,ic->ks3,ic->ksm1,ic->ksm3,
// // 	            ic->kc1,ic->kc2,ic->kc3,ic->kcm1,ic->kcm3);
//   
//   fprintf(fp,"%f %f %f %f %f %f %f ",
//  	            Time, ic->kt,ic->alfat,ic->betat,ic->ke,ic->alfae,ic->betae);
//   
//   
//   fclose(fp);
// }
// 
// //save a cell, but al vars are relative to volume
// // and volume is scaled by target volume
// void SaveDataRelative(double Time,CELL *world)
// {
//   int i;
//   CELL *ic=&world[0];
//   FILE *fp;
//   
//   double cellvol=ic->cellvol;
//   double health=ic->health;
//   double R=(*(ic->Rp) + *(ic->Rr)+ *(ic->R)*2. + *(ic->Crmt)*2. + *(ic->Crme)*2. + *(ic->Crmq)*2. + *(ic->Crmrp)*2.)/cellvol;
//   double Q=*(ic->Q)/cellvol;
//   double T=*(ic->T)/cellvol;// + *(ic->Cst) + *(ic->Cct);
//   double E=*(ic->E)/cellvol;// + *(ic->Cce) + *(ic->Cae);
//   
//   fp=fopen(file_out_var,"a");
//   if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
//   
//   fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", Time,health, cellvol/par_init_tarvol, R,Q,T,E,
// 				    *(ic->mRp)/cellvol,*(ic->mQ)/cellvol,*(ic->mT)/cellvol,*(ic->mE)/cellvol 
// 				    ,*(ic->S), *(ic->C),*(ic->A) 
// 	 );
//   
//   fclose(fp);	
//   
//   fp=fopen(file_out_par,"a");
//   if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
//   
// //   fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f\n",
// // 	            Time, ic->ks1,ic->ks2,ic->ks3,ic->ksm1,ic->ksm3,
// // 	            ic->kc1,ic->kc2,ic->kc3,ic->kcm1,ic->kcm3);
//   
//   fprintf(fp,"%f %f %f %f %f %f %f",
//  	            Time, ic->kt,ic->alfat,ic->betat,
// 		          ic->ke,ic->alfae,ic->betae);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrT[i]);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrE[i]);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrQ[i]);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrRp[i]);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrRr[i]);
//   fprintf(fp,"\n");
//   fclose(fp);
// }
// 
// 
// 
// 
// 
// //save a cell, but al vars are relative to volume
// // and volume is scaled by target volume
// // version 2 splits Rrna contribution from Rproteins
// void SaveDataRelative2(int Time,CELL *world)
// {
//   int i;
//   CELL *ic=&world[0];
//   FILE *fp;
//   
//   double cellvol=ic->cellvol;
//   double health=ic->health;
//   double R=(*(ic->Rp) + *(ic->R) + *(ic->Crmt) + *(ic->Crme) + *(ic->Crmq) + *(ic->Crmrp))/cellvol;
//   double R1=(*(ic->Rr)+ *(ic->R) + *(ic->Crmt) + *(ic->Crme) + *(ic->Crmq) + *(ic->Crmrp))/cellvol;
//   double Q=*(ic->Q)/cellvol;
//   double T=*(ic->T)/cellvol;// + *(ic->Cst) + *(ic->Cct);
//   double E=*(ic->E)/cellvol;// + *(ic->Cce) + *(ic->Cae);
//   
//   fp=fopen(file_out_var,"a");
//   if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
//   
//   fprintf(fp,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", Time, 
// 				    health, cellvol/par_init_tarvol, 
// 				    R, R1, Q, T , E,
// 				    *(ic->mRp)/cellvol,*(ic->mQ)/cellvol,*(ic->mT)/cellvol,*(ic->mE)/cellvol 
// 				    ,*(ic->S), *(ic->C),*(ic->A) 
// 	 );
//   
//   fclose(fp);	
//   
//   fp=fopen(file_out_par,"a");
//   if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
//   
// //   fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f\n",
// // 	            Time, ic->ks1,ic->ks2,ic->ks3,ic->ksm1,ic->ksm3,
// // 	            ic->kc1,ic->kc2,ic->kc3,ic->kcm1,ic->kcm3);
//   
//   fprintf(fp,"%d %f %f %f %f %f %f",
//  	            Time, ic->kt,ic->alfat,ic->betat,
// 		          ic->ke,ic->alfae,ic->betae);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrT[i]);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrE[i]);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrQ[i]);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrRp[i]);
//   for(i=0;i<4;i++ ) fprintf(fp," %f",ic->ktrscrRr[i]);
//   fprintf(fp,"\n");
//   fclose(fp);
// }
// 
// 
// 
// void SaveData_PostReplication(double Time, CELL *world)
// {
//   CELL *ic=&world[0];
//   FILE *fp;
//   
//   double R=*(ic->Rp) + *(ic->Rr)+ *(ic->R)*2. + *(ic->Crmt)*2. + *(ic->Crme)*2. + *(ic->Crmq)*2. + *(ic->Crmrp)*2.;
//   double Q=*(ic->Q);
//   double T=*(ic->T);// + *(ic->Cst) + *(ic->Cct);
//   double E=*(ic->E);// + *(ic->Cce) + *(ic->Cae);
//   double cellvol=ic->cellvol;
//   
//   fp=fopen(file_out_post_repl,"a");
//   if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
//   
//   fprintf(fp,"%f %f %f %f %f %f\n", Time, R/cellvol, Q/cellvol, T/cellvol,E/cellvol , ic->health);
//   
//   fclose(fp);
// }
