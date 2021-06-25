#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "my.h"
#include "parameters.h"
#include "data_output.h"
#include "initialisation.h"
#include "cell.h"
#include "mersenne.h"

int InitialiseFromScratch(CELL *world,double **y)
{  
  int i,j;
  CELL *ic=&world[0];
  printf("Hello init 0.0\n");
  for(i=0;i<par_pop;i++){
    ic=&world[i];
    ic->id=i;
    
    ic->S=&y[i][_S];
    ic->A=&y[i][_A];
    
    ic->T=&y[i][_T];
    ic->Q=&y[i][_Q];
    
    ic->mT=&y[i][_mT];
    ic->mQ=&y[i][_mQ];
    ic->mRp=&y[i][_mRp];
    
    ic->Rr=&y[i][_Rr];
    ic->Rp=&y[i][_Rp];
    
//     ic->Crmt=&y[i][_Crmt];
//     ic->Crmq=&y[i][_Crmq];
    
    // --- PARAMETERS --- //
    
    // -- METABOLIC PARAMETERS
    ic->kt  	= 2.5;
    ic->alfat	= 0.5;
    ic->betat	= 1.0;

    ic->da         =  0.0050   ;
    
    ic->dt         =  0.0050     ;
    ic->dq         =  0.0050   ;
    
    ic->dmq        =  0.05000   ;
    ic->dmt        =  0.05000   ;
    ic->dmrp       =  0.05000   ;
    
    ic->drr        =  0.0050   ;
    ic->drp        =  0.0050   ;
    
    //ic->kr         =  1.00     ; // R+A+{mRNA} -> Ribosomal Complex
    //ic->kcr        =  2.500     ; // Ribosomal Complex -> R + {T,E,Q,Rp} ~~mRNA quickly degraded after translation
    ic->krprr      = 0.9  ;	// fraction of Rp and Rr in complex 
    ic->sM         = 10.   ;	//scaling of protein production: mRNA*Rr*A/( EPSILON + sM*mRNA + sR*Rr + sA*A )
    ic->sR         = 1.   ;
    ic->sA         = 3.   ;	//to make ribosome most important <- why do I think this is so important?
    
    // REGULATION - initialised with evolved parameters -- or something, given how much they fluctuate
    //0:constitutive, 	 	1: A-dep, 	
    ic->ktrscrT[0]  = 0.1; ic->ktrscrT[1]  = 0.01;
    ic->ktrscrQ[0]  = 0.1; ic->ktrscrQ[1]  = 0.01;
    ic->ktrscrRr[0] = 0.1; ic->ktrscrRr[1] = 0.01;
    ic->ktrscrRp[0] = 0.1; ic->ktrscrRp[1] = 0.01;
       
    // MUTATION
    //ic->mut_genes[0]=1./3.;	// duplications
    //ic->mut_genes[1]=1./3.;	// deletions
    //ic->mut_genes[2]=1.-ic->mut_genes[0]-ic->mut_genes[1];	// deleterious mutations
    
    // MUTATION - max rate of transcriptional mutations <- evolvable
    ic->max_ktrmut[0]=par_max_ktrmut/3.;	// duplications
    ic->max_ktrmut[1]=par_max_ktrmut/3.;	// deletions
    ic->max_ktrmut[2]=par_max_ktrmut/3.;	// deleterious mutations
    
    
    //VARIABLES --guessed
    *(ic->S)          =  10.00;
    *(ic->A )         =  0.100;
    
    *(ic->Rr )        =  1.9;
    *(ic->Rp )        =  1.0;
    *(ic->Q )         =  2.10;
    
//     *(ic->Crmt )      =  0.5;
//     *(ic->Crmq   )    =  0.5;
    *(ic->T  )        =  2.000;
    
    *(ic->mT  )       =  0.5;
    *(ic->mQ  )       =  0.5;
    *(ic->mRp  )      =  0.5;
    
    ic->growthrate=0.;	//this is only for statistics
    ic->lifespan=0.;
    ic->av_tr_mut=0.; 
    
    ic->celldiv=0;
    for(j=0;j<N_GENE_TYPES;j++) ic->trload[j]=0.;
    
    // Initialise genome
    for(j=0;j<MAXGENOMESIZE;j++) ic->genome[j]='\0';
    
    //strcpy(ic->genome,init_genome); //
    //random initialisation
    
    //ic->genomesize=10+(int)(10.*genrand_real2()); //random initial genome size between 10 and 20
    
    ic->genomesize=5;//10+(int)(10.*genrand_real2()); //random initial genome size between 10 and 20
    char possible_genes[]={'T','Q','R','P'}; //NOPE ->make 'B' a little more likely at the beginning
    int len_possible_genes=4;
    
    for(j=0;j<ic->genomesize;j++){
      int which = len_possible_genes*genrand_real2();
      ic->genome[j]=possible_genes[which];
    }
    
    
    for(j=0;j<MAXGENOMESIZE;j++){
      //fprintf(stderr,"hello0.1\n");
      if(ic->genome[j]=='\0') break;
    }
    ic->genomesize=j;
    
    CellVolume(ic);	// this should always come before cellhealth
    GenomeStats(ic);	// this takes genome and calculates genesum, delet_genesum and genomesize
    
    CellHealth(ic);
    //printf("Hello init 0.4\n");
    TargetVolume(ic);
    //printf("Volume: %f, Target Volume %f, Health %f\n",ic->cellvol,ic->tarvol,ic->health);
    //InitialiseAncestry(ic,i);
    
    ic->anc=i;	//ancestry
  }
  
//   for(i=0;i<par_pop;i++){
//     ic=&world[i];
//     fprintf(stderr,"Initialisation(): all done, genome=%s, genomesize=%d\n",ic->genome, ic->genomesize);
//     //exit(1);
//   }
  
  if(MUT_REGULATION==0) fprintf(stderr,"Initial(), Warning: No mutations of regulation parameters.\n");
  if(MUT_TRANSCR_LOAD==0) fprintf(stderr,"Initial(), Warning: No mutations from transcriptional load\n");
  if(MUT_BACKGROUND==0) fprintf(stderr,"Initial(), Warning: No backgound mutations\n");
  if(MUT_MUTATION==0) fprintf(stderr,"Initial(), Warning: No mutations of mutation parameters\n");
  
  for(i=0;i<N_MUT_TYPES;i++) if(BGmut_scheme[i] != TRmut_scheme[i]){
    fprintf(stderr,"Initial(), Warning: BGmut_scheme != TRmut_scheme.\n");
    break;
  }
  
  fprintf(stderr,"Initialisation(): complete\n");
  return 0;
}

// If output files already exist...
// you are probably starting the simulation by mistake
void CheckOutputFilesExist(void)
{
  if( access( file_out_var, F_OK ) != -1 ){
    // file exists already, which likely means that the simulation is already going
    // and we should not start another one
    fprintf(stderr,"CheckOutputFilesExist(): Error. Output files already exist\n");
    fprintf(stderr,"Delete or move them before starting a new one\n");
    exit(1);
  }
  
  // else go on with simulation
  //else{
  //  // file doesn't exist
  //}
}

//checks that environmental switches are properly set
void CheckEnvironmentSwitches(void)
{
  int sum_env=ENV_S_RAND+ENV_S_INCR+ENV_S_CONST;
  
  if(sum_env!=1){
    fprintf(stderr,"Initial(): Error. Environment is not set properly\n");
    fprintf(stderr,"Set only ONE of the three possible variables to 1, and the other two to zero.\n");
    fprintf(stderr,"The program will exit now. Goodbye.\n");
    exit(1);
  }

}






// int InitialiseFromBackup(CELL *world, double **y)
// {
//   int i,fi,pi,fpar_pop;
//   char instr[MAXLEN];
//   CELL *ic;
//   FILE *fin;
//   
//   //to initialise variables, and to assign pointers, we will overwrite parameters
//   InitialiseFromScratch(world,y);	
//   
//   //parameters initialisation
//   fin=fopen(file_input,"r");
//   if(fin==NULL) {
//     fprintf(stderr, "InitialiseFromBackup(): Error. Can't open %s\n",file_input);
//     exit(1);
//   }
//   
//   //read heading of backup file: "#Time nr X nc"
//   fgets(instr,MAXLEN,fin);
//   strtok(instr," ");	//strip 'Time'
//   strtok(NULL," ");	//strip actual Time
//   strtok(NULL," ");	//strip actual par_pop
//   
//   fpar_pop = atoi(strtok(NULL," ")); //read par_pop_from file
//   
//   //you must loop around the smaller between par_pop and fpar_pop
//   //but population initialisation can be only until par_pop
//   for(fi=0,pi=0; pi<par_pop; fi++,pi++) {
//     
//     if(fi==fpar_pop && fi<par_pop) {
//       //we put more individuals from the same file...
//       //and we copy elements, and exit afterwards
//       for(i=0; i<par_pop-fpar_pop; i++)	world[i+fpar_pop]=world[i];
//       break;
//     }
//     
//     ic=&world[pi];
//     //read a line
//     if(fgets(instr,MAXLEN,fin)==NULL){
//       fprintf(stderr, "InitialiseFromBackup(): Error at reading file\n");
// 	exit(1);
//     }
//     
//     ic->kt=atof(strtok(instr," "));	//kt
//     ic->alfat=atof(strtok(NULL," "));	//alfat
//     ic->betat=atof(strtok(NULL," "));	
//     ic->ke=atof(strtok(NULL," "));	
//     ic->alfae=atof(strtok(NULL," "));	
//     ic->betae=atof(strtok(NULL," "));	
//     
//     for(i=0;i<4;i++ ) ic->ktrscrT[i]=atof(strtok(NULL," "));
//     for(i=0;i<4;i++ ) ic->ktrscrE[i]=atof(strtok(NULL," "));
//     for(i=0;i<4;i++ ) ic->ktrscrQ[i]=atof(strtok(NULL," "));
//     for(i=0;i<4;i++ ) ic->ktrscrRp[i]=atof(strtok(NULL," "));
//     for(i=0;i<4;i++ ) ic->ktrscrRr[i]=atof(strtok(NULL," "));
//     
//     CellVolume(ic);
//     CellHealth(ic);
//     ic->growthrate=0.;
//   }
//  
//  fprintf(stderr,"InitialiseFromBackup(): Error, no genetics intialised");
//  //add the genes initialisation -----------------------------
//  exit(1);
//  
//  for(i=0;i<par_pop;i++) InitialiseAncestry(ic,i);	//put here because it's safer
//  return 0; 
// }


// void InitialiseAncestry(CELL *ic, int i)
// {
//   sprintf(ic->anc[0],"-1");
//   sprintf(ic->anc[1],"0-%d",i);  
// }










