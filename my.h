#ifndef MY_H
#define MY_H

#include "constants.h"

// ---------------------------------------------------
//  A CELL, rememember to declare CELL variable below 
// ---------------------------------------------------
typedef struct __cell_FOR_LATER{ 
  int id;		//position for retrieval of data from variables' vector y
  
  double cellvol; 	// the sum of all proteins and Rr
  int celldiv;		//a flag for dividing
  double tarvol;	//target volume, depends on genome size and maybe one day on metabolites
  double minvol;	//minvol (depends on genome size), below this cell suffers
  double health;	//how is the cell doing based on all contraints
  double growthrate;	//Delta volume, between timesteps
  double lifespan; //for some statistics, like growth rate, not used for anything
  double av_tr_mut;// also for statistics
  double x; // another variable for printing stuff
  
  int cellmaxgenomesize; //NOT YET
  
  // proteins in various states -- this is summed for cell volume
  double *T;//,*Cst,*Cct;
  
  double *Q;
  double *Rr,*Rp;
  //double *Crmt,*Crmq;
  
  //mRNAs, except ribosomal RNA - Rr, which is above
  double *mT,*mQ,*mRp;
  
  // metabolites
  double *S; //this is not outside
  double *A;
  
  // parameters for all reactions
  // tentatively divided by those to optimise first (and then leave there)
  // and those to evolve later  
  
  // parameters to optimise 
  //double Sin, ds; //this will be outside
  double kt,alfat,betat; 
  //double ke,alfae,betae;
  //double ks1,ks2,ks3,ksm1,ksm3; //transport rates: S+T ->...->C+T
  double dt;	// decay of transporter
  //double dc; 	// decay C
  //double kc1,kc2,kc3,kcm1,kcm3; //metabolic rates: C+E ->...->A+E
  //double de;	// decay of metabolic enzime
  double da; 	// decay of A
  double dq;	//decay of homeostatic proteins
  double drp;  //decay Ribosomal proteins
  
  //double dme,dmq,dmrp,dmt,drp,drr; //decay mRNAs
  double dmq,dmt,dmrp,drr; //decay mRNAs
  //double kcr;	// R+A+{mRNA} -> Ribosomal Complex
  //double kr; // Ribosomal Complex -> R + {T,Q} ~~mRNA quickly degraded after translation
  double krprr; //fraction of Rr and Rp in complex     // OLD STUFF,kmrprr;	// Rp + Rr <=> R
  double sM,sR,sA; //scaling of protein production: mRNA*Rr*A/( EPSILON + sM*mRNA + sR*Rr + sA*A )
  
  // transcription rates from A-activated promoter and C-activated promoter
  // as well as both and constitutive
  // actual transcription is ktrscr[0] + ktrscr[1]*C + ktrscr[2]*A + ktrscr[3]*A*C;
  // it is important that at some point (this can be repressive as well)
  // i.e. all of them can be positive or negative or zero
  // also, I'm not bounding transcription, but mutation rates must grow exponentially with it
  
  double ktrscrT[2];//0: constitutive, 1: A-dependent
  //double ktrscrE[4]; 
  double ktrscrQ[2];
  //double ktrscrRp[4];
  double ktrscrRr[2];
  double ktrscrRp[2];
  
  //ancestor's tag, the position in i 
  // given at the time of previous backup:
  // at next backup I write to file this tag and the current position, 
  // then give the current position as a tag and so on
  int anc;
  
  //-- EXPLICIT genome --//
  // the genome is just a long string of characters:
  // T,E,Q,P,R,B,E
  // P is Rp, R is Rr, B is Breaking point
  // E is end of genome
  // when letters are capitalised that's active genes
  // when they are small, they are inactive
  char genome[MAXGENOMESIZE];	//the genome
  //char *genome;	//the genome
  int genomesize; //tot elements
  
  // Implicit genome, there 5 kinds of genes: 
  //  geneT,geneE,geneRp,geneRr,geneQ;	
  // The genome says how many copies per gene,
  //  from which some formula should yield 
  //  the actual transcription load based mutation rate
  
  int genesum[N_GENE_TYPES]; //the sum of active genes per type
  int delet_genesum[N_GENE_TYPES]; //the sum of inactive genes per type
  
  double trnscr_per_timestep[N_GENE_TYPES]; //mrna produced every integration step
  
  double trload[N_GENE_TYPES]; //total transcription load
  double mut_genes[3]; //probability that a mutation is insertion,deletion,deleterious. Sums to 1
  double max_ktrmut[3]; // max rates of transcr. induced dupl, dels, inact. Does NOT sum up to 1,
                        // it replaces mut_genes[3] which should not be used in this model.
  
} CELL;


//---- FUNCTIONS
CELL InitOneCell();
void Change_Sin();
double TranscrRegulation(CELL *ic,int genetype,double *ktrscr,double A, double integr_step_size);
CELL *AllocateWorld(int par_pop);
double **Allocatey(int par_pop);
void Initial(CELL *world, double **y);//, gsl_odeiv2_driver *d);
void CheckPrintAssignAncestry(int Time,CELL *world);
int UpdateCell(CELL* ic,double time2tarvol);
void Finish(CELL *world,double **y);

#endif


