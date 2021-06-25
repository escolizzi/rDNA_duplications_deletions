/*
  Evolution of transcr. duplications deletions and inact independent from each other.
  deltas are small and rate can go negative (but treated as zero)
*/

/*
 * Corrected minor (?) bug in translation of Rp
 * Life span is printed down the line of descent now!
 */

/*
 * 1) Change regulation so that evol. parameters can evolve between some large 
 *     Max and Min (of course it is still the case that the actual load is bounded 
 *     between 0 and 0.2) -done
 * 2) Print some info on tr_load and mut rate (just average tr. mut. for now)
 * 3) Add option to average environment every integration step -NOT YET
 * 4) Make switching time a parameter in parameters.c -done
 * 5) SingleGeneMut() takes mutational scheme as extra argument
 * 6) bug fixed: regulation did not work because func does not do what I think it does.
 * 		 Regulation now happens outside of func, i.e. in the main loop
 */

/*
 * 1) Background mutational scheme is in parameters
 * 2) 1 step of ancestor tracing is done while the program runs.
 *    The entire population is stil dumped periodically, but less frequently
 */

/*
 * This is a branch, 
 * essentially, if break points are not activated, the result is just an increase of background mutations,
 * which are chosen depending on the scheme (dupl,del,inact)=(1/4,1/4,1/2) or (1/3,1/3,1/3) depending on the simulation.
 * 
 * Technically, we make a more "dynamic" use of the variable "singlegene_mut_flag" (see TranscrMut() in cell.c)
 */

/*
 * In this version we put back in Rp, because yes, ok? 
 * Anyway, this goes that Rr and Rp form a complex, 
 * which is assumed in equilibrium with their relative amounts (sic!)
 * Rr+Rp<=>R; R=K*min(Rr,Rp).. so that you keep stoichiometry
 */

/*
 * This is a tentative branch: 
 * - break points are inserted randomly in the genome proportional to their fraction (see backgroundmutations() in cell.c)
 * - we remove the complex formation variables for the ribosome and introduce a phenomenological term:
 *   protein produciton per time step = mRNA*Rr*A/( EPSILON + sM*mRNA + sR*Rr + sA*A )
 *   with km=10. and kA=kRr=1
 */

/*
 * Here we re-complexify the model a bit (1 parameter, set to 1 :P ) 
 * by introducing saturation of A-usage by ribosome.. this is a test
 */

/*
 * IN THIS VERSION THERE ARE ONLY THREE TYPES gene types: T (general enzime), Q and Rr
 * and 2 molecules, S and A
 * all is the same as before, except that it should be somewhat simplified
 * 
 * - ancestor trace is implemented as a two-tags system (see my.h and data_output.c)
 * - some heavy bugs corrected in background mutations (in cell.c)
 */

/*
 * Split y into many y's
 */

/*
 * In this version cells are given some time to reach target volume,
 * which in turn determines fitness
 * taking too long results in death.
 * Individuals are ranked and so many are replicated in the next round.
 * 
*/
/*
  V.something: 
  linear genome, beads on a string:
  different genes are concatenated downstream of breaking points,
  which encode mutational information, and evolve globally (i.e. there is one per cell)
  -no more mutations for metabolic genes
*/

/*
  par_max_expression not used anymore
  V0.4:
  transcription affects both active and inactive genes, so do mutations,
  the amount of mRNA used by the cell, however, is only the proportion of active ones
*/

/*
  V0.3: genome dynamics, 
  perhaps fix metabolic parameters
  let transcription evolve
  let number of genes evolve
*/

/*
  V0.2: modify S through time 
*/

/*
  V0.1
*/

/*A simple model of metabolism and genome dynamics*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "mersenne.h"
#include "parameters.h"
#include "cell.h"
#include "data_output.h"
#include "initialisation.h"
#include "my.h"

int Time;
double timestep;

//________________________________//
//                                //
//                                //
// ----- THE SYSTEM OF ODEs ----- //
//                                //
//________________________________//

int func(double t, const double y[], double f[], void *params)
{
  
  int j;
  
//   printf("func Hello 0.0\n");
  
    CELL *ic = ( CELL *) params;
//     printf("func Hello 0.1\n");
    if(ic->celldiv==1){
      for(j=0;j<N_VARIABLES;j++) f[j]=0.; //if we are already at target volume, do nothing
      return GSL_SUCCESS;
    }
//     printf("func Hello 0.2\n");

    //Environment
  double S   = y[_S];
  double ds  = par_ds;	//external parameters
  double Sin = par_Sin;
  
  //metabolites
  double A = y[_A];
  
  //metabolism 
  double T = y[_T];
  //housekeeping
  double Q = y[_Q];
  //mRNAs
  double mT = y[_mT];
  double mQ = y[_mQ];
  double mRp = y[_mRp];
  //ribosome and complexes
  double Rr = y[_Rr];
  double Rp = y[_Rp];
//   double Crmt = y[_Crmt];
//   double Crmq = y[_Crmq];

  // parameters
  
  double da = ic->da;  
  double dmq  = ic->dmq;
  double dmt   = ic->dmt;
  double dmrp   = ic->dmrp;
  double dq   = ic->dq;
  double drr   = ic->drr;
  double dt   = ic->dt;
  double drp   = ic->drp;
 
  double kt = ic->kt;
  double alfat= ic->alfat;
  double betat= ic->betat;
 
  //double kcr   = ic->kcr;
  //double kr     = ic->kr;
  double krprr=ic->krprr;
  double sM=ic->sM;
  double sR=ic->sR;
  double sA=ic->sA;
  
  //double *ktrscrT  = ic->ktrscrT;	//taken care of this outside
  //double *ktrscrQ = ic->ktrscrQ;
  //double *ktrscrRr = ic->ktrscrRr;
  //double *ktrscrRp = ic->ktrscrRp;
    
  double R = krprr*fmin(Rp,Rr);	// some equilibrium value and the fact that they are stoichiometrically hard-coupled
  
  // -- Transcription -- //
  double tr_Rr  = ic->trnscr_per_timestep[GENE_Rr];
  double tr_mRp = ic->trnscr_per_timestep[GENE_Rp];
  double tr_mT = ic->trnscr_per_timestep[GENE_T];
  double tr_mQ = ic->trnscr_per_timestep[GENE_Q];
  
  //Up to here only assignments, now cell dynamics: 
  
  
  // -- RNA dynamics (including Rr) -- //
  //          transcr - ribosome association - ribosome diss. - decay  - Rp in ribosome decays
  f[_Rr]= tr_Rr - drr*Rr ;
  //     transcr - translation - decay - ribosome decay
  f[_mRp]= tr_mRp - dmrp*mRp ;
  f[_mT] = tr_mT  - dmt*mT ;
  f[_mQ] = tr_mQ  - dmq*mQ ;
  
  // --Translation -- //
  double Tproduction = mT*R*A/(EPSILON + sM*mT + sR*R + sA*A);	// protein production based on Ribosome, aminoacids and mRNA
  double Qproduction = mQ*R*A/(EPSILON + sM*mQ + sR*R + sA*A);
  double Rpproduction = mRp*R*A/(EPSILON + sM*mRp + sR*R + sA*A);
  
  // -- Metabolism -- //
  f[_S] = Sin - ds*S - kt*T*alfat*S/(alfat*S + betat*A + 1.) ;
  f[_A] = -da*A + kt*T*alfat*S/(alfat*S + betat*A + 1.) - Tproduction - Qproduction - Rpproduction ;
  
  // -- Protein lives -- //
  f[_T] = Tproduction - dt*T;
  f[_Q] = Qproduction - dq*Q;
  f[_Rp] = Rpproduction - drp*Rp;
  
  
  
  
//   //          transcr - ribosome association - ribosome diss. - decay  - Rp in ribosome decays
//   f[_Rr]= tr_Rr - drr*Rr - kr*(A/(1.+A))*Rr*(mT+mQ) + kcr*(Crmt + Crmq);
//   //     transcr - translation - decay - ribosome decay
//   f[_mT] = tr_mT  - kr*(A/(1.+A))*Rr*mT - dmt*mT + drr*Crmt;
//   f[_mQ] = tr_mQ  - kr*(A/(1.+A))*Rr*mQ - dmq*mQ + drr*Crmq;
//   
//   //     Rib complex form  - new proteins  - rib decay
//   f[_Crmt] = kr*(A/(1.+A))*Rr*mT - kcr*Crmt - drr*Crmt;	//assume that only ribosome dacays while in complex (mrna is stabilised?)
//   f[_Crmq] = kr*(A/(1.+A))*Rr*mQ - kcr*Crmq - drr*Crmq;
//   
  return GSL_SUCCESS;
}

//returns transcription
//transcription regulation is very simple here: 
// - If your Q volume is far from target, you transcribe less
// - active and inactive genes are both transcribed and contribute to tr_load, 
//   but only the active ones actually contribute to translation
double TranscrRegulation(CELL *ic,int genetype, double *ktrscr,double A, double integr_step_size)
{
  int activegenes=ic->genesum[genetype];
  double act_plus_inact_genes = activegenes + ic->delet_genesum[genetype];
  
  if(activegenes==0) return 0.;		//no active genes, no transcr
  
  double amount;
  
  amount = ic->health * (  ktrscr[0] + ktrscr[1]*A);
  
  if (amount<=0.) return 0.;
  
  // There is a cap to how much gene expression you can have per gene per time step
  // because you cannot fit infinite RNApolymerases on a gene
  if(amount > par_max_expression_per_gene*act_plus_inact_genes) amount= par_max_expression_per_gene * act_plus_inact_genes;
  
  amount *=integr_step_size; //multiply by integration step size to scale rates.
  ic->trload[genetype] += amount;
  
  double fract_act_actplusinact = activegenes / act_plus_inact_genes;
  return fract_act_actplusinact*amount;

  //return amount;
}

int main(int argc, char **argv)
{
  CELL *world;	//contains parameters and pointers to variables
  CELL *ic;	// pointer to one guy in world
  CELL *prev_world; //this is for ancestor trace
  double **y;	//contains all variables
  
  int i,j;
  double max_growthtime= 100. ;
  double delta_time;
  int slowing_down=0;
  fprintf(stderr,"Hello main 0.0\n");
  
  world=AllocateWorld(par_pop); //Allocate a population of CELLs
  prev_world=AllocateWorld(par_pop); //for ancestor trace
  
  y=Allocatey(par_pop);	// ...and double array of variables y and returns y
  
  for(i=0;i<par_pop;i++)for(j=0;j<N_VARIABLES;j++) y[i][j]=i*N_VARIABLES+j;
  
  Initial(world,y);		// Initialise it
  //printf("Hello main 0.02\n");
  
  for(i=0;i<par_pop;i++) prev_world[i]=world[i]; // initialise ancestry
  
  // This initialises the ode system: its 4 arguments are
  // - func is the function (pointer) where the actual system is specified
  // - jac is the Jacobian of the system, which we do not have and we do not need, so we put NULL
  // - dimension is the total number of equations we deal with= population size * number of variables
  // - params is a pointer to void and should contain parameters, but since it can incorporate anything
  //   we pass a pointer of the whole system
  
  //gsl_odeiv2_system sys = { func, NULL, N_VARIABLES, world };
  
  gsl_odeiv2_system *sys= malloc( par_pop*sizeof( gsl_odeiv2_system ));
  for(i=0;i<par_pop;i++)
    sys[i]= (gsl_odeiv2_system){ func, NULL, N_VARIABLES, &world[i] }; // it works if you cast { func, NULL, N_VARIABLES, world } to gsl_odeiv2_system type... for some reasons.. in C99 and later
  
  
  // This initialises more things, and defines the pointer d which we should free at the end of the sim.
  // Essentially its arguments are the ode system, the integrator, and the accuracy it should keep
  // DEFAULT: gsl_odeiv2_step_rk8pd
  //
  // Other available explict algorithms, benchmarked over 2000 timesteps,
  //Step Type: gsl_odeiv2_step_rk2
  // Explicit embedded Runge-Kutta (2, 3) method. 
  //Step Type: gsl_odeiv2_step_rk4
  // Explicit 4th order (classical) Runge-Kutta. Error estimation is carried out by the step doubling method. For more efficient estimate of the error, use the embedded methods described below. 
  //Step Type: gsl_odeiv2_step_rkf45 - 3m4.893s (second time 3m3.206s)
  //  Explicit embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a good general-purpose integrator. 
  //Step Type: gsl_odeiv2_step_rkck - 4m19.473s
  //  Explicit embedded Runge-Kutta Cash-Karp (4, 5) method. 
  //Step Type: gsl_odeiv2_step_rk8pd - 4m41.648s
  //  Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method. 
  
  //what was there before
  //gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);  
  
  gsl_odeiv2_driver **d;
  d = malloc( par_pop*sizeof( gsl_odeiv2_driver* ));
  for(i=0;i<par_pop;i++) d[i]= gsl_odeiv2_driver_alloc_y_new(&sys[i],gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
  
  // and then you call odes_status = gsl_odeiv2_driver_apply(d[i], &t, timestep, y[i]);	//this intergrates
  
  fprintf(stderr,"Warning: no mutations in the metabolic parameters\n");
  fprintf(stderr,"Simulation initialised (very few checks implemented).\n\t\t Good luck\n");
  // constant time steps
  delta_time = 1./par_init_interval_time;
  
  /******************/
  /*** MAIN CYCLE ***/
  /******************/
  int odes_status,cells_status;
  //printf("Hello main 0\n");
  
  double some_time=0.;
  
  for(Time = par_initTime; Time <= par_maxTime; Time++){
    
    // SAVING
    if( Time % par_time_savedata == 0 ){
      //printf("Hello main, Save data, Time = %d\n",Time);
      SaveDataRelative4(Time, world);
    }
    
    //we do time !=0 because we actually print the past
    if( Time % par_time_anc==0 && Time!=0){
      OutputBackup(Time, world,prev_world); // This also deals with ancestry, partially
      
      for(i=0;i<par_pop;i++) prev_world[i]=world[i]; // update prev_pop, for when it will be used in the future
      
      for(i=0;i<par_pop;i++) world[i].anc=i;	//resets ancestor after you updated prev_pop
      
      //fprintf(stderr,"Hello main, Backup data, Time=%d\n",Time);
    }
    
    //if( Time % par_time_anc==0) CheckPrintAssignAncestry(Time,world);
    
    if(Time % par_time_change_Sin == 0){
      // if(ENV_S_CONST) do nothing, S_in keeps its original value par_Sin
      if(ENV_S_RAND){
        if(par_Sin==1.) par_Sin=10.;
        else par_Sin=1.;
      }else if(ENV_S_INCR){
        par_Sin=10.+0.0005*Time;
      }
      par_time_change_Sin=par_period_time_change + (int)(0.1*par_period_time_change*(genrand_real2()-0.5));
    }
    
    
    //
    //  --  SUBCYCLE  --  one generation
    //
    
    int tot_repl=0;
    
    for(i=0;i<par_pop;i++){
      double t=0.;
      cells_status=0;
      for(timestep=0.; timestep < max_growthtime; timestep+=delta_time){
        
//         if(0==i){
// 	  SaveDataRelative4(some_time, world);
// 	  some_time+=delta_time;
// 	}
        
        // REGULATION //
        ic=&world[i];
        ic->trnscr_per_timestep[GENE_Rr]  = TranscrRegulation(ic,GENE_Rr,ic->ktrscrRr,y[i][_A],delta_time);	//trload is updated inside the function
        ic->trnscr_per_timestep[GENE_Rp]  = TranscrRegulation(ic,GENE_Rp,ic->ktrscrRp,y[i][_A],delta_time);
	ic->trnscr_per_timestep[GENE_T]   = TranscrRegulation(ic,GENE_T, ic->ktrscrT, y[i][_A],delta_time);
	ic->trnscr_per_timestep[GENE_Q]   = TranscrRegulation(ic,GENE_Q, ic->ktrscrQ, y[i][_A],delta_time);
	        
        // - NOT YET
        //The environment may be averaged every integration step, 
        // all individuals receive the same amount of resources 
        // no matter how much they consumed
        //if(AVRG_ENV){
	//  
	//}
        
        // INTEGRATION //
        odes_status = gsl_odeiv2_driver_apply(d[i], &t, timestep, y[i]);	//this intergrates
        
	if(odes_status != GSL_SUCCESS){
          printf("error, return value=%d\n", odes_status);
          break;
        }
        
        cells_status=UpdateCell(&world[i],timestep);// Update Cell.. also volume!
        
        if(cells_status==1){
	  
	  //world[i].x=world[i].trload[GENE_Rr];  // DELETE ME WHEN YOU ARE DONE TESTING
	  
	  tot_repl++;
	  break; //there is no need to continue the loop, the cell made it
	}
	//else world[i].x=-1.;			// DELETE ME WHEN YOU ARE DONE TESTING
	
      }
    }
    
    // Now tot_repl is larger than zero if at least one cell reached tarvol
    if(tot_repl>0){
      
    if(slowing_down!=1){  
      // TEST AND CHECKS //
    // This will likely fail every time replication slows down, so we put it 2 time steps further... don't know if it works
    for(i=0;i<par_pop;i++){
      ic=&world[i];
      //double av_trload_pert=ic->x/(ic->lifespan * ( ic->genesum[GENE_Rr] + ic->delet_genesum[GENE_Rr] ));
      double max_trload= par_max_expression_per_gene* (delta_time + ic->lifespan)*( ic->genesum[GENE_Rr] + ic->delet_genesum[GENE_Rr] );
      //fprintf(stderr,"Time=%d, trload=%f, max_possible_trload=%f, lifespan=%f, active Rr %d, inact Rr %d\n",Time, ic->x, max_trload , ic->lifespan,ic->genesum[GENE_Rr],ic->delet_genesum[GENE_Rr]);
      if(ic->x > EPSILON + max_trload){
	fprintf(stderr,"\nTime=%d, i=%d, trload=%f, max_possible_trload=%f, lifespan=%f, active Rr %d, inact Rr %d\n",Time,i, ic->x, max_trload , ic->lifespan,ic->genesum[GENE_Rr],ic->delet_genesum[GENE_Rr]);
	fprintf(stderr,"JEEEEEEEE\n");
	exit(1);
      }
    }
    }
      
      
      
      //Replication returns max_growthtime based on 
      // how long it took for the slowest guy that reached tarvol to actually reach it.
      max_growthtime=2.*Replication(world,y,cells_status,max_growthtime);
      
      // TEST AND CHECKS //
      if(max_growthtime<0.){
	fprintf(stderr,"How is this possible? max_growthtime less than zero=%f\n", max_growthtime);
	exit(1);
      }
      
      if(max_growthtime>MAX_GROWTH_TIME) max_growthtime=MAX_GROWTH_TIME;
      
      cells_status=0;
      if(slowing_down==1) {
	fprintf(stderr,"Hello main 0.3, succesful replications after longer time span, setting 'slowing_down' to 0\n");
	slowing_down=0;
      }
        
    }else if(slowing_down==0){
      fprintf(stderr,"main(), Generation: %d. Warning: no cell replicated within the allowed time = %f\n",Time,max_growthtime);
      fprintf(stderr,"... trying to give it some more time\n");
      max_growthtime=MAX_GROWTH_TIME;
      slowing_down=1;
      
    }else{
      fprintf(stderr,"main(), Warning: no cell replicated within the allowed time = %f\n",MAX_GROWTH_TIME);
      fprintf(stderr,"The simulation will finish now\n");
      break;
    }
  }
  
  Finish(world,y);
  return 0;
}

//this does cell stuff:
// - re-calculates volumes 
// - determines chances of mutations (?)	- NO
// - and overall fitness?			- YES
// returns zero if nothing is there to do for fitness (no cell needs to divide) -YES
int UpdateCell(CELL* ic,double time2tarvol)
{
  int i,status=0;
  
  if(status){
    fprintf(stderr,"UpdateCell(), Error. Got a cell with status 1\n");
    exit(1);
  }
  
  ic->growthrate=ic->cellvol;
  CellVolume(ic);
  
  ic->growthrate= (ic->cellvol - ic->growthrate)*par_init_interval_time;// actually this is the metabolic flux: Delta Vol/timestep = Delta Vol*Delta Vol
  
  CellHealth(ic);
  
  if(ic->celldiv==0 && ic->cellvol >= ic->tarvol){
    ic->celldiv=1;	// if cell is large enough it can divide
    ic->lifespan=time2tarvol;
    status=1;
  }
  
  return status;
}

void Change_Sin()
{
  par_Sin = 1000*pow(genrand_real1(),10) ;
  if(par_Sin>1000.) par_Sin=2000.-par_Sin;
  if(par_Sin<1.) par_Sin=2.-par_Sin;
}

//variables inside cells are just pointers to the corresponding value in y
// therefore, world[i].S points to y[i*N_VARIABLES + _S]
void Initial(CELL *world,double **y)//, gsl_odeiv2_driver *d)
{
  int i,j;
  //CELL *ic;
  
  init_genrand(par_ulseed);	//initialise random number generator with seed par_ulseed (see parameters.c)
  
  CheckOutputFilesExist();
  CheckEnvironmentSwitches();
  
  if( strlen(file_input)==0 )
    InitialiseFromScratch(world,y);
  else{
    fprintf(stderr,"Initial(), Error. function to initialise from backup is not really working\n");
    exit(1);
    InitialiseFromBackup(world, y);
  }
  
}

// allocate world and y and returns y (to avoid the mess of triple pointers... holy cow I hate triple pointers)
CELL *AllocateWorld(int par_pop)
{
  CELL *a;
  a=malloc(par_pop*sizeof(CELL));
  if(a==NULL){
    fprintf(stderr,"Initial(): Error in memory allocation.\n");
    exit(1);
  }
  return a;
}

double **Allocatey(int par_pop){
  int i;
  double **b;
  
  b= malloc(par_pop*sizeof(double*));
  for(i=0;i<par_pop;i++) b[i]=malloc(N_VARIABLES*sizeof(double));
  
  if(b==NULL){
    fprintf(stderr,"Initial(): Error in memory allocation.\n");
    exit(1);
  }
  
  return b;
}


void CheckPrintAssignAncestry(int Time,CELL *world)
{
  
  ;
  
}

void Finish(CELL *world,double **y)
{
  int i;
  free(world);
  for(i=0;i<par_pop;i++) free(y[i]);
  free(y);
  fprintf(stderr,"\n\n One day please use gsl_odeiv2_driver_free(d) and free stuff properly\n");
  //gsl_odeiv2_driver_free(d);
}



