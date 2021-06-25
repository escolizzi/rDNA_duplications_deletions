#include "parameters.h"
#include "constants.h"
#include "limits.h"

int par_initTime=0;
int par_maxTime=INT_MAX;
double par_init_interval_time=20.;
long int par_ulseed=267247;

int par_time_savedata=10; //make sure this is integer
int par_time_backup=5000;//makes a dump of the whole population
int par_time_anc=10;	//this can mean two things, either how frequently we save tag for ancestry 
			// or how frequently we check everyone's ancestor
double par_period_time_change = 200.; //period after which env. changes (up to stochasticity)

// MUTATIONS - see cell.c, function MutateCell()
// 1- the mutation is on, 0- mutation does not happen
const int MUT_REGULATION= 1;
const int MUT_TRANSCR_LOAD= 1;
const int MUT_BACKGROUND= 1;
const int MUT_MUTATION= 0;

const int ENV_S_RAND =0;	//Environment: random switching between high and low (see main.c)
const int ENV_S_INCR =0;	//Environment: S_in increases steadily over time
const int ENV_S_CONST=1;	//Environment: constant S_in, which never changes

const int AVRG_ENV=0; //NOT IMPLEMENTED YET -average environment every time step?

const int par_pop=1000;
char init_genome[MAXGENOMESIZE]="BTQRPB";
double par_init_tarvol=30.; //vol that triggers division: 
			    // this should scale with growth rate 
			    // and should be larger than min_vol
double par_init_minvol=25.; //below this cells suffer, 
			    // this should scale with genome size
double par_tar_Q=0.3;      // target proportion of Q relative to volume
double transcr_penalty=15.; // transcription slows due to distance from par_tar_Q
			    // by exp( -transcr_penalty*( Q/vol - par_tar_Q )^2 )
			    // notice that the beta for which a difference of 0.05 
			    // leads to a fitness difference of 0.05 is ~20.5
double par_ds=0.1;
double par_Sin=10.;
int par_time_change_Sin=50.;
double par_DNAbeta=1.;  //foldedness of DNA
double par_ktrmut=0.90;	//constant for transcriptional mutations
double par_max_ktrmut=0.10;	//max rate for transcriptional mutations

double par_k_tarvol=3.5;	//a number...
double par_genome_to_volume_scale=0.9;
double par_genemut=0.55; // NOT USED//the chances that a gene mutation cause dupl/del ->moved to a cell depd parameter

// parameters that later are controlled by simulation
double par_max_expression_per_gene=0.2;
double par_mut=0.01; //this is squared ONLY for backgr mut
double par_delta=0.1;	//scaling for mutations: par_new += max_par*delta*(genrand_real1() - 0.5)
double BGmut_scheme[N_MUT_TYPES]={  1./3.  ,  1./3.  ,  1./3.  }; //Dupl, Del, Inact
double TRmut_scheme[N_MUT_TYPES]={  1./2.  ,  1./2.  ,  0.  }; //Dupl, Del, Inact

double par_death=0.1; // minimum fraction of population that is changed every generation

//min max metabolic parameters
double par_min_met_par=0.05;
double par_max_met_par=5.5; 

char file_out_var[MAXLEN]="varMicro.txt";
char file_out_par[MAXLEN]="parMicro.txt";
char file_out_post_repl[MAXLEN]="post_replMicro.txt";
char file_out_all[MAXLEN]="var_allMicro.txt";
char file_outpar_all[MAXLEN]="par_allMicro.txt";
char file_out_anctrace[MAXLEN]="anctraceMicro.txt";

char file_backup[MAXLEN]="backupMicro";
char dirname_backup[MAXLEN]="backupMicro";

char file_input[MAXLEN]="";

