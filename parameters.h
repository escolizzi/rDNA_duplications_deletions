#ifndef PARA_H
#define PARA_H

#include "constants.h"

extern int par_initTime;
extern int par_maxTime;
extern double par_init_interval_time;
extern long int par_ulseed;
extern int par_time_savedata;
extern int par_time_backup;
extern int par_time_anc;
extern double par_period_time_change;

extern const int MUT_REGULATION;
extern const int MUT_TRANSCR_LOAD;
extern const int MUT_BACKGROUND;
extern const int MUT_MUTATION;

extern const int ENV_S_RAND;
extern const int ENV_S_INCR;
extern const int ENV_S_CONST;
extern const int AVRG_ENV;

extern const int par_pop;
extern char init_genome[MAXGENOMESIZE];
extern double par_init_tarvol;
extern double par_init_minvol;
extern double par_tar_Q;
extern double transcr_penalty;
extern double par_ds;
extern double par_Sin;
extern int par_time_change_Sin;

extern double par_DNAbeta;
extern double par_ktrmut;
extern double par_max_ktrmut;
extern double par_k_tarvol;
extern double par_genome_to_volume_scale;
extern double par_genemut;
extern double par_max_expression_per_gene;
extern double par_mut;
extern double par_delta;
extern double BGmut_scheme[N_MUT_TYPES];
extern double TRmut_scheme[N_MUT_TYPES];
extern double par_death;
extern double par_min_met_par;
extern double par_max_met_par;

extern char file_out_var[MAXLEN];
extern char file_out_par[MAXLEN];
extern char file_out_post_repl[MAXLEN];
extern char file_out_all[MAXLEN];
extern char file_outpar_all[MAXLEN];
extern char file_out_anctrace[MAXLEN];

extern char file_backup[MAXLEN];
extern char dirname_backup[MAXLEN];

extern char file_input[MAXLEN];
#endif
