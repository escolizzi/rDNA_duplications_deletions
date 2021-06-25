#ifndef CONST_H
#define CONST_H

// GENERAL CONSTANTS
#define MAXLEN 256
#define MAXGENOMESIZE 10000
#define MAX_GROWTH_TIME 10000.
#define EPSILON 0.0001	// a small number
#define MAX_PAR_EXPRESSION 10.

// VARIABLES
#define N_VARIABLES 9 // number of variables
//resources
#define _S 0
//metabolites
#define _A 1
//metabolism
#define _T 2
// housekeeping
#define _Q 3
//translation
//  mRNAs
#define _mT 4
#define _mQ 5
#define _mRp 6
// ribosome and complexes
#define _Rr 7
#define _Rp 8
// #define _Crmt 7
// #define _Crmq 8


// number of metabolic costants, excluding decay rates.
#define N_METAB_PAR 3 // kt,alfat,betat, 
#define N_GENE_TYPES 4	//number of gene types: 
			// T,Q,Rp,Rr
#define GENE_T 0
#define GENE_Q 1
#define GENE_Rr 2
#define GENE_Rp 3

//Background mutations
#define N_MUT_TYPES 3

#define DUPLIC 0
#define DELET 1
#define INACT 2

#endif
