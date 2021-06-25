#ifndef CELL_H
#define CELL_H

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "parameters.h"
#include "my.h"

void Initialise_gsl(int dimension, CELL *world, gsl_odeiv2_driver *d);

//void IntegrateOneStep(CELL *ic, double delta_time, double *S_pointer, double *ds_pointer, double *Sin_pointer);
void IntegrateOneStep(double Time, double *t,  double **y, gsl_odeiv2_driver * d);
void CellVolume(CELL *icel);
void CellHealth(CELL *icel);
void TargetVolume(CELL *ic);
void GenomeStats(CELL *ic);

void CopyCellHereToThere(CELL *world,double **y, int parent_pos, int child_pos);
void MutateCell(CELL *child);
void HalfCellContent(CELL *world,double **y, int pos);

double Replication(CELL *world,double **y,int cells_status,double max_growthtime);
void TranscrMutations(CELL *child);
void TranscrMutations_diffratesDupDelIn(CELL *child);
void MutMutations(CELL *child);
void MutationsMetabolism(CELL *ic);
void RegMutate(double *regpar);
void MutationsRegulation(CELL *ic);
void MetMutate(double *metpar);
void BackgroundMutations(CELL *ic);
int SingleGeneMut(CELL *ic,int mutpos,double *mut_scheme); //returns the type of mutation (needed by TranscrMutations())

// GRAVEYARD
// void DivideCell(CELL *world,double **y, int parent_pos, int child_pos);
// int Replication_BACKUP(CELL *world,double **y,int cells_status,double Time);
// void DivideCell_BACKUP(CELL *world,double **y, int parent_pos, int child_pos);
// int ReplicationPop2(CELL *world,double **py,int cells_status);
// void DivideCellPop2(CELL *world,double **py, int cell_pos);
// void TranscrMutations_BACKUP(CELL *child,CELL *parent);
// void TranscrMutations_OLD(CELL *child,CELL *parent);

#endif
