#include "my.h"

int InitialiseFromScratch(CELL *world,double **y);
int InitialiseFromBackup(CELL *world, double **y);
void InitialiseAncestry(CELL *ic, int i);
void CheckOutputFilesExist(void);
void CheckEnvironmentSwitches(void);