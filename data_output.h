#ifndef DATAOUT_H
#define DATAOUT_H

#include "my.h"

void PrintACell(int Time, CELL *ic);
void SaveData(double Time,CELL *world);
void SaveData_PostReplication(double Time, CELL *world);
void OutputBackup(int Time, CELL *world, CELL *prev_world);
void SaveBackup(int Time, CELL *world, char *fname, CELL *prev_world);
void SaveDataRelative(double Time,CELL *world);
void SaveDataRelative2(int Time,CELL *world);
void SaveDataRelative3(int Time,CELL *world);
void SaveDataRelative4(double Time,CELL *world);
#endif