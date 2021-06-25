#something of a makefile

PROJECT = pMicro

CC = /usr/bin/gcc

DEPS = my.h constants.h parameters.h cell.h makefile mersenne.h data_output.h initialisation.h
OBJ = main.o cell.o parameters.o mersenne.o data_output.o initialisation.o

#CCOPT = -O0 -ggdb -Wall #-pedantic -Wextra # For debugging.
CCOPT = -O3 -Wall #Most optimised

LDIR = -L/usr/X11R6/lib
LIBS = -lz -lm -lgsl -lgslcblas 

all: $(OBJ)
	$(CC) $(OBJ) $(CCOPT) -o $(PROJECT) $(LIBS) $(LDIR) 

%.o: %.c $(DEPS)
	$(CC) $(CCOPT) -c $< -o $@


source.tar.gz: $(OBJ) $(DEPS)
	tar -zcf source.tar.gz *.c *.h makefile *py README


