# Makefile for the SLANG - MPI - Module
# Written by Thomas Dauser & Fritz Schwarm
#
# Revision:	0 - First
# Version:	1

REV = 0
VER = 1

CC = mpicc
CC_OPTS = -O0 -g -fPIC

SRCDIR = ./
INCDIR = ./
LIBDIR = ./

INC_OPTS = -I${INCDIR}
LIB_OPTS = -L${LIBDIR}

OBJS += $(patsubst %.c, %.o, $(wildcard *.c))
HDRS += $(wildcard *.h)
OPTS += $(INC_OPTS) -lm
MPICC_LOC=/usr/lib64/openmpi/bin/mpicc

MOD_NAME = slmpi

all:		${OBJS} ${HDRS}
		@echo
		@echo BUILDING SLANG MPI MODULE...
		@make --no-print-directory wrapper
		@if [ -f test/test.c ]; then make -s test; fi
		@echo
		@echo TESTING SLANG MPI MODULE...
		@echo
		@make --no-print-directory -f $(MOD_NAME).mf test


.PHONY:		wrapper
wrapper:	
		slirp -make $(MOD_NAME).h $(MOD_NAME).o
		sed -i 's|'`which gcc`'|'$(MPICC_LOC)'|' $(MOD_NAME).mf
		@make --no-print-directory -f $(MOD_NAME).mf

test:		test/test.o
		${CC} ${CC_OPTS} ${INC_OPTS} ${LIB_OPTS} test/test.o -o test/test

%.o:		%.c
		${CC} ${CC_OPTS} -c $< -o $@ $(INC_OPTS)

.PHONY:		clean
clean:
		@OUTPUT=$(shell rm test/test test/*.o test/*~ *.o $(MOD_NAME)_glue.c $(MOD_NAME)-test.sl $(MOD_NAME)-module.so *~ 2> /dev/null)
