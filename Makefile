EXECUTABLE=planetdisk
SOURCES=algo.c fourier.c boundary.c init.c main.c output.c readinputs.c utils.c viscosity.c rk45.c rk45_step.c
HEADER=planetdisk.h defines.h rk45.h

LDFLAGS=-lfftw3 -lm 

CFLAGS=-c -Wall -O3 -Iinc

BIN=bin/
SRC=src/
IN=inputs/
PY=src/pyutils/

CC=mpicc

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(COBJECTS) 
	$(CC)  $(COBJECTS) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER) 
	$(CC) $(CFLAGS) $< -o $@

	
clean:
	rm $(COBJECTS) $(EXECUTABLE) 
