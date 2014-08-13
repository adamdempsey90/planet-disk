EXECUTABLE=planetdisk
SOURCES=algo.c fourier.c boundary.c init.c main.c output.c readinputs.c utils.c viscosity.c
INPUT=params.opt
HEADER1=defines.h
HEADER2=planetdisk.h
LDFLAGS=-lgsl -lgslcblas -lgomp  -lfftw3 -lm 

CFLAGS=-c -Wall -O3 -fopenmp -Iinc

BIN=bin/
SRC=src/
IN=inputs/
PY=src/pyutils/

CC=gcc

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CINPUT=$(addprefix $(IN),$(INPUT))
CHEADER1=$(addprefix $(SRC),$(HEADER1))
CHEADER2=$(addprefix $(SRC),$(HEADER2))


all: $(CSOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(COBJECTS) 
	$(CC)  $(COBJECTS) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER1) $(CHEADER2)
	$(CC) $(CFLAGS) $< -o $@

$(CHEADER1): $(CINPUT)
	python $(PY)defines.py 
	
clean:
	rm $(COBJECTS) $(EXECUTABLE)
