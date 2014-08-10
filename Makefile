EXECUTABLE=planetdisk
SOURCES=algo.c fourier.c init.c main.c output.c readinputs.c utils.c viscosity.c

LDFLAGS=-lgsl -lgslcblas -lgomp -fopenmp -lfftw3 -lm 

CFLAGS=-c -Iinc

BIN=bin/
SRC=src/

CC=gcc

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))

all: $(CSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS)
	$(CC) $(LDFLAGS) $(COBJECTS) -o $@

$(BIN)%.o: $(SRC)%.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(COBJECTS) $(EXECUTABLE)