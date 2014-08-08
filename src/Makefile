CC=gcc-4.9
CFLAGS= 
LDFLAGS=-lgsl -lgslcblas -lgomp -fopenmp -lfftw3 -lm 
SOURCES=algo.c fourier.c init.c main.c output.c readinputs.c utils.c viscosity.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=planetdisk

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS) $(EXECUTABLE)


