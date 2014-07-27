CC=gcc
CFLAGS=-c -Wall -O3 
LDFLAGS=-lgsl -lgslcblas -lm -lgomp -fopenmp
SOURCES=algo.c derivs.c output.c utils.c amf.c init.c readinputs.c viscosity.c boundary.c main.c restart.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=meanwave

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS) $(EXECUTABLE)

