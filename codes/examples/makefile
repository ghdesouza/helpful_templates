DE = main_de
GRAPH = main_graph

CC = g++
CFLAGS = -g -O3 -Wall -fopenmp

de: $(DE)
graph: $(GRAPH)

all: $(DE) $(GRAPH)

$(DE): $(DE).cpp
	$(CC) $(CFLAGS) $(DE).cpp -o bin/$(DE).out

$(GRAPH): $(GRAPH).cpp
	$(CC) $(CFLAGS) $(GRAPH).cpp -o bin/$(GRAPH).out
