DE = main_de
GRAPH = main_graph
CLASSIFIERS = main_classifiers

CC = g++
CFLAGS = -g -Wall -fopenmp

de: $(DE)
graph: $(GRAPH)
classifiers: $(CLASSIFIERS)

all: $(DE) $(GRAPH)

$(DE): $(DE).cpp
	$(CC) $(CFLAGS) $(DE).cpp -o bin/$(DE).out

$(GRAPH): $(GRAPH).cpp
	$(CC) $(CFLAGS) $(GRAPH).cpp -o bin/$(GRAPH).out

$(CLASSIFIERS): $(CLASSIFIERS).cpp
	$(CC) $(CFLAGS) $(CLASSIFIERS).cpp -o bin/$(CLASSIFIERS).out
