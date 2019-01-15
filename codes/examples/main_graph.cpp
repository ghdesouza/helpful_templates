#include "../Graph.tpp"
#include <stdio.h>
#include <iostream>

using namespace std;

void load_graph(Graph<double, double> *graph);
void load_graph_file(string nome_file, Graph<double, double> *graph);

int main(){

	Graph<double, double> *graph = new Graph<double, double>();
	load_graph_file("datas/graph.txt", graph);
	//load_graph(graph);
	
	graph->remove_edge(1, 2);
	graph->remove_edge(4);
	graph->remove_node(3);
	graph->set_node(88.02, 6);
	graph->set_edge(12.1, 2, 5);
	graph->set_edge(7.77, 1);
	graph->print_graph();

	delete graph;
    return 0;
}

void load_graph(Graph<double, double> *graph){

	graph->add_node(1.5);
	graph->add_node(2.5);
	graph->add_node(3.5);
	graph->add_node(4.5);
	graph->add_node(5.5);
	graph->add_node(6.5);
	
	graph->add_edge(1, 1, 1); graph->add_edge(2, 1, 2); graph->add_edge(3, 1, 3); graph->add_edge(4, 1, 4); graph->add_edge(5, 1, 5); graph->add_edge(6, 1, 6);
	graph->add_edge(7, 2, 1); graph->add_edge(8, 2, 2); graph->add_edge(9, 2, 3); graph->add_edge(10, 2, 4); graph->add_edge(11, 2, 5); graph->add_edge(12, 2, 6);
	graph->add_edge(13, 3, 1); graph->add_edge(14, 3, 2); graph->add_edge(15, 3, 3); graph->add_edge(16, 3, 4); graph->add_edge(17, 3, 5); graph->add_edge(18, 3, 6);
	graph->add_edge(19, 4, 1); graph->add_edge(20, 4, 2); graph->add_edge(21, 4, 3); graph->add_edge(22, 4, 4); graph->add_edge(23, 4, 5); graph->add_edge(24, 4, 6);
	graph->add_edge(25, 5, 1); graph->add_edge(26, 5, 2); graph->add_edge(27, 5, 3); graph->add_edge(28, 5, 4); graph->add_edge(29, 5, 5); graph->add_edge(30, 5, 6);
	graph->add_edge(31, 6, 1); graph->add_edge(32, 6, 2); graph->add_edge(33, 6, 3); graph->add_edge(34, 6, 4); graph->add_edge(35, 6, 5); graph->add_edge(36, 6, 6);
}

void load_graph_file(string nome_file, Graph<double, double> *graph){
	FILE* arq = fopen(nome_file.c_str(), "r");
	int nodes_amount;
	double weigth;
	int id_node, id_destiny, amount_edges;
	if(!fscanf(arq, "%d", &nodes_amount)) return;
	graph->extend_nodes(nodes_amount);
	for(int i = 0; i < nodes_amount; i++){
		if(!fscanf(arq, "%lf", &weigth)) return;
		graph->add_node(weigth);
	}
	if(!fscanf(arq, "%d", &id_node)) return;
	while(id_node != 0){
		if(!fscanf(arq, "%d", &amount_edges)) return;
		graph->extend_edges(id_node, amount_edges);
		for(int i = 0; i < amount_edges; i++){
			if(!fscanf(arq, "%d", &id_destiny) || 
			!fscanf(arq, "%lf", &weigth)) return;
			graph->add_edge(weigth, id_node, id_destiny);
		}
		if(!fscanf(arq, "%d", &id_node)) return;
	}
	fclose(arq);
}

