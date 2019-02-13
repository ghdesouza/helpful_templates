/**
 * @file graph.tpp
 *
 * @brief Generic Graph Implementation
 *
 * @details This file implement the generic graph optimized for graph with pew modifications in topology during the code execution.
 *
 * @author Gabriel Henrique de Souza (ghdesouza@gmail.com)
 *
 * @date february 13, 2019
 *
 * @copyright Distributed under the Mozilla Public License 2.0 ( https://opensource.org/licenses/MPL-2.0 )
 *
 * @see https://github.com/ghdesouza/helpful_templates
 *
 * Created on: january 15, 2019
 *
 */

#ifndef GRAPH_TPP
#define GRAPH_TPP

#include <stdio.h>

/**
 * @class Graph
 * @brief Generic Graph Implementation
 * @details This class implement the generic graph optimized for graph with pew modifications in topology during the code execution.
 */
template<typename N, typename E>
class Graph{
	
	protected:
		int amount_nodes; /**< Number of nodes in graph */
		
		int id_number; /**< ID of the next node to be created */
		
		int id_edge_number; /**< ID of the next edge to be created */
		
		int *amount_edges; /**< vector with the amount of edge of each node */
		
		N *nodes; /**< Nodes of graph */
		
		E **edges; /**< Edges of graph */
		
		int *id_nodes; /**< list of ID nodes of graph */
		
		int **id_edges; /**< matrix of ID edges of graph */
		
		int **edges_destiny; /**< matrix of destiny of edges in edges matrix */
		
		int nodes_slack; /**< amount of nodes that can be add without change the struct of graph */
		
		int *edges_slack; /**< vector with amount of edges that can be add in each node without change the struct of graph */
		
		int increase_nodes_slack; /**< quantity node increase in graph */
		
		int increase_edges_slack; /**< quantity edges increase in node */
		
		/**
		* @brief Function that extend the max number of nodes in graph
		* @return none
		*/
		void extend_nodes();
		
		/**
		* @brief Function that extend the max number of edges in a node
		* @param id_node - id of node that will increase the max number of edges
		* @return none
		*/
		void extend_edges(int id_node);
		
	public:
	
		/**
		* @brief Constructor
		* @details Create the base strut of start the graph
		* @return none
		*/
		Graph();
		
		/**
		* @brief Destructor
		* @return none
		*/
		~Graph();
		
		/**
		* @brief Function that extend the max number of nodes in graph
		* @param amount_increase - increase size of free space for new nodes in graph
		* @return none
		*/
		void extend_nodes(int amount_increase);
		
		/**
		* @brief Function that add a new node in graph
		* @param node - new node to be add in graph
		* @return id of the add node
		*/
		int add_node(N node);
		
		/**
		* @brief Change a node in graph
		* @param node - new node in graph
		* @param id - id of node that will substituted
		* @return bool that show if find the id node
		*/
		bool set_node(N node, int id);
		
		/**
		* @brief Remove a node in graph
		* @param id - id of node that will remove in graph
		* @return bool that show if find the id node
		*/
		bool remove_node(int id);
		
		/**
		* @brief Function that extend the max number of edges in a node
		* @param id_node - id of node that will increase the max number of edges
		* @param amount_increase - increase size of free space for new edges in node
		* @return none
		*/
		void extend_edges(int id_node, int amount_increase);
		
		/**
		* @brief Function that add a new edge in node
		* @param edge - new edge to be add in node
		* @param id1 - id of start node of edge
		* @param id2 - id of finish node of edge
		* @warning if exist 2 edges with the same start and finish nodes this function will return the first case find
		* @return id of the add edge
		*/
		int add_edge(E edge, int id1, int id2);
		
		/**
		* @brief Change a edge in graph
		* @param edge - new edge in node
		* @param id1 - id of start node of edge
		* @param id2 - id of finish node of edge
		* @warning if exist 2 edges with the same start and finish nodes this function will change the first case find
		* @return bool that show if find the edge
		*/
		bool set_edge(E edge, int id1, int id2);
		
		/**
		* @brief Remove a edge in node
		* @param id1 - id of start node of edge
		* @param id2 - id of finish node of edge
		* @warning if exist 2 edges with the same start and finish nodes this function will remove the first case find
		* @return bool that show if find the id node
		*/
		bool remove_edge(int id1, int id2);
		
		/**
		* @brief Change a edge in graph
		* @param edge - new edge in graph
		* @param id - id of edge that will substituted
		* @return bool that show if find the id edge
		*/
		bool set_edge(E edge, int id_edge);
		
		/**
		* @brief Remove a edge in graph
		* @param id - id of edge that will removed in graph
		* @return bool that show if find the id edge
		*/
		bool remove_edge(int id_edge);
		
		/**
		* @brief Function that return a node
		* @param id - id of node that will be returned
		* @return node finded
		*/
		N get_node(int id);
		
		/**
		* @brief Function that return a edge
		* @param id1 - id of start node of edge
		* @param id2 - id of finish node of edge
		* @warning if exist 2 edges with the same start and finish nodes this function will return the first case find
		* @return edge finded
		*/
		E get_edge(int id1, int id2);
		
		/**
		* @brief Function that return a edge
		* @param id - id of edge that will be returned
		* @return edge finded
		*/
		E get_edge(int id_edge);
		
		/**
		* @brief Function that return the amount of neighbors of a node
		* @param id - id of node that will return the amount of neighbors
		* @return amount neighbors of node
		*/
		int get_amount_neighbors(int id);
		
		/**
		* @brief Function that show the struct of graph for basic tests
		* @warning the type of nodes and edges should be int, float or double
		* @return none
		*/
		void print_graph();

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename N, typename E>
Graph<N, E>::Graph(){

	this->amount_nodes = 0;
	this->increase_nodes_slack = 2;
	this->increase_edges_slack = 2;
	this->id_number = 1;
	this->id_edge_number = 1;
	this->nodes_slack = this->increase_nodes_slack;
	
	this->edges_slack = new int[this->amount_nodes+this->increase_nodes_slack];
	this->amount_edges = new int[this->amount_nodes+this->increase_nodes_slack];
	this->id_nodes = new int[this->amount_nodes+this->increase_nodes_slack];
	this->nodes = new N[this->amount_nodes+this->increase_nodes_slack];
	this->edges_destiny = new int*[this->amount_nodes+this->increase_nodes_slack];
	this->id_edges = new int*[this->amount_nodes+this->increase_nodes_slack];
	this->edges = new E*[this->amount_nodes+this->increase_nodes_slack];
	for(int i = 0; i < this->amount_nodes+this->increase_nodes_slack; i++){
		this->edges_slack[i] = this->increase_edges_slack;
		this->amount_edges[i] = 0;
		this->id_nodes[i] = -1;
		this->edges[i] = new E[this->amount_edges[i]+this->edges_slack[i]];
		this->edges_destiny[i] = new int[this->amount_edges[i]+this->edges_slack[i]];
		this->id_edges[i] = new int[this->amount_edges[i]+this->edges_slack[i]];
		for(int j = 0; j < this->amount_edges[i]+this->edges_slack[i]; j++){
			this->edges_destiny[i][j] = -1;
			this->id_edges[i][j] = -1;
		}
	}

}

template<typename N, typename E>
Graph<N, E>::~Graph(){
	
	for(int i = 0; i < this->amount_nodes+this->nodes_slack; i++){
		delete[] this->edges[i];
		delete[] this->edges_destiny[i];
		delete[] this->id_edges[i];
	}	delete[] this->edges;
		delete[] this->edges_destiny;
		delete[] this->id_edges;
		delete[] this->edges_slack;
		delete[] this->amount_edges;
		delete[] this->id_nodes;
		delete[] this->nodes;
	
}

template<typename N, typename E>
void Graph<N, E>::extend_nodes(int amount_increase){
	
	int* new_edges_slack = new int[this->amount_nodes+amount_increase];
	int* new_amount_edges = new int[this->amount_nodes+amount_increase];
	int* new_id_nodes = new int[this->amount_nodes+amount_increase];
	N* new_nodes = new N[this->amount_nodes+amount_increase];
	int** new_edges_destiny = new int*[this->amount_nodes+amount_increase];
	int** new_id_edges = new int*[this->amount_nodes+amount_increase];
	E** new_edges = new E*[this->amount_nodes+amount_increase];
	
	for(int i = 0; i < this->amount_nodes; i++){
		new_edges_slack[i] = this->increase_edges_slack;
		new_amount_edges[i] = this->amount_edges[i];
		new_id_nodes[i] = this->id_nodes[i];
		new_nodes[i] = this->nodes[i];
		new_edges[i] = new E[this->amount_edges[i]+this->edges_slack[i]];
		new_edges_destiny[i] = new int[this->amount_edges[i]+this->edges_slack[i]];
		new_id_edges[i] = new int[this->amount_edges[i]+this->edges_slack[i]];
		for(int j = 0; j < this->amount_edges[i]; j++){
			new_edges[i][j] = this->edges[i][j];
			new_edges_destiny[i][j] = this->edges_destiny[i][j];
			new_id_edges[i][j] = this->id_edges[i][j];
		}
		for(int j = this->amount_edges[i]; j < this->amount_edges[i]+this->edges_slack[i]; j++){
			new_edges_destiny[i][j] = -1;
			new_id_edges[i][j] = -1;
		}
	}
	
	for(int i = this->amount_nodes; i < this->amount_nodes+amount_increase; i++){
		new_edges_slack[i] = this->increase_edges_slack;
		new_amount_edges[i] = 0;
		new_id_nodes[i] = -1;
		new_edges[i] = new E[this->increase_edges_slack];
		new_edges_destiny[i] = new int[this->increase_edges_slack];
		new_id_edges[i] = new int[this->increase_edges_slack];
		for(int j = 0; j < this->increase_edges_slack; j++){
			new_edges_destiny[i][j] = -1;
			new_id_edges[i][j] = -1;
		}
	}
	
	for(int i = 0; i < this->amount_nodes+this->nodes_slack; i++){
		delete[] this->edges[i];
		delete[] this->edges_destiny[i];
		delete[] this->id_edges[i];
	}	delete[] this->edges;
		delete[] this->edges_destiny;
		delete[] this->id_edges;
		delete[] this->edges_slack;
		delete[] this->amount_edges;
		delete[] this->id_nodes;
		delete[] this->nodes;
	
	this->nodes_slack = amount_increase;
	this->edges = new_edges;
	this->edges_destiny = new_edges_destiny;
	this->id_edges = new_id_edges;
	this->edges_slack = new_edges_slack;
	this->amount_edges = new_amount_edges;
	this->id_nodes = new_id_nodes;
	this->nodes = new_nodes;
}

template<typename N, typename E>
void Graph<N, E>::extend_nodes(){
	this->extend_nodes(this->increase_nodes_slack);
}

template<typename N, typename E>
void Graph<N, E>::extend_edges(int id_node, int amount_increase){
	
	int position = -1;
	for(int i = 0; i < this->amount_nodes; i++) if(this->id_nodes[i] == id_node) position = i;
	if(position < 0) return;
	 
	this->edges_slack[position] = amount_increase;
	E* new_edges = new E[this->amount_edges[position]+this->edges_slack[position]];
	int* new_edges_destiny = new int[this->amount_edges[position]+this->edges_slack[position]];
	int* new_id_edges = new int[this->amount_edges[position]+this->edges_slack[position]];
	for(int j = 0; j < this->amount_edges[position]; j++){
		new_edges[j] = this->edges[position][j];
		new_edges_destiny[j] = this->edges_destiny[position][j];
		new_id_edges[j] = this->id_edges[position][j];
	}
	for(int j = this->amount_edges[position]; j < this->amount_edges[position]+this->edges_slack[position]; j++){
		new_edges_destiny[j] = -1;
		new_id_edges[j] = -1;
	}
	delete[] this->edges[position];
	delete[] this->edges_destiny[position];
	delete[] this->id_edges[position];
	this->edges[position] = new_edges;
	this->edges_destiny[position] = new_edges_destiny;
	this->id_edges[position] = new_id_edges;
}

template<typename N, typename E>
void Graph<N, E>::extend_edges(int id_node){
	this->extend_edges(id_node, this->increase_edges_slack);
}

template<typename N, typename E>
int Graph<N, E>::add_node(N node){
	
	if(this->nodes_slack == 0) this->extend_nodes();
	
	this->id_nodes[this->amount_nodes] = this->id_number;
	this->nodes[this->amount_nodes] = node;
	this->amount_nodes++;
	this->nodes_slack--;
	this->id_number++;
	
	return (this->id_number-1);
}

template<typename N, typename E>
bool Graph<N, E>::set_node(N node, int id){
	
	int initial_position;
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id; initial_position++);
	if(this->id_nodes[initial_position] != id) return false;
	this->nodes[initial_position] = node;
	return true;
}

template<typename N, typename E>
int Graph<N, E>::add_edge(E edge, int id1, int id2){

	int initial_position;
	
	// verify the existence of destiny
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id2; initial_position++);
	if(this->id_nodes[initial_position] != id2) return false;
	
	// find the start of edge
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id1; initial_position++);
	if(this->id_nodes[initial_position] != id1) return false;
	
	if(this->edges_slack[initial_position] == 0) this->extend_edges(initial_position+1);
	
	this->edges[initial_position][this->amount_edges[initial_position]] = edge;
	this->edges_destiny[initial_position][this->amount_edges[initial_position]] = id2;
	this->id_edges[initial_position][this->amount_edges[initial_position]] = this->id_edge_number;
	this->amount_edges[initial_position]++;
	this->edges_slack[initial_position]--;
	this->id_edge_number++;
	
	return (this->id_edge_number-1);
}

template<typename N, typename E>
bool Graph<N, E>::set_edge(E edge, int id1, int id2){
	
	int initial_position, last_position;
		
	// find the start of edge
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id1; initial_position++);
	if(this->id_nodes[initial_position] != id1) return false;
	
	// find the finish of edge
	for(last_position = 0; last_position < this->amount_edges[initial_position] && this->edges_destiny[initial_position][last_position] != id2; last_position++);
	if(this->edges_destiny[initial_position][last_position] != id2) return false;
	
	this->edges[initial_position][last_position] = edge;
	return true;
}

template<typename N, typename E>
bool Graph<N, E>::set_edge(E edge, int id_edge){
	
	int initial_position, last_position;
	
	for(initial_position = 0; initial_position < this->amount_nodes; initial_position++){
		for(last_position = 0; last_position < this->amount_edges[initial_position]; last_position++){
			if(this->id_edges[initial_position][last_position] == id_edge){
				this->edges[initial_position][last_position] = edge;
				return true;
			}
		}
	}
	return false;
}

template<typename N, typename E>
bool Graph<N, E>::remove_edge(int id1, int id2){

	int initial_position, last_position;
		
	// find the start of edge
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id1; initial_position++);
	if(this->id_nodes[initial_position] != id1) return false;
	
	// find the finish of edge
	for(last_position = 0; last_position < this->amount_edges[initial_position] && this->edges_destiny[initial_position][last_position] != id2; last_position++);
	if(this->edges_destiny[initial_position][last_position] != id2) return false;
	
	E temp = this->edges[initial_position][last_position];
	this->edges[initial_position][last_position] = this->edges[initial_position][this->amount_edges[initial_position]-1];
	this->edges[initial_position][this->amount_edges[initial_position]-1] = temp;
	this->edges_destiny[initial_position][last_position] = this->edges_destiny[initial_position][this->amount_edges[initial_position]-1];
	this->edges_destiny[initial_position][this->amount_edges[initial_position]-1] = -1;
	this->id_edges[initial_position][last_position] = this->id_edges[initial_position][this->amount_edges[initial_position]-1];
	this->id_edges[initial_position][this->amount_edges[initial_position]-1] = -1;
	this->amount_edges[initial_position]--;
	this->edges_slack[initial_position]++;
	
	return true;
}

template<typename N, typename E>
bool Graph<N, E>::remove_node(int id){
	
	int initial_position;
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id; initial_position++);
	if(this->id_nodes[initial_position] != id) return false;
	
	for(int i = 0; i < this->amount_nodes; i++){
		while(this->remove_edge(this->id_nodes[i], id));
	}
	
	N temp = this->nodes[this->amount_nodes-1];
	this->nodes[this->amount_nodes-1] = this->nodes[initial_position];
	this->nodes[initial_position] = temp;
	
	E* temp2 = this->edges[this->amount_nodes-1];
	this->edges[this->amount_nodes-1] = this->edges[initial_position];
	this->edges[initial_position] = temp2;
	
	int temp3 = this->id_nodes[this->amount_nodes-1];
	this->id_nodes[this->amount_nodes-1] = this->id_nodes[initial_position];
	this->id_nodes[initial_position] = temp3;
	this->id_nodes[this->amount_nodes-1] = -1;
	
	this->amount_nodes--;
	this->nodes_slack++;
	return true;
}

template<typename N, typename E>
bool Graph<N, E>::remove_edge(int id_edge){
	
	int initial_position, last_position;
		
	for(initial_position = 0; initial_position < this->amount_nodes; initial_position++){
		for(last_position = 0; last_position < this->amount_edges[initial_position]; last_position++){
			if(this->id_edges[initial_position][last_position] == id_edge){

				E temp = this->edges[initial_position][last_position];
				this->edges[initial_position][last_position] = this->edges[initial_position][this->amount_edges[initial_position]-1];
				this->edges[initial_position][this->amount_edges[initial_position]-1] = temp;
				this->edges_destiny[initial_position][last_position] = this->edges_destiny[initial_position][this->amount_edges[initial_position]-1];
				this->edges_destiny[initial_position][this->amount_edges[initial_position]-1] = -1;
				this->id_edges[initial_position][last_position] = this->id_edges[initial_position][this->amount_edges[initial_position]-1];
				this->id_edges[initial_position][this->amount_edges[initial_position]-1] = -1;
				this->amount_edges[initial_position]--;
				this->edges_slack[initial_position]++;
				return true;
			}
		}
	}
	
	return false;
}

template<typename N, typename E>
N Graph<N, E>::get_node(int id){
	
	int initial_position;
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id; initial_position++);
	if(this->id_nodes[initial_position] != id){
		printf("\n\nERROR: ID of node not found!\n\n");
		return this->nodes[0];
	}
	return this->nodes[initial_position];
}

template<typename N, typename E>
E Graph<N, E>::get_edge(int id1, int id2){
	
	int initial_position, last_position;
		
	// find the start of edge
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id1; initial_position++);
	if(this->id_nodes[initial_position] != id1){
		printf("\n\nERROR: ID of node not found!\n\n");		
	}
	
	// find the finish of edge
	for(last_position = 0; last_position < this->amount_edges[initial_position] && this->edges_destiny[initial_position][last_position] != id2; last_position++);
	if(this->edges_destiny[initial_position][last_position] != id2){
		printf("\n\nERROR: ID of node not found!\n\n");		
	}
	
	return this->edges[initial_position][last_position];
}

template<typename N, typename E>
E Graph<N, E>::get_edge(int id_edge){

	int initial_position, last_position;
	
	for(initial_position = 0; initial_position < this->amount_nodes; initial_position++){
		for(last_position = 0; last_position < this->amount_edges[initial_position]; last_position++){
			if(this->id_edges[initial_position][last_position] == id_edge) return this->edges[initial_position][last_position];
		}
	}
	printf("\n\nERROR: ID of node not found!\n\n");
	return NULL;
}

template<typename N, typename E>
int Graph<N, E>::get_amount_neighbors(int id){
	
	int initial_position;		
	for(initial_position = 0; initial_position < this->amount_nodes && this->id_nodes[initial_position] != id; initial_position++);
	if(this->id_nodes[initial_position] != id) return 0;
	
	return this->amount_edges[initial_position];
}

template<typename N, typename E>
void Graph<N, E>::print_graph(){

	for(int i = 0; i < this->amount_nodes; i++){
		printf("ID node: %d (%.3e)\n\tEdges: ", this->id_nodes[i], (double) this->nodes[i]);
		for(int j = 0; j < this->amount_edges[i]; j++){
			printf("(%d, %d, %.2e) ", this->id_nodes[i], this->edges_destiny[i][j], (double) this->edges[i][j]);
		}
		printf("\n");
	}
	
}

#endif // GRAPH_TPP

