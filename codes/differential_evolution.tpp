/**
 * @file differential_evolution.tpp
 *
 * @brief Differential Evolution Algorithm
 *
 * @details This file implement the Differential Evolution Algorithm for generic use includes constraint and local search problems.
 *
 * @author Gabriel Henrique de Souza (ghdesouza@gmail.com)
 *
 * @date September 17, 2019
 *
 * @copyright Distributed under the Mozilla Public License 2.0 ( https://opensource.org/licenses/MPL-2.0 )
 *
 * @see https://github.com/ghdesouza/helpful_templates
 *
 * Created on: january 15, 2019
 *
 */

#ifndef DIFFERENTIAL_EVOLUTION_TPP
#define DIFFERENTIAL_EVOLUTION_TPP

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

#include <iostream>
using namespace std;

/**
 * @class Differential_Evolution
 * @brief Differential Evolution Algorithm
 * @details This class implement the Differential Evolution Algorithm for generic use includes constraint and local search problems.
 */
template<typename T = float>
class Differential_Evolution{

    protected:

		int dimension; /**< Dimension of agent */
		
		int population_size; /**< Population Size */
		
		int best_agent_position; /**< Indice of the best agent of population */
		
		int amount_fitness_evaluation; /**< Amount of Fitness Evaluation */
		
		T** population; /**< Agents of population */
		
		T** temporaly_population; /**< Population created before the fitness avaliation */
		
		T* fitness; /**< Fitness values of population */
		
		T CR; /**< Crossover Probability */
		
		T F; /**< Differential Weight */
		
		T LS; /**< Local Search Probability */
		
		int mutation_type; /**< Type of mutation that is used in DE: 0 - rand_1, 1 - best_1, 2 - target_to_rand, 3 - target_to_best, 4 - rand_2, 5 - best_2 */
		
		bool agent_repair; /**< show if exist agent_repair_function */
		
		bool parallel; /**< show if next_generation is parallelized */
		
		T (*fitness_function)(T* agent); /**< Pointer for Fitness Function */
		
		void (*agent_repair_function)(T* agent); /**< Pointer for Agent Repair Function */

		void (*local_search_function)(T* agent, T* new_agent); /**< Pointer for Local Search Function */
		
		/**
		* @brief Function that create a base agent
		* @details initializes the agent with randon values in the range of -1 to 1 
		* @param agent - pointer for a new agent
		* @return none
		*/
		void create_agent_basic(T* agent);
		
		/**
		* @brief Function that create a mutated agent based in rand 1 operator
		* @details Move the agent based in other 3 differentes agents using the equation:
		* Y = (A + F x (B-C)) where A, B, and C is a 3 different random agents
		* The new agent will saved in a private vector of class "temporaly_population" in same position that agent_x.
		* @param agent_x - Agent position that will moved
		* @return none
		*/
		void mutation_rand_1(int agent_x);
		void mutation_best_1(int agent_x);		
		void mutation_target_rand(int agent_x);		
		void mutation_target_best(int agent_x);		
		void mutation_rand_2(int agent_x);		
		void mutation_best_2(int agent_x);		
		
		/**
		* @brief Function to do a crossover between the population and the mutation population
		* @details Assembly between the populations
		* @param agent_x - Agent position
		* @return none
		*/
		void crossover(int agent_x);
		
		/**
		* @brief Try to move the agent of temporaly population to population
		* @details Verify if the fitness function is improved in temporaly population so move the agent for population
		* and if is the new best fitness update the best_agent_position.
		* @param agent_x - Agent position that will moved
		* @return none
		*/
		void selection(int agent_x);
		
    public:

		/**
		* @brief Constructor
		* @param dimension - Size of agent in DE
		* @param population_size - number of agents in population
		* @param fitness_function - pointer to function that avaliable an agent
		* @param F - differential weigth
		* @param CR - crossover probability
		* @return none
		*/
        Differential_Evolution(int dimension, int population_size, T (*fitness_function)(T* agent), T F, T CR);
        Differential_Evolution(){}
		
		/**
		* @brief Destructor
		* @return none
		*/
        ~Differential_Evolution();
		
		/**
		* @brief Function that return the size of agent
		* @return size of agent
		*/
		int get_dimension(){return this->dimension;}

		/**
		* @brief Function that return the number of agents in population
		* @return population size
		*/
		int get_population_size(){return this->population_size;}
		
		/**
		* @brief Function that return the crossover probability
		* @return crossover probability
		*/
		T get_CR(){return this->CR;}
		
		/**
		* @brief Function that return the differential weigth
		* @return differential weigth
		*/
		T get_F(){return this->F;}
		
		/**
		* @brief Function that return the best fitness value in population
		* @return best fitness in population
		*/
		T get_best_fitness(){return this->fitness[this->best_agent_position];}
		
		/**
		* @brief Function that return the mean fitness values of population
		* @return mean fitness
		*/
		T get_mean_fitness();		
		
		/**
		* @brief Function that return the standard deviation of fitness values of population
		* @return standard deviation of fitness values
		*/
		T get_convergence_fitness();
		
		/**
		* @brief Function that return the standard deviation of ith dimension
		* @return standard deviation of ith dimension
		*/
		T get_convergence_dimension(int dimension);
		
		/**
		* @brief Function that return the mean of standard deviation of dimensions
		* @return mean of standard deviation of dimensions
		*/
		T get_convergence_population();
		
		/**
		* @brief Function that return the amount of fitness function is used
		* @return amount of fitness function is used
		*/
		int get_amount_fitness_evaluation(){return this->amount_fitness_evaluation;}
		
		/**
		* @brief Function that set amount of fitness function in zero
		* @return none
		*/
		void clean_amount_fitness_evaluation(){this->amount_fitness_evaluation = 0;}
		
		/**
		* @brief Set the function that will repair a invalid agent
		* @param agent_repair_function - pointer to function that will repair a invalid agent
		* @return none
		*/
		void set_agent_repair(void (*agent_repair_function)(T* agent));		
		
		/**
		* @brief Set the local search function
		* @param local_search_function - pointer to local search function
		* @param LS - probability of local search
		* @return none
		*/
		void set_local_search(void (*local_search_function)(T* agent, T* new_agent), T LS);
		
		/**
		* @brief Set the crossover probability
		* @param CR - new crossover probability
		* @return none
		*/
		void set_CR(T CR){this->CR = CR;}
		
		/**
		* @brief Set the differential weigth
		* @param F - new differential weigth
		* @return none
		*/
		void set_F(T F){this->F = F;}
		
		/**
		* @brief Set the mutation type
		* @param type - type of mutation: rand1, best1, target_rand, target_best, rand2, best2
		* @return none
		*/
		void set_mutation_type(string type);
		
		/**
		* @brief Set if next generation is parallelized
		* @param parallel - bool that show if next_generation is parallelized
		* @return none
		*/
		void set_parallel(bool parallel);
		
		/**
		* @brief Create a new population
		* @details Create a new population based in a function
		* @param create_agent_function - pointer to function that create an agent
		* @return none
		*/
		void new_population(void (*create_agent_function)(T* agent));
		
		/**
		* @brief Create a new population
		* @details Create a new population with agents initializes with randon values in the range of -1 to 1
		* @return none
		*/
		void new_population();
		
		/**
		* @brief Differential Generation
		* @details Create the next generation
		* @return none
		*/
		void next_generation();

		/**
		* @brief Differential Evolution
		* @details Population Evolution many genarations
		* @param stoped_value - value of fitness function that stop the evolution
		* @param max_stoped_generation - max generation without improve the mean fitness function in population
		* @param max_generation - max number of generations
		* @return none
		*/		
		void evolution(T stop_value, int max_stoped_generation, int max_generation);
		
		/**
		* @brief Clone the best agent in population
		* @param agent - pointer for to save the best agent
		* @return none
		*/
		void clone_best_agent(T* agent);

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
Differential_Evolution<T>::Differential_Evolution(int dimension, int population_size, T (*fitness_function)(T* agent), T F, T CR){

	this->dimension = dimension;
	this->population_size = population_size;
	if(this->population_size < 4) this->population_size = 4;
	this->fitness_function = fitness_function;
	this->F = F;
	this->CR = CR;
	this->LS = 0;
	this->mutation_type = 0;
	this->amount_fitness_evaluation = 0;
	this->agent_repair = false;
	this->parallel = false;
	
	this->fitness = new T[this->population_size];
	this->population = new T*[this->population_size];
	this->temporaly_population = new T*[this->population_size];
	for(int i = 0; i < this->population_size; i++){
		this->population[i] = new T[this->dimension];
		this->create_agent_basic(this->population[i]);
		this->temporaly_population[i] = new T[this->dimension];
	}
	this->fitness[0] = this->fitness_function(this->population[0]);
	this->amount_fitness_evaluation++;
	this->best_agent_position = 0;
	for(int i = 1; i < this->population_size; i++){
		this->fitness[i] = this->fitness_function(this->population[i]);
		this->amount_fitness_evaluation++;
		if(this->fitness[i] < this->fitness[this->best_agent_position]) this->best_agent_position = i;		
	}
	
}

template<typename T>
Differential_Evolution<T>::~Differential_Evolution(){

	for(int i = 0; i < this->population_size; i++){
		delete[] this->population[i];
		delete[] this->temporaly_population[i];
	}	delete[] this->population;
		delete[] this->temporaly_population;
		delete[] this->fitness;

}

template<typename T>
T Differential_Evolution<T>::get_mean_fitness(){
	T exit = 0;
	for(int i = 0; i < this->population_size; i++) exit += this->fitness[i];
	return exit/this->population_size;
}

template<typename T>
T Differential_Evolution<T>::get_convergence_fitness(){
	T mean = this->get_mean_fitness();
	T exit = 0;
	for(int i = 0; i < this->population_size; i++) exit += pow(this->fitness[i]-mean, 2);
	return sqrt(exit/this->population_size);
}

template<typename T>
T Differential_Evolution<T>::get_convergence_dimension(int dimension){
	T mean = 0;
	for(int i = 0; i < this->population_size; i++) mean += this->population[i][dimension];
	mean /= this->population_size;
	
	T exit = 0;
	for(int i = 0; i < this->population_size; i++) exit += pow(this->population[i][dimension]-mean, 2);
	return sqrt(exit/this->population_size);
}

template<typename T>
T Differential_Evolution<T>::get_convergence_population(){
	T mean = 0;
	for(int i = 0; i < this->dimension; i++) mean += this->get_convergence_dimension(i);
	mean /= this->dimension;
	
	return mean;
}

template<typename T>
void Differential_Evolution<T>::set_parallel(bool parallel){
	
	#ifdef _OPENMP
		this->parallel = parallel;
	#else
		this->parallel = false;
	#endif
	
}

template<typename T>
void Differential_Evolution<T>::set_agent_repair(void (*agent_repair_function)(T* agent)){
	this->agent_repair = true;
	this->agent_repair_function = agent_repair_function;
}

template<typename T>
void Differential_Evolution<T>::set_local_search(void (*local_search_function)(T* agent, T* new_agent), T LS){
	this->LS = LS;
	this->local_search_function = local_search_function;
}

template<typename T>
void Differential_Evolution<T>::set_mutation_type(string type){
	if(type == "rand1")	this->mutation_type = 0;
	else if(type == "best1") this->mutation_type = 1;
	else if(type == "target_rand") this->mutation_type = 2;
	else if(type == "target_best") this->mutation_type = 3;
	else if(type == "rand2") this->mutation_type = 4;
	else if(type == "best2") this->mutation_type = 5;
	else printf("ERROR (DE): invalid type of mutation");
	return;
}

template<typename T>
void Differential_Evolution<T>::mutation_rand_1(int agent_x){

	int a, b, c;
	a = rand()%this->population_size;
	b = rand()%this->population_size;
	c = rand()%this->population_size;
	while(b == a) b = (b+1)%this->population_size;
	while(c == a || c == b) c = (c+1)%this->population_size;
	for(int i = 0; i < this->dimension; i++){
		this->temporaly_population[agent_x][i] = this->population[a][i]
												 + this->F*(this->population[b][i]-this->population[c][i]);
	}	
}

template<typename T>
void Differential_Evolution<T>::mutation_best_1(int agent_x){

	int a, b;
	a = rand()%this->population_size;
	b = rand()%this->population_size;
	while(b == a) b = (b+1)%this->population_size;
	for(int i = 0; i < this->dimension; i++){
		this->temporaly_population[agent_x][i] = this->population[this->best_agent_position][i] 
												 + this->F*(this->population[a][i]-this->population[b][i]);
	}
}

template<typename T>
void Differential_Evolution<T>::mutation_target_rand(int agent_x){

	int a, b, c;
	a = rand()%this->population_size;
	b = rand()%this->population_size;
	c = rand()%this->population_size;
	while(b == a) b = (b+1)%this->population_size;
	while(c == a || c == b) c = (c+1)%this->population_size;
	for(int i = 0; i < this->dimension; i++){
		this->temporaly_population[agent_x][i] = this->population[agent_x][i]
												 + this->F*(this->population[a][i]-this->population[agent_x][i])
												 + this->F*(this->population[b][i]-this->population[c][i]);
	}	
}

template<typename T>
void Differential_Evolution<T>::mutation_target_best(int agent_x){

	int a, b;
	a = rand()%this->population_size;
	b = rand()%this->population_size;
	while(b == a) b = (b+1)%this->population_size;
	for(int i = 0; i < this->dimension; i++){
		this->temporaly_population[agent_x][i] = this->population[agent_x][i]
												 + this->F*(this->population[this->best_agent_position][i]-this->population[agent_x][i])
												 + this->F*(this->population[a][i]-this->population[b][i]);
	}	
}

template<typename T>
void Differential_Evolution<T>::mutation_rand_2(int agent_x){

	int a, b, c, d, e;
	a = rand()%this->population_size;
	b = rand()%this->population_size;
	c = rand()%this->population_size;
	d = rand()%this->population_size;
	e = rand()%this->population_size;
	while(b == a) b = (b+1)%this->population_size;
	while(c == a || c == b) c = (c+1)%this->population_size;
	while(d == a || d == b || d == c) d = (d+1)%this->population_size;;
	while(e == a || e == b || e == c || e == d) e = (e+1)%this->population_size;
	for(int i = 0; i < this->dimension; i++){
		this->temporaly_population[agent_x][i] = this->population[a][i]
												 + this->F*(this->population[b][i]-this->population[c][i])
												 + this->F*(this->population[d][i]-this->population[e][i]);
	}	
}

template<typename T>
void Differential_Evolution<T>::mutation_best_2(int agent_x){

	int a, b, c, d;
	a = rand()%this->population_size;
	b = rand()%this->population_size;
	c = rand()%this->population_size;
	d = rand()%this->population_size;
	while(b == a) b = (b+1)%this->population_size;
	while(c == a || c == b) c = (c+1)%this->population_size;
	while(d == a || d == b || d == c) d = (d+1)%this->population_size;
	for(int i = 0; i < this->dimension; i++){
		this->temporaly_population[agent_x][i] = this->population[this->best_agent_position][i]
												 + this->F*(this->population[a][i]-this->population[b][i])
												 + this->F*(this->population[c][i]-this->population[d][i]);
	}	
}

template<typename T>
void Differential_Evolution<T>::crossover(int agent_x){

	int R = rand()%this->dimension;
	for(int i = 0; i < this->dimension; i++){
		if((((T)rand())/RAND_MAX) < this->CR || i == R);
		else this->temporaly_population[agent_x][i] = this->population[agent_x][i];
	}
	
}

template<typename T>
void Differential_Evolution<T>::selection(int agent_x){
	
	T new_fitness = this->fitness_function(this->temporaly_population[agent_x]);
	T *swap_agent;
	
	if(new_fitness < this->fitness[agent_x]){
		this->fitness[agent_x] = new_fitness;
		if(new_fitness < this->fitness[this->best_agent_position]){
			this->best_agent_position = agent_x;
		}
		swap_agent = this->population[agent_x];
		this->population[agent_x] = this->temporaly_population[agent_x];
		this->temporaly_population[agent_x] = swap_agent;
	}
	
}

template<typename T>
void Differential_Evolution<T>::next_generation(){
	
	switch(this->mutation_type){
		case 0:
			for(int agent_x = 0; agent_x < this->population_size; agent_x++){
				this->mutation_rand_1(agent_x);
			}
		break;
		case 1:
			for(int agent_x = 0; agent_x < this->population_size; agent_x++){
				this->mutation_best_1(agent_x);
			}
		break;
		case 2:
			for(int agent_x = 0; agent_x < this->population_size; agent_x++){
				this->mutation_target_rand(agent_x);
			}
		break;
		case 3:
			for(int agent_x = 0; agent_x < this->population_size; agent_x++){
				this->mutation_target_best(agent_x);
			}
		break;
		case 4:
			for(int agent_x = 0; agent_x < this->population_size; agent_x++){
				this->mutation_rand_2(agent_x);
			}
		break;
		case 5:
			for(int agent_x = 0; agent_x < this->population_size; agent_x++){
				this->mutation_best_2(agent_x);
			}
		break;
		
	}

	for(int agent_x = 0; agent_x < this->population_size; agent_x++){
		if((((T)rand())/RAND_MAX) < this->LS){
			this->local_search_function(this->population[agent_x], this->temporaly_population[agent_x]);
		}else{
			this->crossover(agent_x);
		}
		if(this->agent_repair) this->agent_repair_function(this->temporaly_population[agent_x]);
	}

	this->amount_fitness_evaluation += this->population_size;
	#ifdef _OPENMP
		if(this->parallel){
			#pragma omp parallel for
			for(int i = 0; i < this->population_size; i++){
				this->selection(i);
			}
		}else{
			for(int i = 0; i < this->population_size; i++){
				this->selection(i);
			}	
		}
	#else
		for(int i = 0; i < this->population_size; i++){
			this->selection(i);
		}
	#endif

}

template<typename T>
void Differential_Evolution<T>::evolution(T stop_value, int max_stoped_generation, int max_generation){
	
	clock_t time_base;
	time_base = clock();
	
	int amount_generations_stoped = 0;
	T last_fitness = this->get_mean_fitness();
	printf("(Running): best: %e  mean:%e  generation:0\r", this->get_best_fitness(), this->get_mean_fitness());
	fflush(stdout);
	for(int i = 0 ; i < max_generation && amount_generations_stoped < max_stoped_generation && this->fitness[this->best_agent_position] > stop_value; i++){

		if(((clock()-time_base)/((double)CLOCKS_PER_SEC)) > 2){
			time_base = clock();
			printf("(Running): best: %e  mean:%e  generation:%d\r", this->get_best_fitness(), this->get_mean_fitness(), i);
			fflush(stdout);
		}

		this->next_generation();
		if(this->get_mean_fitness() < last_fitness){
			amount_generations_stoped = 0;
			last_fitness = this->get_mean_fitness();
		}else amount_generations_stoped++;
	}
	printf("                                                                      \r");
	fflush(stdout);
	
}

template<typename T>
void Differential_Evolution<T>::clone_best_agent(T* agent){
	for(int i = 0; i < this->dimension; i++) agent[i] = this->population[this->best_agent_position][i];
}

template<typename T>
void Differential_Evolution<T>::create_agent_basic(T* agent){
	for(int j = 0; j < this->dimension; j++){
		agent[j] = (2*(((T)rand())/RAND_MAX))-1;
	}
}

template<typename T>
void Differential_Evolution<T>::new_population(void (*create_agent_function)(T* agent)){
	
	for(int i = 0; i < this->population_size; i++){
		create_agent_function(this->population[i]);
	}
	this->fitness[0] = this->fitness_function(population[0]);
	this->amount_fitness_evaluation++;
	this->best_agent_position = 0;
	for(int i = 1; i < this->population_size; i++){
		this->fitness[i] = this->fitness_function(population[i]);
		this->amount_fitness_evaluation++;
		if(this->fitness[i] < this->fitness[this->best_agent_position]) this->best_agent_position = i;		
	}
	
}

template<typename T>
void Differential_Evolution<T>::new_population(){
	
	for(int i = 0; i < this->population_size; i++){
		this->create_agent_basic(this->population[i]);
	}
	this->fitness[0] = this->fitness_function(population[0]);
	this->amount_fitness_evaluation++;
	this->best_agent_position = 0;
	for(int i = 1; i < this->population_size; i++){
		this->fitness[i] = this->fitness_function(population[i]);
		this->amount_fitness_evaluation++;
		if(this->fitness[i] < this->fitness[this->best_agent_position]) this->best_agent_position = i;		
	}
	
}

#endif //DIFFERENTIAL_EVOLUTION_TPP

