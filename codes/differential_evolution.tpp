/**
 * @file differential_evolution.tpp
 *
 * @brief Differential Evolution Algorithm
 *
 * @details This file implement the Differential Evolution Algorithm for generic use includes constraint and local search problems.
 *
 * @author Gabriel Henrique de Souza (ghdesouza@gmail.com)
 *
 * @date january 15, 2019
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
		* @brief Function that create a new agent based in DE moviment or Local Search
		* @details Move the agent based in other 3 differentes agents using the equation:
		* Y = (R) x (A + F x (B-C)) + (!R) x X where R is a boolenan vector created randomicaly with the CR probability.
		* The new agent will saved in a private vector of class "temporaly_population" in same position that agent_x.
		* If the Local Search probability is bigger then zero so this case is avaliable before move the agent.
		* @param agent_x - Agent position that will moved
		* @param agent_a - First agent position used for move the agent_x
		* @param agent_b - Second agent position used for move the agent_x
		* @param agent_c - Third agent position used for move the agent_x
		* @return none
		*/
		void create_movement(int agent_x, int agent_a, int agent_b, int agent_c);
		
		/**
		* @brief Try to move the agent of temporaly population to population
		* @details Verify if the fitness function is improved in temporaly population so move the agent for population
		* and if is the new best fitness update the best_agent_position.
		* @param agent_x - Agent position that will moved
		* @return none
		*/
		void try_move(int agent_x);
		
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
Differential_Evolution<T>::~Differential_Evolution(){

	for(int i = 0; i < this->population_size; i++){
		delete[] this->population[i] = new T[this->dimension];
		delete[] this->temporaly_population[i] = new T[this->dimension];
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
void Differential_Evolution<T>::create_movement(int agent_x, int agent_a, int agent_b, int agent_c){
	
	if((((T)rand())/RAND_MAX) < this->LS){

		this->local_search_function(this->population[agent_x], this->temporaly_population[agent_x]);

	}else{
			
		int R = rand()%this->dimension;
		for(int i = 0; i < this->dimension; i++){
			if((((T)rand())/RAND_MAX) < CR || i == R){
				this->temporaly_population[agent_x][i] = this->population[agent_a][i] + F*(this->population[agent_b][i]-this->population[agent_c][i]);
			}else{
				this->temporaly_population[agent_x][i] = this->population[agent_x][i];
			}
		}
		
	}
	
	if(this->agent_repair) this->agent_repair_function(this->temporaly_population[agent_x]);
	
}

template<typename T>
void Differential_Evolution<T>::try_move(int agent_x){
	
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
	
	int a, b, c;
	for(int i = 0; i < this->population_size; i++){
		a = rand()%this->population_size;
		b = rand()%this->population_size;
		c = rand()%this->population_size;
		while(a == i) a = rand()%this->population_size;
		while(b == i || b == a) b = rand()%this->population_size;
		while(c == i || c == a || c == b) c = rand()%this->population_size;
		this->create_movement(i, a, b, c);
	}
	this->amount_fitness_evaluation += this->population_size;
	#ifdef _OPENMP
		if(this->parallel){
			#pragma omp parallel for
			for(int i = 0; i < this->population_size; i++){
				this->try_move(i);
			}
		}else{
			for(int i = 0; i < this->population_size; i++){
				this->try_move(i);
			}	
		}
	#else
		for(int i = 0; i < this->population_size; i++){
			this->try_move(i);
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

		if(((clock()-time_base)/((double)CLOCKS_PER_SEC)) > 4){
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

