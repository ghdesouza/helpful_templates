#include "../differential_evolution.tpp"

#include <stdio.h>
#include <time.h>
#include <iostream>
using namespace std;

int dim = 25;
double cr = 0.9;
double f = 0.7;
int NP = 1000;
double bl = 0.001;

Differential_Evolution<double>* de;

double fitness_function(double* agent);
void repair_function(double* agent);
void local_search_function(double* agent, double* new_agent);
void create_agent_function(double* agent);

int main(){

	srand(time(NULL));
	de = new Differential_Evolution<double>(dim, NP, fitness_function, f, cr);
	de->set_agent_repair(repair_function);
	de->set_local_search(local_search_function, bl);
	de->clean_amount_fitness_evaluation();
	de->new_population(create_agent_function);
	de->set_parallel(true);

	de->evolution(1e-5, 1e+1, 1e+9);
	
	printf("Best Fitness: %.10f\n", de->get_best_fitness());
	printf("Number of Evaluation: %d\n", de->get_amount_fitness_evaluation());
	
	return 0;
}

double fitness_function(double* agent){	
	// ACKLEY FUNCTION
	
	double a = 20.0;
	double b = 0.2;
	double c = 2.0*M_PI;
	
	double sum_1 = 0;
	double sum_2 = 0;
	
	for(int i = 0; i < dim; i++) sum_1 += agent[i]*agent[i];
	sum_1 = -a*exp(-b*sqrt(sum_1/dim));
	
	for(int i = 0; i < dim; i++) sum_2 += cos(c*agent[i]);
	sum_2 = -exp(sum_2/dim);
	
	return sum_1+sum_2+a+exp(1);
}

void repair_function(double* agent){
	for(int i = 0; i < dim; i++){
		if(agent[i] > 30) agent[i] = 30;
		else if(agent[i] < -30) agent[i] = -30;
	}
}

void local_search_function(double* agent, double* new_agent){ // in this case the local_search_function is used to create a new mutation
	double *temp_agent = new double[dim];
	de->clone_best_agent(temp_agent);
	for(int i = 0; i < dim; i++) new_agent[i] = agent[i]+1.5*(((double)rand())/RAND_MAX)*(temp_agent[i]-agent[i]);
	delete[] temp_agent;
}

void create_agent_function(double* agent){
	
	for(int j = 0; j < dim; j++) agent[j] = ((2.0*(((double)rand())/RAND_MAX))-1.0)*30.0;
}

