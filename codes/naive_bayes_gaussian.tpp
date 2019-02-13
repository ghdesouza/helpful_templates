/**
 * @file naive_bayes_gaussian.tpp
 *
 * @brief Naive Bayes Gaussian Classifier
 *
 * @details This file implement a Naive Bayes Gaussian Classifier
 *
 * @author Gabriel Henrique de Souza (ghdesouza@gmail.com)
 *
 * @date february 13, 2019
 *
 * @copyright Distributed under the Mozilla Public License 2.0 ( https://opensource.org/licenses/MPL-2.0 )
 *
 * @see https://github.com/ghdesouza/helpful_templates
 *
 * Created on: february 13, 2019
 *
 */

#ifndef NAIVE_BAYES_GAUSSIAN_TPP
#define NAIVE_BAYES_GAUSSIAN_TPP

#include <math.h>
#include "classifiers.tpp"

using namespace std;

template<typename T = float>
class Naive_Bayes_Gaussian : public Classifiers<T>{
	
	protected:
	
		T* py;
		T** means;
		T** variances;
		T* fitness;

		T calculate_px_y(int posicao, T valor, int classe); // P( X_j | y )
		T calculate_pX_y(T *X, int y); // P( X | y )
		void calculate_py_X(T *X); // P( y | X )

	public:
		Naive_Bayes_Gaussian(int dimension, int amount_labels);
		~Naive_Bayes_Gaussian();
		
		void fit(Dataset<T>* data);
		int predict(T *X);
		T get_fitness(int label){}
};

template<typename T>
Naive_Bayes_Gaussian<T>::Naive_Bayes_Gaussian(int dimension, int amount_labels){
	
	this->dimension = dimension;
	this->amount_labels = amount_labels;
	
	this->py = new T[this->amount_labels];
	this->fitness = new T[this->amount_labels];
	for(int i = 0; i < this->amount_labels; i++) this->fitness[i] = 0.0;
	
	this->means = new T*[this->amount_labels];
	this->variances = new T*[this->amount_labels];
	for(int i = 0; i < this->amount_labels; i++){
		this->means[i] = new T[this->dimension];
		this->variances[i] = new T[this->dimension];
	}

}

template<typename T>
Naive_Bayes_Gaussian<T>::~Naive_Bayes_Gaussian(){
	
	delete[] this->py;
	delete[] this->fitness;
	for(int i = 0; i < this->amount_labels; i++){
		delete[] this->means[i];
		delete[] this->variances[i];
	}	delete[] this->means;
		delete[] this->variances;
}

template<typename T>
T Naive_Bayes_Gaussian<T>::calculate_px_y(int position, T value, int label){
	
	// gaussian kernel
	return ((1.0/(sqrt(2.0*M_PI*this->variances[label-1][position])))*
	exp(-pow(value-this->means[label-1][position], 2)/(2*this->variances[label-1][position])));
}

template<typename T>
T Naive_Bayes_Gaussian<T>::calculate_pX_y(T *X, int y){
	
	T pX_y = 1;

	for(int i = 0; i < this->dimension; i++){
		pX_y *= this->calculate_px_y(i, X[i], y);
	}
	
	return pX_y;
}

template<typename T>
void Naive_Bayes_Gaussian<T>::calculate_py_X(T *X){

	T sum = 0;
	
	for(int i = 0; i < this->amount_labels; i++){
		this->fitness[i] = py[i]*this->calculate_pX_y(X, i+1);
		sum += this->fitness[i];
	}
	for(int i = 0; i < this->amount_labels; i++) this->fitness[i] /= sum;	
	
	return;
}

template<typename T>
void Naive_Bayes_Gaussian<T>::fit(Dataset<T>* data){

	int amount_trials = data->get_train_size();
	T* temp = new T[amount_trials];
	int amount_trials_of_label;
	
	for(int i = 0; i < this->dimension; i++){
		for(int k = 0; k < this->amount_labels; k++){
			amount_trials_of_label = 0;
			for(int j = 0; j < amount_trials; j++){
				if(k+1 == data->get_label_train(j)){
					temp[amount_trials_of_label] = data->get_trial_train(j)[i];
					amount_trials_of_label++;
				}
			}
			this->means[k][i] = mean(temp, amount_trials_of_label);
			this->variances[k][i] = variance(temp, amount_trials_of_label);
			if(this->variances[k][i] < 1e-8){
				this->variances[k][i] = RAND_MAX;
				printf("\n\nERROR: variance of %dth dimension of %dth label is vary small!\n\n", i, k);
			}
		}
	}
	
	delete[] temp;

	for(int j = 0; j < this->amount_labels; j++) this->py[j] = 0;	
	for(int i = 0; i < amount_trials; i++){
		this->py[data->get_label_train(i)-1]++;
	}
	for(int j = 0; j < this->amount_labels; j++){
		this->py[j] /= amount_trials;
	}
		
}

template<typename T>
int Naive_Bayes_Gaussian<T>::predict(T *X){
	
	this->calculate_py_X(X);
	
	int label = 0;
	T prob = this->fitness[0];
	for(int i = 1; i < this->amount_labels; i++){
		if(this->fitness[i] > prob){
			label = i;
			prob = this->fitness[i];
		}
	}
	
	return (label+1);
}


#endif // NAIVE_BAYES_GAUSSIAN_TPP

