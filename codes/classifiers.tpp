/**
 * @file dataset.tpp
 *
 * @brief Abstract class for general classifiers
 *
 * @details This file implement a abstract class with avaliations approach.
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

#ifndef CLASSIFIERS_TPP
#define CLASSIFIERS_TPP

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "basic_stats.tpp"
#include "dataset.tpp"

using namespace std;

template<typename T = float>
class Classifiers{
	
    protected:
		int dimension;
		int amount_labels;
		
    public:
		virtual ~Classifiers(){};
		
		virtual int get_dimension(){return this->dimension;}
		virtual int get_amount_labels(){return this->amount_labels;}
		
		virtual void fit(Dataset<T>* data) = 0;
		virtual int predict(T *X) = 0;
		virtual T get_fitness(int label) = 0;
		
		virtual float accurace(Dataset<T>* data){
			float correct = 0.0;
			int amount_trials = data->get_test_size();
			for(int i = 0; i < amount_trials; i++){
				if(this->predict(data->get_trial_test(i)) == data->get_label_test(i)){
					correct += 1;
				}
			}
			return correct/amount_trials;
		}
		
		virtual float kappa(Dataset<T>* data){
			int** confusion_matrix = new int*[this->amount_labels];
			for(int i = 0; i < this->amount_labels; i++){
				confusion_matrix[i] = new int[this->amount_labels];
				for(int j = 0; j < this->amount_labels; j++){
					confusion_matrix[i][j] = 0;
				}
			}
			int amount_trials = data->get_test_size();
			for(int i = 0; i < amount_trials; i++){
				confusion_matrix[data->get_label_test(i)-1][this->predict(data->get_trial_test(i))-1]++;
			}
			
			float p_0, p_e = 0.0, total = 0.0, correct = 0.0;
			float *expected = new float[this->amount_labels];
			float *finded = new float[this->amount_labels];
			for(int i = 0; i < this->amount_labels; i++){
				expected[i] = 0;
				finded[i] = 0;
				correct += confusion_matrix[i][i];
				for(int j = 0; j < this->amount_labels; j++){
					expected[i] += confusion_matrix[j][i];
					finded[i] += confusion_matrix[i][j];
					total += confusion_matrix[i][j];
				}
			}

			p_0 = correct/total;
			for(int i = 0; i < this->amount_labels; i++) p_e += (finded[i]/total)*(expected[i]/total);

			delete[] expected;
			delete[] finded;
			
			for(int i = 0; i < this->amount_labels; i++){
				delete[] confusion_matrix[i];
			}	delete[] confusion_matrix;
			
			return (p_0-p_e)/(1.0-p_e);
			
		}
		
		virtual float cross_entropy(Dataset<T>* data){
			T crossentropy = 0;
			int temp_class;

			for(int i = 0; i < data->get_test_size(); i++){
				temp_class = this->predict(data->get_trial_test(i));
				for(int c = 1; c <= this->amount_labels; c++){
					if(c == data->get_label_test(i)) crossentropy -= log(this->get_fitness(c));
					else crossentropy -= log(1-this->get_fitness(c));
				}
			}
			return crossentropy/(data->get_test_size());
		}

		virtual float mean_squared_error(Dataset<T>* data){
			T error = 0;
			int temp_class;

			for(int i = 0; i < data->get_test_size(); i++){
				temp_class = this->predict(data->get_trial_test(i));
				error += pow((1-this->get_fitness(data->get_label_test(i))), 2);
			}
			return error/data->get_test_size();
		}

};

#endif // CLASSIFIERS_TPP

