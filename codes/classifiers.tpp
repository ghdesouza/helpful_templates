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

};

#endif // CLASSIFIERS_TPP

