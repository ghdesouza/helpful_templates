/**
 * @file dataset.tpp
 *
 * @brief Struct to save datasets for classifiers
 *
 * @details This file implement a Struct that is used to training and aviliable classifiers
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

#ifndef DATASET_TPP
#define DATASET_TPP

#include <fstream>

using namespace std;

template<typename T=float>
class Dataset{

	private:
	
		int data_size;
		int trial_size;
		T** data;
		int* labels;
		
		int kfold_size;
		int test_kfold_id;
		int train_size;
		int test_size;
		
		int start_test_position;
		
		void load_dataset(string name_file);

	public:
	
		Dataset(T** data, int* labels, int data_size, int trial_size, int kfold_size);
		Dataset(string name_file, int data_size, int trial_size, int kfold_size);
		~Dataset();
		
		void set_test_kfold_id(int number);
		void set_kfold_size(int kfold_size);
		
		int get_data_size(){ return this->data_size; }
		int get_trial_size(){ return this->trial_size; }
		int get_kfold_size(){ return this->kfold_size; }
		int get_test_kfold_id(){ return this->test_kfold_id; }
		int get_train_size(){ return this->train_size; }
		int get_test_size(){ return this->test_size; }
		
		T* get_trial_train(int number);
		T* get_trial_test(int number);
		
		int get_label_train(int number);
		int get_label_test(int number);
		
		void shuffle_data();

};

template<typename T>
Dataset<T>::Dataset(T** data, int* labels, int data_size, int trial_size, int kfold_size){
	this->data_size = data_size;
	this->trial_size = trial_size;
	
	this->data = new T*[this->data_size];
	this->labels = new int[this->data_size];
	for(int i = 0; i < this->data_size; i++){
		this->data[i] = new T[this->trial_size];
		this->labels[i] = labels[i];
		for(int j = 0; j < this->trial_size; j++){
			this->data[i][j] = data[i][j];
		}
	}
	
	this->set_kfold_size(kfold_size);
}

template<typename T>
Dataset<T>::Dataset(string name_file, int data_size, int trial_size, int kfold_size){
	this->data_size = data_size;
	this->trial_size = trial_size;
	
	this->load_dataset(name_file);
	
	this->set_kfold_size(kfold_size);
}

template<typename T>
Dataset<T>::~Dataset(){
	for(int i = 0; i < this->data_size; i++){
		delete[] this->data[i];
	}	delete[] this->data;
		delete[] this->labels;
}

template<typename T>
void Dataset<T>::load_dataset(string name_file){

	FILE* arq = fopen(name_file.c_str(), "r");
    if(!arq){
        printf("ERROR (DATASET): File not found!\n");
        return;
    }
	this->data = new T*[this->data_size];
	this->labels = new int[this->data_size];
	for(int i = 0; i < this->data_size; i++){
		this->data[i] = new T[this->trial_size];
		for(int j = 0; j < this->trial_size; j++){
			fscanf(arq, "%f", &this->data[i][j]);
		}
		fscanf(arq, "%d", &this->labels[i]);
	}
	fclose(arq);
}

template<typename T>
void Dataset<T>::set_test_kfold_id(int number){
	if(this->kfold_size == 1){
		this->train_size = this->data_size;
		this->test_size = 0;
		this->start_test_position = -1;
		return;
	}
	
	if(number >= 0 && number < this->kfold_size){
		this->test_kfold_id = number;
		this->test_size = (int) (this->data_size/this->kfold_size);
		this->start_test_position = this->test_kfold_id*this->test_size;
		if(this->test_kfold_id < (this->data_size % this->kfold_size)){
			this->test_size += 1;
			this->start_test_position += this->test_kfold_id;
		}else{
			this->start_test_position += (this->data_size % this->kfold_size);
		}
		this->train_size = this->data_size-this->test_size;		
	}else{
		printf("ERROR (DATASET CLASS): invalid test_kfold_id!\n");
	}
}

template<typename T>
T* Dataset<T>::get_trial_test(int number){
	if(number >= 0 && number < this->test_size){
		return this->data[this->start_test_position+number];
	}else{
		printf("ERROR (DATASET CLASS): invalid id for trial on test!\n");
		return this->data[0];
	}
}

template<typename T>
T* Dataset<T>::get_trial_train(int number){
	if(number >= 0 && number < this->train_size){
		if(number < this->start_test_position){
			return this->data[number];
		}else{
			return this->data[this->test_size+number];	
		}
	}else{
		printf("ERROR (DATASET CLASS): invalid id for trial on train!\n");
		return this->data[0];
	}
}

template<typename T>
int Dataset<T>::get_label_test(int number){
	if(number >= 0 && number < this->test_size){
		return this->labels[this->start_test_position+number];
	}else{
		printf("ERROR (DATASET CLASS): invalid id for label on test!\n");
		return this->labels[0];
	}
}

template<typename T>
int Dataset<T>::get_label_train(int number){
	if(number >= 0 && number < this->train_size){
		if(number < this->start_test_position){
			return this->labels[number];
		}else{
			return this->labels[this->test_size+number];	
		}
	}else{
		printf("ERROR (DATASET CLASS): invalid id for label on train!\n");
		return this->labels[0];
	}
}

template<typename T>
void Dataset<T>::set_kfold_size(int kfold_size){
	
	if(kfold_size < 1 || kfold_size > this->data_size) kfold_size = 1;
	this->kfold_size = kfold_size;
	this->set_test_kfold_id(0);
	this->shuffle_data();
}

template<typename T>
void Dataset<T>::shuffle_data(){
	
	int* order = new int[this->data_size];
	int pos_max;
	int max_val;
	T* temp_T;
	int temp_int; // for label and pivot
		
	for(int i = 0; i < this->data_size; i++) order[i] = (int)rand();

	for(int i = 0; i < this->data_size; i++){
		pos_max = i;
		max_val = order[i];
		for(int j = i+1; j < this->data_size; j++){
			if(order[j] > max_val){
				pos_max = j;
				max_val = order[j];
			}
		}
		
		{ // swaps
			// trial
			temp_T = data[i];
			data[i] = data[pos_max];
			data[pos_max] = temp_T;
			// label
			temp_int = labels[i];
			labels[i] = labels[pos_max];
			labels[pos_max] = temp_int;
			// randon vector
			temp_int = order[i];
			order[i] = order[pos_max];
			order[pos_max] = temp_int;
		}
	}
	
	delete[] order;
}

#endif // DATASET_TPP

