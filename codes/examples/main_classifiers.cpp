#include "../naive_bayes_gaussian.tpp"
#include "../naive_bayes_parzen_window.tpp"
#include "../dataset.tpp"

#include <time.h>

const int kfold_size = 10;
const int amount_kfolds = 1;

float **data_full, **data_train, **data_test;
int *labels_full, *labels_train, *labels_test;
const int data_full_size = 150;
const int block_size = (int)(data_full_size/kfold_size);
const int data_train_size = data_full_size-block_size;
const int data_test_size = block_size;

const int trial_size = 3;
const int amount_labels = 4;

int main(){

	srand(time(NULL));
	FILE* file;
	string name_file = "datas/random_dataset.txt";
    Classifiers<>* classifier;
	
	{ // creating datasets
			labels_full = new int[data_full_size];
			data_full = new float*[data_full_size];	
			for(int i = 0; i < data_full_size; i++) data_full[i] = new float[trial_size];
			data_train = new float*[data_train_size]; labels_train = new int[data_train_size];
			data_test = new float*[data_test_size];   labels_test = new int[data_test_size];
	}
	
	classifier = new Naive_Bayes_Gaussian<>(trial_size, amount_labels);
	//classifier = new Naive_Bayes_Parzen_Window<>(trial_size, amount_labels);
    
    
	{ // k-fold cross-validation
		load_dataset(data_full, labels_full, trial_size, data_full_size, name_file.c_str()); // loading dataset
		for(int j = 0; j < amount_kfolds; j++){
			shuffle_data(data_full, labels_full, data_full_size);
			for(int i = 0; i < kfold_size; i++){ 
                printf("\n%d %d\n", j, i);
				setfold_data(data_full, data_train, data_test, 
							 labels_full, labels_train, labels_test, 
							 data_full_size, kfold_size, i); // next fold
				classifier->fit(data_train, labels_train, data_train_size); // training
                printf("%.2f\n", classifier->accurace(data_test, labels_test, data_test_size));
			}
		}
	}
	
	{ // deleting dataset
		delete[] labels_full;
		for(int i = 0; i < data_full_size; i++) delete[] data_full[i]; 
		delete[] data_full;
		delete[] data_train; delete[] labels_train;
		delete[] data_test; delete[] labels_test;	
	}
    
	delete classifier;
	return 0;
}

