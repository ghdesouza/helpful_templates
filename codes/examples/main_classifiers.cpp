#include "../naive_bayes_gaussian.tpp"
#include "../naive_bayes_parzen_window.tpp"

#include <time.h>

int main(){

	srand(time(NULL));

	int data_size = 150;
	int trial_size = 3;
	int amount_labels = 4;
	int kfold_size = 10;
	int amount_kfolds = 10;
	float* accuraces = new float[kfold_size*amount_kfolds];
	Classifiers<>* classifier;
	Dataset<>* dset;
	
	//classifier = new Naive_Bayes_Gaussian<>(trial_size, amount_labels);
	classifier = new Naive_Bayes_Parzen_Window<>(trial_size, amount_labels);
	
	{ // dataset load (this dataset was created based in normal distributions)
		dset = new Dataset<>("datas/random_dataset.txt", data_size, trial_size, kfold_size);
	}

	{ // k-fold cross-validation
		for(int j = 0; j < amount_kfolds; j++){
			dset->shuffle_data(); // shuffle
			for(int i = 0; i < kfold_size; i++){
				dset->set_test_kfold_id(i); // next fold
				classifier->fit(dset); // training
				accuraces[j*kfold_size+i] = classifier->accurace(dset);
			}
		}
		printf("(%f +- %f) %%\n", 100*mean(accuraces, kfold_size*amount_kfolds), 100*sqrt(variance(accuraces, kfold_size*amount_kfolds)));
	}

	delete dset;
	delete classifier;
	delete[] accuraces;
	return 0;
}

