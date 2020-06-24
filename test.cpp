#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <sys/time.h>

#include "mapreduce.h"
#include "keyvalue.h"
using namespace MAPREDUCE_NS;

struct ultron{
	int n_pages;
	double* p_ranks;
	map<int, vector<int> > map;
}

// GLOBAL VARIABLES
double limit = 0.00001;
int max_iterations = 1000;
double damping_factor = 0.85;

void write_data(double* matrix, int n, string name){
	double sum = 0.0;
	ofstream file;
	file.precision(16);
	file.open("./" + name);
	
	if(file.is_open()){
		for(int i = 0;i<n;i++){
			file << i << " = " << *(matrix + i) << " \n";
			sum+=(*(matrix + i));
		}
		file << "sum " << sum;
	}

	file.close();
}

void initialize(ultron &ds, string file_path){
	int n = -1;
	int a,b;

	ifstream data_file(file_path);
	while(data_file >> a >> b){
		ds.map[a].push_back(b);
		n = max(n,max(a,b));
	}
	n+=1;

	ds.p_ranks = (double*)malloc(n*sizeof(double));
	for(int i = 0;i<n;i++){
		*(ds.p_ranks + i) = 1.0/n;
	}

	ds.n_pages = n;
	data_file.close();
}

int main(int argc, char const *argv[]){
	return 0;
}