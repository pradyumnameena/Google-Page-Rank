#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sys/time.h>

using namespace std;

struct ultron{
	int n_pages;
	double* p_ranks;
	map<int, vector<int> > map;
};

double limit = 0.00001;
int max_iterations = 1000;
double dampening_fac = 0.85;

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

void pagerank(ultron ds, string output_file_path){
	double* res = (double*)malloc(ds.n_pages * sizeof(double));
	double dangling_contri = 0.0;
	int iteration_count = 0;
	double diff = 10.0;
	int n = ds.n_pages;
	double val,sum;

	struct timeval start, end;
   	double time_taken;
   	gettimeofday(&start, NULL);
   	
	while(diff>limit && iteration_count<max_iterations){
		sum = 0.0;
		diff = 0.0;
		dangling_contri = 0.0;

		for(int i = 0;i<n;i++){
			*(res+i) = (1.0 - dampening_fac)/n;
		}

		for(int i = 0;i<n;i++){
			if(ds.map.find(i)!=ds.map.end()){
				val = (*(ds.p_ranks + i))/ds.map[i].size();
				for(int j = 0;j<ds.map[i].size();j++){
					*(res + ds.map[i][j]) += (dampening_fac * val);
				}
			}else{
				dangling_contri+=(*(ds.p_ranks + i));
			}
		}

		dangling_contri/=n;
		dangling_contri*=dampening_fac;
		
		for(int i = 0;i<n;i++){
			*(res + i) += dangling_contri;
			diff+=fabs(*(res+i) - *(ds.p_ranks + i));
			*(ds.p_ranks + i) = *(res + i);
			sum+=(*(ds.p_ranks + i));
		}

		cout << iteration_count << "->" << diff << endl;

		iteration_count++;
	}

	gettimeofday(&end, NULL);
	time_taken = (end.tv_sec - start.tv_sec) * 1e6;
	time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

	write_data(ds.p_ranks,n,output_file_path);
	cout << "Iteration Count- " << iteration_count << endl;
	cout << "Time taken " << time_taken << " sec" << endl; 
}

int main(int argc, char const *argv[]){
	string input_file_path = argv[1];
	string output_file_path = argv[2];
	max_iterations = atoi(argv[3]);
	
	ultron ds;
	initialize(ds,input_file_path);
	pagerank(ds,output_file_path);
	return 0;
}