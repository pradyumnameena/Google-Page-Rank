#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sys/time.h>

using namespace std;

vector<double*> initialize(int n){
	vector<double*> matrix;

	for(int i = 0;i<n;i++){
		double* address = (double*) malloc(sizeof(double));
		matrix.push_back(address);
		*(matrix[i]) = 1.0/n;
	}

	return matrix;
}

vector<double*> read_data(string file_path,double fac1){
	int a,b;
	int n = -1;
	double fac2 = 1.0 - fac1;
	vector<double*> matrix;
	map<int, vector<int> > map;
	
	ifstream data_file(file_path);

	while(data_file >> a >> b){
		map[a].push_back(b);
		n = max(n,max(a,b));
	}

	n+=1;
	data_file.close();

	for(int i = 0;i<n;i++){
		double* address = (double*) malloc(n*sizeof(double));
		matrix.push_back(address);
	}
	
	for(int i = 0;i<n;i++){
		for(int j = 0;j<n;j++){
			*(matrix[i] + j) = 0;
		}
	}

	for(int i = 0;i<n;i++){
		if(map.find(i)!=map.end()){
			for(int j = 0;j<map[i].size();j++){
				*(matrix[map[i][j]] + i) += (fac1/map[i].size());
			}
		}else{
			for(int j = 0;j<n;j++){
				*(matrix[j] + i) += (fac1/n);	
			}
		}
	}

	for(int i = 0;i<n;i++){
		for(int j = 0;j<n;j++){
			*(matrix[i]+j) += fac2/n;
		}
	}

	return matrix;
}

void print(vector<double*> matrix, int m, int n){
	cout.precision(16);
	double sum = 0;
	for(int i = 0;i<m;i++){
		for(int j = 0;j<n;j++){
			cout << *(matrix[i] + j) << ", ";
			sum+=(*(matrix[i] + j));
		}
	}
	cout << endl;
	cout << "****************       THE END       ****************" << endl;
}

void write_data(vector<double*> matrix,string name){
	int n = matrix.size();
	double sum = 0;
	ofstream file;
	file.precision(16);
	file.open("./" + name);
	
	if(file.is_open()){
		for(int i = 0;i<n;i++){
			file << i << " = " << *(matrix[i]) << " \n";
			sum+=(*(matrix[i]));
		}
		file << "sum " << sum;
	}

	file.close();
}

double calcDiff(vector<double*> m1,vector<double*> m2){
	int n = m1.size();
	double ans = 0;

	for(int i = 0;i<n;i++){
		ans+=fabs(*(m1[i]) - *(m2[i]));
		*(m2[i]) = *(m1[i]);
	}

	return ans;
}

void matrixMult(vector<double*> m1, vector<double*> m2,vector<double*> m3){
	int n = m1.size();
	double val = 0;
	for(int i = 0;i<n;i++){
		val = 0;
		for(int j = 0;j<n;j++){
			val+=(*(m1[i] + j) * (*(m2[j])));
		}
		*(m3[i]) = val;
	}
}

void pagerank(vector<double*> googleM, vector<double*> initM, double limit,string name){
	vector<double*> res;
	double diff = 1;
	int n = initM.size();
	int iteration_count = 0;

	for(int i = 0;i<n;i++){
		double* address = (double*) malloc(sizeof(double));
		res.push_back(address);
	}

	while(diff>limit && iteration_count<200){
		matrixMult(googleM,initM,res);
		diff = calcDiff(initM,res);
		cout << iteration_count << " iteration : " << diff << endl;
		iteration_count++;
	}

	cout << "Iteration Count- " << iteration_count << endl;
}

int main(int argc, char const *argv[]){
	string input_file_path = argv[1];
	string output_file_path = argv[2];
	double dampening_fac = 0.85;
	double limit = 0.00001;

	vector<double*> googleM = read_data(input_file_path, dampening_fac);
	vector<double*> initM = initialize(googleM.size());
	
	struct timeval start, end;
   	double time_taken;
   
   	gettimeofday(&start, NULL);
	pagerank(googleM,initM,limit,output_file_path);
	gettimeofday(&end, NULL);

	time_taken = (end.tv_sec - start.tv_sec) * 1e6;
	time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
	cout << "Time taken is " << time_taken << endl;

	write_data(initM,output_file_path);
	return 0;
}