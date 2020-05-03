#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <sys/time.h>

using namespace std;

// global variables
double* rank_array;
int num_procs = 2;
int num_nodes = 0;
double limit = 0.00001;
int max_iterations = 10;
double damping_factor = 0.85;

MPI_Datatype hybrid_data_type;

struct ultron{
	double* p_ranks;
	int n_msgs,n_pages;
	vector<int> noEntryNodes;
	vector<int> procs_msg_count;
	map<int, vector<int> > map;
};

struct hybrid_data{
	double value;
	int id;
};

void print_graph(ultron ds){
	cout << "Number of messages " << ds.n_msgs << endl;
	cout << "Number of pages " << ds.n_pages << endl;
	cout << "noEntryNodes: ";
	for(int i = 0;i<ds.noEntryNodes.size();i++){
		cout << ds.noEntryNodes[i] << " ";
	}

	cout << "Map: " << endl;
	for(int i = 0;i<ds.n_pages;i++){
		cout << i << "-> ";
		for(int j = 0;j<ds.map[i].size();j++){
			cout << ds.map[i][j] << ", ";
		}
		cout << endl;
	}
	cout << endl;
}

void write_data(double* matrix, int n, string name){
	double sum = 0.0;
	ofstream file;
	file.precision(16);
	file.open("./" + name);
	
	if(file.is_open()){
		for(int i = 0;i<n;i++){
			file << i << " = " << *(matrix + i) << " \n";
			// file << *(matrix + i) << " \n";
			sum+=(*(matrix + i));
		}
		file << "sum " << sum;
	}

	file.close();
}

void initialize2(ultron &ds, string file_path){
	int n = -1;
	int count = 0;
	int a,b,proc_id;

	for(int i = 0;i<=num_procs;i++){
		ds.procs_msg_count.push_back(0);
	}

	ifstream data_file(file_path);
	while(data_file >> a >> b){
		ds.map[a].push_back(b);
		n = max(n,max(a,b));
		count+=1;

		proc_id = (b%num_procs)+1;
		ds.procs_msg_count[proc_id]++;
	}
	n+=1;

	ds.p_ranks = (double*)malloc(n*sizeof(double));
	for(int i = 0;i<n;i++){
		*(ds.p_ranks + i) = 1.0/n;
		if(ds.map.find(i)==ds.map.end()){
			for(int j = 0;j<n;j++){
				// if(j!=i){
					ds.map[i].push_back(j);
					proc_id = (j%num_procs)+1;
					ds.procs_msg_count[proc_id]++;
					count++;
				// }
			}
		}
	}

	ds.n_pages = n;
	ds.n_msgs = count;
	data_file.close();
}

void mapper(ultron ds, MPI_Request reqArray[], MPI_Status statusArray[]){
	int proc_id = 0;
	int count = 0;

	MPI_Request tempReqs[num_procs];
	MPI_Status tempStats[num_procs];
	for(int i = 1;i<=num_procs;i++){
		MPI_Isend(&ds.procs_msg_count[i],1,MPI_INT,i,0,MPI_COMM_WORLD,&tempReqs[i-1]);
	}
	MPI_Waitall(num_procs,tempReqs,tempStats);
	
	// for(int i = 0;i<ds.noEntryNodes.size();i++){
	// 	proc_id = ds.noEntryNodes[i]%num_procs + 1;
	// 	struct hybrid_data send_packet;
	// 	send_packet.id = ds.noEntryNodes[i];
	// 	send_packet.value = 0.0;
	// 	MPI_Isend(&send_packet,1,hybrid_data_type,proc_id,0,MPI_COMM_WORLD,&reqArray[count]);
	// 	count++;
	// }

	for(int i = 0;i<ds.n_pages;i++){
		if(ds.map.find(i)!=ds.map.end()){
			for(int j = 0;j<ds.map[i].size();j++){
				proc_id = (ds.map[i][j]%num_procs) + 1;
				struct hybrid_data send_packet;
				send_packet.id = ds.map[i][j];
				send_packet.value = (*(ds.p_ranks + i))/ds.map[i].size();
				cout << i << "->" << ds.map[i][j] << "->" << send_packet.value << endl;
				MPI_Isend(&send_packet,1,hybrid_data_type,proc_id,0,MPI_COMM_WORLD,&reqArray[count]);
				count++;
			}
		}
	}

	MPI_Waitall(ds.n_msgs,reqArray,statusArray);
}

void reduce(){
	int num_messages,idx,rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Request tempReq1,tempReq2;
	MPI_Status tempStat1, tempStat2;

	MPI_Irecv(&num_messages,1,MPI_INT,0,0,MPI_COMM_WORLD,&tempReq1);
	MPI_Wait(&tempReq1, &tempStat1);
	
	MPI_Request recvArray[num_messages];
	MPI_Status statusArray[num_messages];
	struct hybrid_data datArray[num_messages];

	idx = 0;
	while(idx<num_messages){
		MPI_Irecv(&datArray[idx],1,hybrid_data_type,0,0,MPI_COMM_WORLD,&recvArray[idx]);
		idx++;
	}

	MPI_Waitall(num_messages,recvArray,statusArray);

	for(int i = 0;i<num_nodes;i++){
		*(rank_array + i) = (1.0 - damping_factor)/(num_procs*num_nodes);
	}

	for(int i = 0;i<num_messages;i++){
		idx = datArray[i].id - (rank-1)*num_nodes;
		*(rank_array + idx)+=(damping_factor * datArray[i].value);
	}

	cout.precision(16);
	for(int k = 0;k<num_nodes;k++){
		cout << *(rank_array + k) << ", ";
	}
	cout << endl;
	cout << "reduce waala print done" << endl;

	MPI_Isend(rank_array,num_nodes,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&tempReq2);
	MPI_Wait(&tempReq2,&tempStat2);
}

int main(int argc, char *argv[]){
	string input_file_path = argv[1];
	string output_file_path = argv[2];
	num_procs = atoi(argv[4]) - 1;
	
	ultron ds;
	initialize2(ds,input_file_path);
	num_nodes = ds.n_pages/num_procs;
	
	MPI_Init(&argc,&argv);

   	int lengths[2] = {1,1};
   	const MPI_Aint displacement[2] = {0, sizeof(double)};
   	MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
   	MPI_Type_create_struct(2,lengths,displacement,types,&hybrid_data_type);
   	MPI_Type_commit(&hybrid_data_type);

   	int rank,i;
   	double curr_diff = 1.0;
   	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   	rank_array = (double*)malloc(ds.n_pages*sizeof(double));

   	while(i<max_iterations && curr_diff>limit){
   		if(rank==0){
   			// print_graph(ds);
   			MPI_Request reqArray[ds.n_msgs];
   			MPI_Status statusArray[ds.n_msgs];

   			cout.precision(16);
   			for(int k = 0;k<ds.n_pages;k++){
   				cout << *(ds.p_ranks + k) << ", ";
   			}
   			cout << "main function waala print done" << endl;
   			mapper(ds, reqArray, statusArray);

   			MPI_Request reqArray2[num_procs];
   			MPI_Status statArray2[num_procs];
   			
   			for(int proc_id = 1;proc_id<=num_procs;proc_id++){
   				MPI_Irecv(rank_array + num_nodes*(proc_id-1),num_nodes,MPI_DOUBLE,proc_id,0,MPI_COMM_WORLD,&reqArray2[proc_id-1]);
   			}
   			MPI_Waitall(num_procs,reqArray2,statArray2);

   			curr_diff = 0.0;
   			double sum = 0.0;
   			for(int index = 0;index<ds.n_pages;index++){
   				curr_diff+=fabs(*(ds.p_ranks+index) - *(rank_array+index));
   				*(ds.p_ranks + index) = (*(rank_array + index));
   				sum+=(*(rank_array + index));
   			}

   			for(int index = 0;index<ds.n_pages;index++){
   				*(ds.p_ranks + index) = (*(ds.p_ranks + index))/sum;
   			}
   			
   			cout << i << "->" << curr_diff << endl;
   			MPI_Request reqArray3[num_procs];
   			MPI_Status statArray3[num_procs];
   			for(int proc_id = 1;proc_id<=num_procs;proc_id++){
   				MPI_Isend(&curr_diff,1,MPI_DOUBLE,proc_id,0,MPI_COMM_WORLD,&reqArray3[proc_id-1]);
   			}

   			MPI_Waitall(num_procs,reqArray3,statArray3);
   			
   			if(curr_diff<limit || i==max_iterations-1){
   				write_data(ds.p_ranks,ds.n_pages,output_file_path);
   			}
   		}else{
   			reduce();
   			MPI_Request tempReq3;
   			MPI_Status tempStat3;
   			MPI_Irecv(&curr_diff,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&tempReq3);
   			MPI_Wait(&tempReq3, &tempStat3);
   		}
   		i++;
   	}

	MPI_Finalize();
	return 0;
}