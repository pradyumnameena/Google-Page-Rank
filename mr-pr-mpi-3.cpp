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

MPI_Datatype data_packet_mpi;

// global variables
int num_procs = 2;
int num_nodes = 0;
double limit = 0.00001;
int max_iterations = 1000;
double damping_factor = 0.85;

struct data_packet{
	double value;
	int id;
};

struct graph_container_ds{
	double* p_ranks;
	int n_msgs,n_pages;
	vector<int> procs_msg_count;
	map<int, vector<int> > map;
};

int num_messages;
double* rank_array;
double dangling_contri;
data_packet* datArray2;

void write_data(double* matrix, int n, string name){
	int i;
	double sum;
	ofstream file;
	file.open("./" + name);
	file.precision(16);
	
	if(file.is_open()){
		sum = 0.0;
		while(i<n){
			file << i;
			file << " = ";
			file << *(matrix + i);
			file << " \n";
			i++;
			sum+=(*(matrix + i));
		}
		file << "sum ";
		file << sum;
	}

	file.close();
}

void initialize(graph_container_ds &ds, string file_path){
	int i = 0;
	int n = -1;
	int count = 0;
	int from = 0;
	int to = 0;
	int temp;
	int processor_id = 0;

	while(i<num_procs){
		ds.procs_msg_count.push_back(0);
		i++;
	}

	ifstream data_file(file_path);
	while(data_file >> from >> to){
		temp = max(from,to);
		n = max(temp,n);
		ds.map[from].push_back(to);
	}
	num_nodes = (n+1)/num_procs;
	n+=1;

	ds.p_ranks = (double*)malloc(n*sizeof(double));
	
	i = 0;
	while(i<n){
		*(ds.p_ranks + i) = 1.0/n;
		int j = 0;
		if(ds.map.find(i)!=ds.map.end()){
			while(j<ds.map[i].size()){
				processor_id = ds.map[i][j];
				processor_id/=num_nodes;

				ds.procs_msg_count[processor_id]++;

				if(processor_id>0){
					count++;
				}

				j++;
			}
		}
		i++;
	}

	ds.n_msgs = count + (num_procs-1);
	data_file.close();
	ds.n_pages = n;
}

void reducer_function_2(){
	int i,idx;

	i = 0;
	while(i<num_nodes){
		*(rank_array + i) = (dangling_contri*damping_factor);
		*(rank_array + i)+=(1.0 - damping_factor)/(num_procs*num_nodes);
		i++;
	}

	i = 0;
	struct data_packet recv_packet;
	while(i<num_messages){
		recv_packet = *(datArray2 + i);
		*(rank_array + recv_packet.id)+=(damping_factor * recv_packet.value);
		i++;
	}
}

void map_function(graph_container_ds ds, MPI_Request reqArray[], MPI_Status statusArray[]){
	int i,idx,count,proc_id;

	MPI_Request tempReqs[num_procs-1];
	MPI_Status tempStats[num_procs-1];

	i = 0;
	while(i<num_procs-1){
		MPI_Isend(&ds.procs_msg_count[i+1],1,MPI_INT,i+1,0,MPI_COMM_WORLD,&tempReqs[i]);
		i++;
	}

	MPI_Waitall(num_procs-1,tempReqs,tempStats);
	datArray2 = (data_packet*)malloc(ds.procs_msg_count[0]*sizeof(data_packet));
	num_messages = ds.procs_msg_count[0];
	
	i = 0;
	idx = 0;
	count = 0;
	dangling_contri = 0.0;
	while(i<ds.n_pages){
		if(ds.map.find(i)==ds.map.end()){
			dangling_contri+=(*(ds.p_ranks + i));
		}else{
			int j = 0;
			while(j<ds.map[i].size()){
				struct data_packet send_packet;
				proc_id = ds.map[i][j];
				proc_id/=num_nodes;

				send_packet.value = (*(ds.p_ranks + i))/ds.map[i].size();
				send_packet.id = ds.map[i][j];
				if(proc_id!=0){
					MPI_Isend(&send_packet,1,data_packet_mpi,proc_id,0,MPI_COMM_WORLD,&reqArray[count]);
					count++;
				}else{
					*(datArray2 + idx) = send_packet;
					idx++;
				}
				j++;
			}
		}
		i++;
	}

	i = 1;
	dangling_contri/=ds.n_pages;
	while(i<num_procs){
		MPI_Isend(&dangling_contri,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqArray[count]);
		i++;
		count++;
	}
	MPI_Waitall(ds.n_msgs,reqArray,statusArray);
}

void reducer_function_(){
	MPI_Request tempReq1;
	MPI_Status tempStat1;
	MPI_Irecv(&num_messages,1,MPI_INT,0,0,MPI_COMM_WORLD,&tempReq1);
	
	int idx,rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Request tempReq2;
	MPI_Status tempStat2;

	MPI_Wait(&tempReq1, &tempStat1);
	
	struct data_packet data;
	datArray2 = (data_packet*)malloc(num_messages*sizeof(data_packet));
	MPI_Request recvArray[num_messages+1];
	MPI_Status statusArray[num_messages+1];

	idx = 0;
	for(idx = 0;idx<num_messages;idx++){
		MPI_Irecv(datArray2 + idx,1,data_packet_mpi,0,0,MPI_COMM_WORLD,&recvArray[idx]);
	}

	MPI_Irecv(&dangling_contri,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&recvArray[num_messages]);
	MPI_Waitall(num_messages+1,recvArray,statusArray);

	int i = 0;
	while(i<num_nodes){
		*(rank_array + i) = (damping_factor*dangling_contri);
		*(rank_array + i)+=(1.0 - damping_factor)/(num_procs*num_nodes);
		i++;
	}

	i = 0;
	while(i<num_messages){
		data = *(datArray2 + i);
		idx = data.id;
		idx-=rank*num_nodes;
		*(rank_array + idx)+=(damping_factor * data.value);
		i++;
	}

	MPI_Isend(rank_array,num_nodes,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&tempReq2);
	MPI_Wait(&tempReq2,&tempStat2);
}

int main(int argc, char *argv[]){
	string input_file_path = argv[1];
	string output_file_path = argv[3];
	
	graph_container_ds ds;
	initialize(ds,input_file_path);
	
	MPI_Init(&argc,&argv);

   	MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
   	const MPI_Aint displacement[2] = {0, sizeof(double)};
   	int lengths[2] = {1,1};
   	MPI_Type_create_struct(2,lengths,displacement,types,&data_packet_mpi);
   	MPI_Type_commit(&data_packet_mpi);

   	int i = 0;
   	int rank;
   	double curr_diff = 1.0;
   	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   	// Timer functionalities
   	double time_taken,time_taken1,time_taken2;
   	struct timeval start, end;
   	rank_array = (double*)malloc(ds.n_pages*sizeof(double));
   	gettimeofday(&start, NULL);

   	do{
   		if(rank!=0){
   			MPI_Request tempReq3;
   			MPI_Status tempStat3;
   			reducer_function_();
   			MPI_Irecv(&curr_diff,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&tempReq3);
   			MPI_Wait(&tempReq3, &tempStat3);
   		}else{
   			MPI_Request reqArray[ds.n_msgs];
   			MPI_Request reqArray2[num_procs-1];
   			MPI_Request reqArray3[num_procs-1];
   			
   			MPI_Status statusArray[ds.n_msgs];
   			MPI_Status statArray2[num_procs-1];
   			MPI_Status statArray3[num_procs-1];

   			map_function(ds, reqArray, statusArray);
   			reducer_function_2();
   			
   			int proc_id = 1;
   			while(proc_id<num_procs){
   				MPI_Irecv(rank_array + num_nodes*proc_id,num_nodes,MPI_DOUBLE,proc_id,0,MPI_COMM_WORLD,&reqArray2[proc_id-1]);
   				proc_id++;
   			}
   			MPI_Waitall(num_procs-1,reqArray2,statArray2);

   			curr_diff = 0.0;
   			double sum = 0.0;
   			int index = 0;
   			while(index<ds.n_pages){
   				sum+=(*(rank_array + index));
   				curr_diff+=fabs(*(ds.p_ranks+index) - *(rank_array+index));
   				*(ds.p_ranks + index) = (*(rank_array + index));
   				index++;
   			}

   			cout << i << "->" << curr_diff << ", " << sum << endl;
   			proc_id = 1;
   			while(proc_id<num_procs){
   				MPI_Isend(&curr_diff,1,MPI_DOUBLE,proc_id,0,MPI_COMM_WORLD,&reqArray3[proc_id-1]);
   				proc_id++;
   			}
   			
   			if(curr_diff<limit || i==max_iterations-1){
   				gettimeofday(&end, NULL);
   				write_data(ds.p_ranks,ds.n_pages,output_file_path);
   				time_taken1 = (end.tv_sec - start.tv_sec);
   				time_taken2 = (end.tv_usec - start.tv_usec) * 1e-6;
   				time_taken = time_taken1 + time_taken2;
   				cout << "Time taken " << time_taken << " sec" << endl; 
   			}else{
   				// do nothing
   			}

   			MPI_Waitall(num_procs-1,reqArray3,statArray3);
   		}
   		i++;
   	}while(i<max_iterations && curr_diff>limit);

	MPI_Finalize();
	return 0;
}