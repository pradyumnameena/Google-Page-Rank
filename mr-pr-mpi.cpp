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

// global variables
int num_procs = 2;
int num_nodes = 0;
double limit = 0.00001;
int max_iterations = 1000;
double damping_factor = 0.85;

int num_messages;
double* rank_array;
double dangling_contri;
hybrid_data* datArray2;

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
			sum+=(*(matrix + i));
		}
		file << "sum " << sum;
	}

	file.close();
}

void initialize(ultron &ds, string file_path){
	int n = -1;
	int count = 0;
	int a,b,proc_id;

	for(int i = 0;i<num_procs;i++){
		ds.procs_msg_count.push_back(0);
	}

	ifstream data_file(file_path);
	while(data_file >> a >> b){
		ds.map[a].push_back(b);
		n = max(n,max(a,b));
	}
	n+=1;
	num_nodes = n/num_procs;

	ds.p_ranks = (double*)malloc(n*sizeof(double));
	for(int i = 0;i<n;i++){
		*(ds.p_ranks + i) = 1.0/n;
		if(ds.map.find(i)!=ds.map.end()){
			for(int j = 0;j<ds.map[i].size();j++){
				proc_id = (ds.map[i][j])/num_nodes;
				if(ds.map[i][j]%num_nodes==0 && proc_id>1){
					proc_id--;
				}
				ds.procs_msg_count[proc_id]++;
				if(proc_id!=0){
					count++;
				}
			}
		}
	}

	ds.n_pages = n;
	ds.n_msgs = count + (num_procs-1);
	data_file.close();
}

void mapper(ultron ds, MPI_Request reqArray[], MPI_Status statusArray[]){
	int proc_id = 0;
	int count = 0;
	int idx = 0;
	dangling_contri = 0.0;

	MPI_Request tempReqs[num_procs-1];
	MPI_Status tempStats[num_procs-1];
	for(int i = 0;i<num_procs-1;i++){
		MPI_Isend(&ds.procs_msg_count[i+1],1,MPI_INT,i+1,0,MPI_COMM_WORLD,&tempReqs[i]);
	}

	num_messages = ds.procs_msg_count[0];
	datArray2 = (hybrid_data*)malloc(ds.procs_msg_count[0]*sizeof(hybrid_data));
	MPI_Waitall(num_procs-1,tempReqs,tempStats);
	
	for(int i = 0;i<ds.n_pages;i++){
		if(ds.map.find(i)!=ds.map.end()){
			for(int j = 0;j<ds.map[i].size();j++){
				proc_id = (ds.map[i][j]/num_nodes);

				struct hybrid_data send_packet;
				send_packet.id = ds.map[i][j];
				send_packet.value = (*(ds.p_ranks + i))/ds.map[i].size();
				if(proc_id==0){
					*(datArray2 + idx) = send_packet;
					idx++;
				}else{
					MPI_Isend(&send_packet,1,hybrid_data_type,proc_id,0,MPI_COMM_WORLD,&reqArray[count]);
					count++;
				}
			}
		}else{
			dangling_contri+=(*(ds.p_ranks + i));
		}
	}
	dangling_contri/=ds.n_pages;

	for(int i = 1;i<num_procs;i++){
		MPI_Isend(&dangling_contri,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqArray[count]);
		count++;
	}
	MPI_Waitall(ds.n_msgs,reqArray,statusArray);
}

void reduce(){
	int idx,rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Request tempReq1,tempReq2;
	MPI_Status tempStat1, tempStat2;

	MPI_Irecv(&num_messages,1,MPI_INT,0,0,MPI_COMM_WORLD,&tempReq1);
	MPI_Wait(&tempReq1, &tempStat1);
	
	MPI_Request recvArray[num_messages+1];
	MPI_Status statusArray[num_messages+1];
	datArray2 = (hybrid_data*)malloc(num_messages*sizeof(hybrid_data));
	struct hybrid_data data;

	idx = 0;
	while(idx<num_messages){
		MPI_Irecv(datArray2 + idx,1,hybrid_data_type,0,0,MPI_COMM_WORLD,&recvArray[idx]);
		idx++;
	}

	MPI_Irecv(&dangling_contri,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&recvArray[idx]);
	MPI_Waitall(num_messages+1,recvArray,statusArray);

	for(int i = 0;i<num_nodes;i++){
		*(rank_array + i) = (damping_factor*dangling_contri) + (1.0 - damping_factor)/(num_procs*num_nodes);
	}

	for(int i = 0;i<num_messages;i++){
		data = *(datArray2 + i);
		idx = data.id - rank*num_nodes;
		*(rank_array + idx)+=(damping_factor * data.value);
	}

	MPI_Isend(rank_array,num_nodes,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&tempReq2);
	MPI_Wait(&tempReq2,&tempStat2);
}

void reduce2(){
	struct hybrid_data recv_packet;
	int idx;

	for(int i = 0;i<num_nodes;i++){
		*(rank_array + i) = (dangling_contri*damping_factor) + (1.0 - damping_factor)/(num_procs*num_nodes);
	}

	for(int i = 0;i<num_messages;i++){
		recv_packet = *(datArray2 + i);
		idx = recv_packet.id;
		*(rank_array + idx)+=(damping_factor * recv_packet.value);
	}
}

int main(int argc, char *argv[]){
	// command is ./mr-pr-mpi.o ./input-file.txt -o ./output-file.txt -np number_of_processors
	string input_file_path = argv[1];
	string output_file_path = argv[3];
	if(argc>4){
		// num_procs = atoi(argv[5]);
		max_iterations = atoi(argv[4]);
	}
	
	ultron ds;
	initialize(ds,input_file_path);
	
	MPI_Init(&argc,&argv);

   	int lengths[2] = {1,1};
   	const MPI_Aint displacement[2] = {0, sizeof(double)};
   	MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
   	MPI_Type_create_struct(2,lengths,displacement,types,&hybrid_data_type);
   	MPI_Type_commit(&hybrid_data_type);

   	int i = 0;
   	int rank,size;
   	double curr_diff = 1.0;
   	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   	MPI_Comm_size(MPI_COMM_WORLD,&size);
   	rank_array = (double*)malloc(ds.n_pages*sizeof(double));

   	// Timer functionalities
   	struct timeval start, end;
   	double time_taken;
   	gettimeofday(&start, NULL);

   	while(i<max_iterations && curr_diff>limit){
   		if(rank==0){
   			MPI_Request reqArray[ds.n_msgs];
   			MPI_Status statusArray[ds.n_msgs];

   			mapper(ds, reqArray, statusArray);
   			reduce2();

   			MPI_Request reqArray2[num_procs-1];
   			MPI_Status statArray2[num_procs-1];
   			
   			for(int proc_id = 1;proc_id<num_procs;proc_id++){
   				MPI_Irecv(rank_array + num_nodes*proc_id,num_nodes,MPI_DOUBLE,proc_id,0,MPI_COMM_WORLD,&reqArray2[proc_id-1]);
   			}
   			MPI_Waitall(num_procs-1,reqArray2,statArray2);

   			curr_diff = 0.0;
   			double sum = 0.0;
   			for(int index = 0;index<ds.n_pages;index++){
   				curr_diff+=fabs(*(ds.p_ranks+index) - *(rank_array+index));
   				*(ds.p_ranks + index) = (*(rank_array + index));
   				sum+=(*(rank_array + index));
   			}

   			// for(int index = 0;index<ds.n_pages;index++){
   			// 	*(ds.p_ranks + index) = (*(ds.p_ranks + index))/sum;
   			// }
   			
   			cout << i << "->" << curr_diff << ", " << sum << endl;
   			MPI_Request reqArray3[num_procs-1];
   			MPI_Status statArray3[num_procs-1];
   			for(int proc_id = 1;proc_id<num_procs;proc_id++){
   				// cout << proc_id << ",";
   				MPI_Isend(&curr_diff,1,MPI_DOUBLE,proc_id,0,MPI_COMM_WORLD,&reqArray3[proc_id-1]);
   			}

   			MPI_Waitall(num_procs-1,reqArray3,statArray3);
   			
   			if(curr_diff<limit || i==max_iterations-1){
   				gettimeofday(&end, NULL);
   				time_taken = (end.tv_sec - start.tv_sec) * 1e6;
   				time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
   				cout << "Time taken " << time_taken << " sec" << endl; 
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