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
int num_procs = 2;
int num_nodes = 0;
double limit = 0.00001;
int max_iterations = 1000;
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
	set<int> secSet;

	for(int i = 0;i<=num_procs;i++){
		ds.procs_msg_count.push_back(0);
	}

	ifstream data_file(file_path);
	while(data_file >> a >> b){
		ds.map[a].push_back(b);
		secSet.insert(b);
		n = max(n,max(a,b));
		count+=1;

		proc_id = (b%num_procs)+1;
		ds.procs_msg_count[proc_id]++;
	}
	n+=1;

	ds.p_ranks = (double*)malloc(n*sizeof(double));
	for(int i = 0;i<n;i++){
		*(ds.p_ranks + i) = 1.0/n;
		if(secSet.find(i)==secSet.end()){
			proc_id = (i%num_procs) + 1;
			ds.noEntryNodes.push_back(i);
			ds.procs_msg_count[proc_id]++;
		}
	}

	ds.n_pages = n;
	ds.n_msgs = count;
	data_file.close();
}

void mapper(ultron ds, MPI_Request reqArray[]){
	int proc_id = 0;
	int count = 0;

	for(int i = 1;i<=num_procs;i++){
		MPI_Send(&ds.procs_msg_count[i],1,MPI_INT,i,0,MPI_COMM_WORLD);
	}
	
	for(int i = 0;i<ds.noEntryNodes.size();i++){
		proc_id = ds.noEntryNodes[i]%num_procs + 1;
		struct hybrid_data send_packet;
		send_packet.id = ds.noEntryNodes[i];
		send_packet.value = 0.0;
		MPI_Isend(&send_packet,1,hybrid_data_type,proc_id,0,MPI_COMM_WORLD,&reqArray[count]);
		count++;
	}

	for(int i = 0;i<ds.n_pages;i++){
		if(ds.map.find(i)!=ds.map.end()){
			for(int j = 0;j<ds.map[i].size();j++){
				proc_id = ds.map[i][j]%num_procs + 1;
				struct hybrid_data send_packet;
				send_packet.id = ds.map[i][j];
				send_packet.value = (*(ds.p_ranks + i))/ds.map[i].size();
				MPI_Isend(&send_packet,1,hybrid_data_type,proc_id,0,MPI_COMM_WORLD,&reqArray[count]);
				count++;
			}
		}
	}
}

void reduce(){
	int num_messages,num_nodes,idx,rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Recv(&num_messages,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	
	MPI_Request recvArray[num_messages];
	MPI_Status statusArray[num_messages];
	struct hybrid_data datArray[num_messages];

	while(idx<num_messages){
		MPI_Irecv(&datArray[idx],1,hybrid_data_type,0,0,MPI_COMM_WORLD,&recvArray[idx]);
		idx++;
	}

	MPI_Waitall(num_messages,recvArray,statusArray);

	for(int i = 0;i<num_nodes;i++){
		*(rank_array + i) = 0.0;
	}

	for(int i = 0;i<num_messages;i++){
		idx = datArray[i].id - rank*num_nodes;
		*(rank_array + idx)+=datArray[i].value;
	}

	MPI_Send(rank_array,num_nodes,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
}

int main(int argc, char *argv[]){
	string input_file_path = argv[1];
	string output_file_path = argv[2];
	num_procs = atoi(argv[4]) - 1;
	
	ultron ds;
	initialize(ds,input_file_path);
	num_nodes = ds.n_pages/num_procs;
	
	MPI_Init(&argc,&argv);

   	int lengths[2] = {1,1};
   	const MPI_Aint displacement[2] = {0, sizeof(double)};
   	MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
   	MPI_Type_create_struct(2,lengths,displacement,types,&hybrid_data_type);
   	MPI_Type_commit(&hybrid_data_type);

   	// ALGORITHM
   	int rank,i;
   	double curr_diff = 1.0;
   	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   	double* rank_array = (double*)malloc(ds.n_pages*sizeof(double));

   	while(i<max_iterations && curr_diff>limit){
   		if(rank==0){
   			MPI_Request reqArray[ds.n_msgs];
   			MPI_Status statusArray[ds.n_msgs];
   			mapper(ds, reqArray);
   			MPI_Waitall(ds.n_msgs,reqArray,statusArray);

   			for(int proc_id = 1;proc_id<=num_procs;proc_id++){
   				MPI_Recv(rank_array + (proc_id-1)*num_nodes,num_nodes,MPI_DOUBLE,proc_id,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   			}

   			curr_diff = 0.0;
   			for(int index = 0;index<ds.n_pages;index++){
   				curr_diff+=fabs(*(ds.p_ranks+index) - *(rank_array+index));
   				*(ds.p_ranks + index) = *(rank_array + index);
   			}
   			
   			for(int proc_id = 1;proc_id<=num_procs;proc_id++){
   				MPI_Send(&curr_diff,1,MPI_DOUBLE,proc_id,0,MPI_COMM_WORLD);
   			}
   			
   			if(curr_diff<limit || i==max_iterations-1){
   				write_data(rank_array,ds.n_pages,output_file_path);
   			}

   		}else{
   			reduce();
   			MPI_Recv(&curr_diff,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   		}
   		i++;
   	}

	MPI_Finalize();
	return 0;
}