// #include "string.h"
// #include "sys/stat.h"
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


int global_n = -1;
#define LAMBDA 0.85

std::fstream file;
std::string filename;

//node struct
struct node
{
  std::vector<int> neighbours;
  float rank;
  int id;
};

std::map<int, node> graph;
std::vector<int> sink_nodes;
std::vector<int> source_nodes;
std::vector<float> ranks;

// functions
void mymap(int, KeyValue *, void *);
void myreduce(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

int main(int argc,char** argv)
{
    MPI_Init(&argc, &argv);
    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;
    if (argc <= 1) {
      if (rank == 0) printf("Syntax: pagerank filename ...\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    if(rank == 0){
      // calculating n
      int a, b;
      filename = argv[1];
      file.open(filename.c_str());
      while(file >> a >> b){
        global_n = std::max(a,global_n);
        global_n = std::max(b,global_n);
      }
      file.close();

      // assigning default pagerank values to each node
      for(int i = 0; i<global_n+1; i++){
        node temp;
        temp.id = i;
        if (i==0){
          temp.rank = 1.0;
        }
        else{
          temp.rank = 0.0;
        }
        graph[i] = temp;
      }

      // creating adjacency list
      file.open(filename.c_str());
      while(file >> a >> b){
        graph[a].neighbours.push_back(b);
      }
      file.close();

      // calculating sink and source nodes
      int* temp_count;
      temp_count = (int*)malloc((global_n+1)*sizeof(int));
      for(int i = 0;i<=global_n;i++){
        *(temp_count + i) = 0;
      }
      // int temp_count[global_n + 1] = {0};

      // sink_nodes
      std::map<int, node>::iterator it1 = graph.begin();
      while(it1 != graph.end()){
        node temp = it1->second;
        if(temp.neighbours.size() == 0){
            sink_nodes.push_back(it1->first);
        }
        for(int i =0; i< temp.neighbours.size(); i++){
          *(temp_count + temp.neighbours.at(i))+=1;
          // temp_count[temp.neighbours.at(i)] += 1;
        }
        it1++;
      }

      // source_nodes
      for(int i=0; i<global_n+1; i++){
        if(*(temp_count + i) == 0){
          source_nodes.push_back(i);
        }
      }
    }



    MapReduce *mr = new MapReduce(MPI_COMM_WORLD);
    mr->verbosity = 2;
    mr->timer = 1;

    // convergence params
    bool flag = true;
    float error_rank = 0;
    float temp_error = 0;
    float limit = 0.00001;

    for(int i=0; i<global_n+1; i++){
      if (i==0){
        ranks.push_back(1);
      }
      else{
        ranks.push_back(0);
      }
    }

    // handling dangling nodes
    float sum_ranks, sink_ranks;

    // loop counter
    int count = 0;


    MPI_Barrier(MPI_COMM_WORLD);
    float time = 0;


    // LOOP
    while(flag && count < 200){
      sum_ranks = 0;
      sink_ranks = 0;

      // calculating the mass lost due to sink nodes
      for(int k = 0; k < global_n+1; k++){
        float cpr = ranks[k];
        sum_ranks += cpr;
        if (graph[k].neighbours.size() == 0) {
            sink_ranks += cpr;
        }
      }

      // handling dangling nodes
      float term1 = LAMBDA*sink_ranks/(float)(global_n+1);
      // random walker
      float term2 = (1-LAMBDA)/(float)(global_n+1);

      // differnce for convergence
      error_rank = 0;



      double tstart = MPI_Wtime();
      int nnodes = mr->map(argc-1, mymap, NULL);
      mr->collate(NULL);
      int n_unique = mr->reduce(myreduce, NULL);
      MPI_Barrier(MPI_COMM_WORLD);
      double tstop = MPI_Wtime();
      delete mr;

      if (rank == 0){
        time+= (tstop - tstart);
        for(int i=0; i<global_n+1; i++){
          ranks[i] = ranks[i] + term1 + term2;
        }
        for(int i=0; i<global_n+1; i++){
          error_rank += std::fabs(ranks[i] - graph[i].rank);
          graph[i].rank = ranks[i];
        }
        // convergence
        if(fabs(error_rank - temp_error) < limit){
          flag = false;
        }
        else{
          temp_error = error_rank;
        }
        count++;
      }
    }
    if(rank==0){
      std::cout<< "\n Time ->>"<<time<< " Iterations "<<count <<" Number of nodes-->"<<global_n+1<<'\n';


      std::ofstream myfile;
      std::string answer_temp = "-pr-cpp.txt";
      std::string answer;
      answer = filename.substr(0,filename.length()-4) + answer_temp;
      myfile.open(answer);

      float sum = 0;
      // print final ranks
      for(int i = 0; i<global_n+1; i++){
        myfile << i << "-->" << graph[i].rank << std::endl;
        sum+=ranks[i];
      }
      myfile<<"sum-->"<<sum<<std::endl;
      myfile.close();
    }

    MPI_Finalize();
}

void mymap(int itask, KeyValue *kv, void *ptr){
  std::map<int, node>::iterator it2 = graph.begin();
  while(it2 != graph.end()){
    std::vector<int> temp_neighbours = (it2->second).neighbours;
    float temp_rank = (it2->second).rank;
    for(int i =0; i< temp_neighbours.size(); i++){
      std::string key = std::to_string(temp_neighbours.at(i));
      float temp_val = temp_rank/(float)temp_neighbours.size();
      std::string val = std::to_string(temp_val);
      kv->add((char *)key.c_str(), (int)sizeof(key), (char *)val.c_str(), (int)sizeof(val));
    }
    if (it2->first ==0){
      for(int i = 0; i<source_nodes.size(); i++){
        std::string key = std::to_string(source_nodes.at(i));
        float temp_val = 0.0;
        std::string val = std::to_string(temp_val);
        kv->add((char *)key.c_str(), (int)sizeof(key), (char *)val.c_str(), (int)sizeof(val));
      }
    }
    it2++;
  }
}


void myreduce(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr){
    float final = 0;
    for(int i=0; i<nvalues;i++){
      final+= std::atof((const char *)(multivalue + i));
    }
    ranks[std::atoi(key)] = final;
}
