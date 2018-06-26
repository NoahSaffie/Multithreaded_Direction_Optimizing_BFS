#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include "wtime.h"

//Timing
/*
  Failed Ideas/Methods:
  std::bitset can not be initialized with a variable (only constant)
*/
/*
  Current Issues:
  1 Still Reachable memory leak
  6 "Memory leaks" from OpenMP
  Future Issues/Questions:
*/
enum NextToRun { Top_Down, Bottom_Up };
typedef graph<long, long, int, long, long, char> Graph;
inline NextToRun bfs_top_down(Graph *ginst, int &edgesFrontier, int &edgesExplored, int* level, std::vector<int> &edges, std::vector<int> &nodes, int &current_level, const int number_threads, double alpha);
inline NextToRun bfs_bottom_up(Graph *ginst, int &edgesFrontier, int &edgesExplored, unsigned num_vertices, int* level, std::vector<int> &edges, std::vector<int> &nodes, int &current_level, const int number_threads, double beta);
int main(int argc, char* argv[])
{
  if(!omp_get_cancellation())
    {
      char *env = (char*)"OMP_CANCELLATION=true";
      printf("Cancellations were not enabled, enabling cancellation and rerunning program\n");
      putenv(env);
      execv(argv[0], argv);
    }
  //Graph Init
  std::cout << "Input: ./exe beg csr weight number_threads root(-1 if not sure) (optional)alpha (optional)beta\n";
  if (argc != 6 && argc != 8) { std::cout << "Wrong input\n"; return -1; }
  const char *beg_file = argv[1]; //Vertices -- Like RowPtr
  const char *csr_file = argv[2]; //Edges -- That lead to the Node described from the Node (index+1) in Offsets array (Beg_file)
  const char *weight_file = argv[3]; //Weight -- Unused
  const int number_threads = atoi(argv[4]);
  Graph *ginst = new Graph(beg_file, csr_file, weight_file);
  int root = atoi(argv[5]);
  double alpha = .15;
  double beta = 20;
  if(argc == 8)
  {
    alpha = atof(argv[6]);
    beta = atof(argv[7]);
  }

  //General Init
  std::vector<int> nodes;
  std::vector<int> edges;
  nodes.push_back(0); edges.push_back(0); //First Level with only root
  unsigned num_vertices = 1;
  int edgesExplored = 0, edgesFrontier = 0;
  int* level = new int[ginst->vert_count];
  memset(level, -1, sizeof(int)*ginst->vert_count);
  if(root == -1)
  {
    //Finds a root
    for(int i = 0; i < ginst->vert_count; i++)
        {
        //Once we have found the index in beg_pos that is no longer looking at the first index of edge list (->csr), then we are currently sitting on the true index that the adjency list
        //begining at ->csr[0] belongs to
            if(ginst->beg_pos[i+1] != 0) //We skip any "missing" vertexs (no edges) by skipping all but the last vertex that has beg_pos[i] of 0
            {
                level[i] = 0; root =i;
                break;
            }
        }
    }
  printf("Root: %d\n", root);
  nodes.at(0) = 1;
  edges.at(0) = (ginst->beg_pos[root+1] - ginst->beg_pos[root]);
  edgesFrontier += edges.at(0);

  //Find true number of verticies
  int last_num = 0;
  for(int j = 1; j < ginst->vert_count; ++j)
    {
      if(ginst->beg_pos[j] > last_num)
        {
          last_num = ginst->beg_pos[j];
          ++num_vertices;
        }
    }
  //Stats - Init
  int edges_TD = 0, edges_BU = 0, edges_base = 0;
  double time_TD = 0.0, time_BU = 0.0;
  double start = wtime();
  NextToRun next = Top_Down;
  int current_level = 0;
  unsigned total_nodes = 0;
  //Main Loop
  while (num_vertices > total_nodes)
    {
      if(next == Top_Down)
        {
          double start_TD = wtime();
          next = bfs_top_down(ginst, edgesFrontier, edgesExplored, level, edges, nodes, current_level, number_threads, alpha);
          double end_TD = wtime();
          edges_TD += (edgesExplored-edges_base); edges_base = edgesExplored;
          time_TD += end_TD-start_TD;
        }
      if(next == Bottom_Up)
        {
          double start_BU = wtime();
          next = bfs_bottom_up(ginst, edgesFrontier, edgesExplored, num_vertices, level, edges, nodes, current_level, number_threads, beta);
          double end_BU = wtime();
          edges_BU += (edgesExplored-edges_base); edges_base = edgesExplored;
          time_BU += end_BU-start_BU;
        }
      total_nodes = 0;
      for(unsigned i = 0; i < nodes.size(); ++i) { total_nodes+= nodes.at(i); } //printf("Total nodes: %d\n", total_nodes);
    }
  double end = wtime();
  double time_total = end-start;
  for(unsigned i = 0; i < nodes.size(); ++i)
    {
      if(nodes.at(i) != 0)
	    {
	        printf("Level: %d\tFrontier Size: %d\tEdges: %d\n", i, nodes.at(i), edges.at(i));
	    }
    }
  std::cout << "Edges: " << "\t\t%Top_Down: " << 100*(edges_TD*1.0/edgesExplored) << "\t%Bottom_Up: " << 100*(edges_BU*1.0/edgesExplored) << "\tTotal: " << edgesExplored << std::endl;
  std::cout << "Timing(seconds): " << "\t%Top_Down: " << 100*(time_TD/time_total) << "\t%Bottom_Up: " << 100*(time_BU/time_total) << "\tTotal: " << time_total << std::endl;

  // Cleanup
  free(ginst->beg_pos);
  free(ginst->csr);
  free(ginst->weight);
  delete ginst;
  delete[] level;
  return 0;
}






inline NextToRun bfs_top_down(Graph *ginst, int &edgesFrontier, int &edgesExplored, int* level, std::vector<int> &edges, std::vector<int> &nodes, int &current_level, const int number_threads, double alpha)
{
  int* partial_edgesExplored = new int[number_threads];  assert(partial_edgesExplored != NULL);
  int* partial_edgesFrontier = new int[number_threads]; assert(partial_edgesFrontier != NULL);
  int* partial_nodes = new int[number_threads]; assert(partial_nodes != NULL);
  memset(partial_edgesExplored, 0, sizeof(int)*number_threads);
  memset(partial_edgesFrontier, 0, sizeof(int)*number_threads);
  memset(partial_nodes, 0, sizeof(int)*number_threads);
  #pragma omp parallel num_threads(number_threads) default(none) firstprivate(ginst, alpha, number_threads, partial_edgesExplored, partial_nodes, partial_edgesFrontier) shared(level, current_level, edges, nodes, edgesFrontier, edgesExplored)
  {
    const int tid = omp_get_thread_num();
    const int section_size = ginst->vert_count/number_threads;
    const int section_beg = tid*section_size;
    const int section_end = (tid==number_threads-1 ? ginst->vert_count:section_beg+section_size);
    while(true)
      {
        for (int vert_index = section_beg; vert_index < section_end; vert_index++) //Go through your section
          {
            if (level[vert_index]==current_level) // If unexplored
              {
                int neighbor_beg = ginst->beg_pos[vert_index];
                int neighbor_end = ginst->beg_pos[vert_index+1];
                for(int neighbor_index = neighbor_beg; neighbor_index < neighbor_end; neighbor_index++) //Go through neighbors
                  {
                    int node = ginst->csr[neighbor_index];
                    partial_edgesFrontier[tid] += (ginst->beg_pos[node+1] - ginst->beg_pos[node]);
                    ++(partial_edgesExplored[tid]);
                    #pragma omp critical
                    if(level[node] == -1) //If neighbor is unexplored
                      {
                        level[node] = current_level+1; //Add to next frontier
                        ++(partial_nodes[tid]);
                      }
                  }
              }
          }
        #pragma omp barrier //Wait for other threads to finish their section
        #pragma omp single //Stats updating and level update
        {
          for(int j = 0; j < number_threads; j++)
            {
              edges.at(current_level) += partial_edgesExplored[j];
              edgesFrontier += partial_edgesFrontier[j];
              nodes.at(current_level) += partial_nodes[j];
            }
          edgesExplored += edges.at(current_level);
          memset(partial_edgesExplored, 0, sizeof(int)*number_threads);
          memset(partial_edgesFrontier, 0, sizeof(int)*number_threads);
          memset(partial_nodes, 0, sizeof(int)*number_threads);
          ++current_level; edges.push_back(0); nodes.push_back(0);
        }
        if(nodes.at(current_level-1) == 0)
          {
            #pragma omp cancel parallel
          }
        if(edgesFrontier/ginst->edge_count > alpha)
          {
            #pragma omp cancel parallel
          }
      }
  }
  //No longer in parallel region so just be one thread
  delete[] partial_edgesExplored;
  delete[] partial_edgesFrontier;
  delete[] partial_nodes;
  return (edgesFrontier/ginst->edge_count > alpha) ? Bottom_Up:Top_Down;
}
inline NextToRun bfs_bottom_up(Graph *ginst, int &edgesFrontier, int& edgesExplored, unsigned num_vertices, int* level, std::vector<int> &edges, std::vector<int> &nodes, int &current_level, const int number_threads, double beta)
{
  edgesFrontier = 0;
  int* partial_edgesExplored = new int[number_threads]; assert(partial_edgesExplored != NULL);
  int* partial_edgesFrontier = new int[number_threads]; assert(partial_edgesFrontier != NULL);
  int* partial_nodes = new int[number_threads]; assert(partial_nodes != NULL);
  //Init
  memset(partial_nodes, 0, sizeof(int)*number_threads);
  memset(partial_edgesExplored, 0, sizeof(int)*number_threads);
  memset(partial_edgesFrontier, 0, sizeof(int)*number_threads);
  #pragma omp parallel num_threads(number_threads) default(none) firstprivate(ginst, number_threads, partial_edgesExplored, partial_edgesFrontier, partial_nodes) shared(level, current_level)
  {
    const int tid = omp_get_thread_num();
    const int section_size = ginst->vert_count/number_threads;
    const int section_beg = tid*section_size;
    const int section_end = (tid==number_threads-1 ? ginst->vert_count:section_beg+section_size);
    for (int vert_index = section_beg; vert_index < section_end; ++vert_index)
      {
        if (level[vert_index]==-1) //Unexplored (Only want to look through potential children of frontier)
          {
            long int neighbor_beg = ginst->beg_pos[vert_index];
            long int neighbor_end = ginst->beg_pos[vert_index+1];
            for (long int neighbor_index = neighbor_beg; neighbor_index < neighbor_end; ++neighbor_index) //Indexes in CSR
              {
                long int node = ginst->csr[neighbor_index]; //Neighbor of vert_index
                ++partial_edgesExplored[tid];
                if (level[node] == current_level) //Visited neighbor limited to current_level
                  {
                    level[vert_index] = current_level+1;
                    ++(partial_nodes[tid]); //Should be 1 for next_level and 0 for previous level (==current_level)
                    partial_edgesFrontier[tid] += (ginst->beg_pos[vert_index+1] - ginst->beg_pos[vert_index]);
                    break;
                  }
              }
          }
      }
  }
  for(int j = 0; j < number_threads; ++j)
    {
      nodes.at(current_level) += partial_nodes[j];
      edges.at(current_level) += partial_edgesExplored[j];
      edgesFrontier += partial_edgesFrontier[j];
    }
  edgesExplored += edges.at(current_level);
  current_level++; nodes.push_back(0); edges.push_back(0);
  delete[] partial_nodes;
  delete[] partial_edgesExplored;
  delete[] partial_edgesFrontier;
  return nodes.at(current_level)/num_vertices <  beta ? Top_Down:Bottom_Up;
}
