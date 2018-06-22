#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include <chrono>
#define alpha 2.77  // Was 12 but found it favored Bottom-Up way too much - 1 keeps it from switching but is always close to switching
// -- Change -- After neeeding to use Reversing of Edges, bumping this up became beneficial
#define beta 24

//Timing
typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()
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
inline NextToRun bfs_top_down(Graph *ginst, int &edgesFrontier, int &edgesExplored, int* level, std::vector<int> &edges, std::vector<int> &nodes, int &current_level);
inline NextToRun bfs_bottom_up(Graph *ginst, int &edgesFrontier, int &edgesExplored, unsigned num_vertices, int* level, std::vector<int> &edges, std::vector<int> &nodes, int &current_level);
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
  std::cout << "Input: ./exe beg csr weight\n";
  if (argc != 4) { std::cout << "Wrong input\n"; return -1; }
  const char *beg_file = argv[1]; //Vertices -- Like RowPtr
  const char *csr_file = argv[2]; //Edges -- That lead to the Node described from the Node (index+1) in Offsets array (Beg_file)
  const char *weight_file = argv[3]; //Weight -- Unused
  Graph *ginst = new Graph(beg_file, csr_file, weight_file);
  //General Init

  std::vector<int> nodes;
  std::vector<int> edges; //Do we want this to entirely replace edgesExplored? If we do that then we would need to sum through the vector every time we wanted edgesExplored (But this isn't very costly)
  nodes.push_back(0); edges.push_back(0); //First Level with only root
  nodes.push_back(0); edges.push_back(0); //Where first frontier begins in Top_Down
  unsigned num_vertices = 1;
  int edgesExplored = 0, edgesFrontier = 0;
  int* level = new int[ginst->vert_count];
  memset(level, -1, sizeof(int)*ginst->vert_count);
  //Frontier - Holds Node/Vertex # which also act as an Index in beg_pos
  //Init the frontier - Find the lowest Node/Vertex that has an adjacency
  int root = -1;
  for(int i = 0; i < ginst->vert_count; i++)
    {
      //Once we have found the index in beg_pos that is no longer looking at the first index of edge list (->csr), then we are currently sitting on the true index that the adjency list
      //begining at ->csr[0] belongs to
      if(ginst->beg_pos[i+1] != 0) //We skip any "missing" vertexs (no edges) by skipping all but the last vertex that has beg_pos[i] of 0
        {
          printf("Root: %d\n", i);
          level[i] = 0; root =i;
          break;
        }
    } assert(root >= 0);
  nodes.at(0) = 1;
  edges.at(0) = edgesFrontier += (ginst->beg_pos[root+1] - ginst->beg_pos[root]);
  //Find true # of vertices
  //Can loop through and count how many different numbers we see for a accurate value (Confirmed accurate)
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
  TimeVar start = timeNow();
  NextToRun next = Top_Down;
  int current_level = 0;
  unsigned total_nodes = 0;
  while (num_vertices > total_nodes)
    {
      if(next == Top_Down)
        {
          TimeVar start_TD = timeNow();
          next = bfs_top_down(ginst, edgesFrontier, edgesExplored, level, edges, nodes, current_level);
          TimeVar end_TD = timeNow();
          edges_TD += (edgesExplored-edges_base); edges_base = edgesExplored;
          time_TD += duration(end_TD-start_TD);
        }
      if(next == Bottom_Up)
        {
          TimeVar start_BU = timeNow();
          next = bfs_bottom_up(ginst, edgesFrontier, edgesExplored, num_vertices, level, edges, nodes, current_level);
          TimeVar end_BU = timeNow();
          edges_BU += (edgesExplored-edges_base); edges_base = edgesExplored;
          time_BU += duration(end_BU-start_BU);
        }
      total_nodes = 0;
      for(unsigned i = 0; i < nodes.size(); ++i){total_nodes+= nodes.at(i);}
      //printf("Nodes: %d\n", total_nodes);
    }
  TimeVar end = timeNow();
  double time_total = duration(end-start);
  for(unsigned i = 0; i < nodes.size(); ++i)
    {
      printf("Level: %d\tFrontier Size: %d\tEdges: %d\n", i, nodes.at(i), edges.at(i));
    }
  std::cout << "Edges: " << "\t\t%Top_Down: " << 100*(edges_TD*1.0/edgesExplored) << "\t%Bottom_Up: " << 100*(edges_BU*1.0/edgesExplored) << "\tTotal: " << edgesExplored << std::endl;
  std::cout << "Timing(ms): " << "\t%Top_Down: " << 100*(time_TD/time_total) << "\t%Bottom_Up: " << 100*(time_BU/time_total) << "\tTotal: " << (int)time_total << std::endl;

  // Cleanup
  free(ginst->beg_pos);
  free(ginst->csr);
  free(ginst->weight);
  delete ginst;
  delete[] level;
  return 0;
}
inline NextToRun bfs_top_down(Graph *ginst, int &edgesFrontier, int &edgesExplored, int* level, std::vector<int> &edges, std::vector<int> &nodes, int &current_level)
{
  //printf("Current_Level: %d\n", current_level);
  int number_threads = omp_get_num_threads(); //printf("Threads TD: %d\n", number_threads); //4 for this system {0,1,2,3}
  int* partial_edgesExplored = new int[number_threads];  assert(partial_edgesExplored != NULL);
  int* partial_edgesFrontier = new int[number_threads]; assert(partial_edgesFrontier != NULL);
  int* partial_nodes = new int[number_threads]; assert(partial_nodes != NULL);
#pragma omp parallel num_threads(number_threads) default(none) firstprivate(ginst, number_threads) shared(level, partial_edgesExplored, partial_nodes, current_level, partial_edgesFrontier, edges, nodes, edgesFrontier, edgesExplored)
  {
    const int tid = omp_get_thread_num();
    const int section_size = ginst->vert_count/number_threads;
    const int section_beg = tid*section_size;
    const int section_end = tid==number_threads-1 ? ginst->vert_count:section_beg+section_size;
    while(true)
      {
        memset(partial_edgesExplored, 0, sizeof(int)*number_threads);
        memset(partial_edgesFrontier, 0, sizeof(int)*number_threads);
        memset(partial_nodes, 0, sizeof(int)*number_threads);
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
                    ++partial_edgesExplored[tid];
                    if(level[node] == -1) //If neighbor is unexplored
                      {
                        level[node] = current_level+1; //Add to next frontier
                        ++partial_nodes[tid];
                      }
                  }
              }
          }
        #pragma omp barrier
        for(int j = 0; j < number_threads; j++)
          {
            edges.at(current_level) += partial_edgesExplored[j];
            edgesFrontier += partial_edgesFrontier[j];
            nodes.at(current_level) += partial_nodes[j];
          }
        edgesExplored += edges.at(current_level);
        if(nodes.at(current_level) == 0) break; //Absolute end condition (found all nodes)
        #pragma omp single
        {
          ++current_level; edges.push_back(0); nodes.push_back(0);
	}
        if(edgesFrontier > (ginst->edge_count-edgesExplored)/alpha)
          {
             #pragma omp cancel parallel
          }
        //Another barrier to ensure everyone is done here before moving back
	#pragma omp cancellation point parallel
        #pragma omp barrier
      }
  }
  delete[] partial_edgesExplored;
  delete[] partial_edgesFrontier;
  delete[] partial_nodes;
  return (edgesFrontier > (ginst->edge_count-edgesExplored)/alpha) ? Bottom_Up:Top_Down;
}
inline NextToRun bfs_bottom_up(Graph *ginst, int &edgesFrontier, int& edgesExplored, unsigned num_vertices, int* level, std::vector<int> &edges, std::vector<int> &nodes, int &current_level)
{
  //Potential inconsistent results at loop level, since if a loop finishes a incomplete level, and then comes across some of the neighbors to that incomplete part, then it will see them as parents
  //as well. BUt sometimes a loop may have already passed this point so it never the previously imcomplete part of previous level as parents (in the end answer should be consistent)
  //printf("Current_Level: %d\n", current_level);
  edgesFrontier = 0;
  int number_threads = omp_get_num_threads(); //printf("Threads BU: %d\n", omp_get_num_threads());
  int* partial_edgesExplored = new int[number_threads]; assert(partial_edgesExplored != NULL);
  int* partial_edgesFrontier = new int[number_threads]; assert(partial_edgesFrontier != NULL);
  //int* partial_nodes = new int[number_threads]; assert(partial_nodes != NULL); //1D
  int** partial_nodes = new int*[number_threads]; assert(partial_nodes != NULL);
  for(int i = 0; i < number_threads; ++i) //The second array is 2 wide because we only want to find nodes from the previous incomplete level (if that the case) or for the current level we are forming
    {
      partial_nodes[i] = new int[2]; assert(partial_nodes[i] != NULL);
      memset(partial_nodes[i], 0, sizeof(int)*2);
    }
  //Init
  memset(partial_edgesExplored, 0, sizeof(int)*number_threads);
  memset(partial_edgesFrontier, 0, sizeof(int)*number_threads);
#pragma omp parallel num_threads(number_threads) default(none) firstprivate(ginst, number_threads, current_level) shared(level, partial_edgesExplored, partial_edgesFrontier, partial_nodes)
  {
    const int tid = omp_get_thread_num();
    const int section_size = ginst->vert_count/number_threads;
    const int section_beg = tid*section_size;
    const int section_end = tid==number_threads-1 ? ginst->vert_count:section_beg+section_size;
    for (int vert_index = section_beg; vert_index < section_end; ++vert_index)
      {
        if (level[vert_index]==-1) //Unexplored (Only want to look through potential children of frontier)
          {
            int neighbor_beg = ginst->beg_pos[vert_index];
            int neighbor_end = ginst->beg_pos[vert_index+1];
            for (int neighbor_index = neighbor_beg; neighbor_beg < neighbor_end; ++neighbor_index) //Indexes in CSR
              {
		int node = ginst->csr[neighbor_index]; //Neighbor of vert_index
                ++partial_edgesExplored[tid];
                if (level[node] != -1 && (level[node] == current_level || level[node] == current_level-1)) //Visited neighbor limited to current_level or 1 before (if previous level is incomplete)
                  { 
                    level[vert_index] = level[node]+1;
                    ++(partial_nodes[tid][(current_level-level[vert_index]+1)]); //Should be 1 for next_level and 0 for previous level (==current_level)
                    //if(level[vert_index] == current_level || level[vert_index] == current_level-1)
		    partial_edgesFrontier[tid] += (ginst->beg_pos[vert_index+1] - ginst->beg_pos[vert_index]);
                    break;
                  }
              }
          }
      }
  }
  for(int j = 0; j < number_threads; ++j)
    {
      for(int i = 0; i < 2; ++i)
        {
          nodes.at(current_level+i) += partial_nodes[j][i];
        }
      edges.at(current_level) += partial_edgesExplored[j];
      edgesFrontier += partial_edgesFrontier[j];
    }
  edgesExplored += edges.at(current_level);
  for(int j = 0; j < number_threads; ++j)
    {
      delete [] partial_nodes[j];
    }
  current_level++; nodes.push_back(0); edges.push_back(0);
  delete[] partial_nodes;
  delete[] partial_edgesExplored;
  delete[] partial_edgesFrontier;
  return nodes.at(current_level) < ((int)num_vertices / beta) ? Top_Down:Bottom_Up;
}
