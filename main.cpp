#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include <chrono>
#define alpha 2.77  // Was 12 but found it favored Bottom-Up way too much - 1 keeps it from switching but is always close to switching, so a little >1 is perfect
//@ 1.15:: Explored nodes 2789664/3072627, 12071273/117185083 edges, and Top-Down found 60 % of edges, taking up 52% of runtime (1838461 ms total)
//@ 1.17:: Explored nodes 2781249/3072627, 11776471/117185083 edges, and Top-Down found 56 % of edges, taking up 50% of runtime (1833991 ms total)
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
  Potential Room for Improvement:
  Finding # of vertexs,  can likely be improved
  Finding edgesFrontier
*/
/*
  Current Issues:
  1 Still Reachable memory leak
  6 "Memory leaks" from OpenMP
  Future Issues/Questions:
*/
enum NextToRun { Top_Down, Bottom_Up };
typedef graph<long, long, int, long, long, char> Graph;
inline NextToRun bfs_top_down(Graph *ginst, int &edgesFrontier,  int* level, std::vector<int> edges, std::vector<int> nodes, int &current_level);
inline NextToRun bfs_bottom_up(Graph *ginst, int &edgesFrontier,  unsigned num_vertices, int* level, std::vector<int> edges, std::vector<int> nodes, int &current_level);
int main(int argc, char* argv[])
{
  //Graph Init
  std::cout << "Input: ./exe beg csr weight\n";
  if (argc != 4) { std::cout << "Wrong input\n"; return -1; }
  const char *beg_file = argv[1]; //Vertices -- Like RowPtr
  const char *csr_file = argv[2]; //Edges -- That lead to the Node described from the Node (index+1) in Offsets array (Beg_file)
  const char *weight_file = argv[3]; //Weight -- Unused
  Graph *ginst = new Graph(beg_file, csr_file, weight_file);

  //General Init
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
  edgesFrontier += (ginst->beg_pos[root+1] - ginst->beg_pos[root]);

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
  //NextToRun next = Top_Down;
  std::vector<int> nodes;
  std::vector<int> edges; //Do we want this to entirely replace edgesExplored? If we do that then we would need to sum through the vector every time we wanted edgesExplored (But this isn't very costly)
  int current_level = 0;
  TimeVar start_TD = timeNow();
  bfs_top_down(ginst, edgesFrontier, level, edges, nodes, current_level);
  TimeVar end_TD = timeNow();
  edges_TD += (edgesExplored-edges_base); edges_base = edgesExplored;
  time_TD += duration(end_TD-start_TD);
  /*
    while (true)
    {
    // Will switch to the other when the previous one returns (which it will do depending on the value it recieves from switch function)
    // Requires one check extra check somewhere if we finish program on top_down, we need to know to skip bottom_up and break
    if(next == Top_Down)
    {
    TimeVar start_TD = timeNow();
    next = bfs_top_down(ginst, edgesExplored, edgesFrontier, level);
    TimeVar end_TD = timeNow();
    edges_TD += (edgesExplored-edges_base); edges_base = edgesExplored;
    time_TD += duration(end_TD-start_TD);
    }
    if(next == Bottom_Up)
    {
    TimeVar start_BU = timeNow();
    next = bfs_bottom_up(ginst, edgesExplored, edgesFrontier, num_vertices, level);
    TimeVar end_BU = timeNow();
    edges_BU += (edgesExplored-edges_base); edges_base = edgesExplored;

    time_BU += duration(end_BU-start_BU);
    }
    }*/
  TimeVar end = timeNow();
  double time_total = duration(end-start);
  //Nodes: 3072441
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
/*
  Use a vector for edges and nodes, since that information needs to be passed between functions and may be incomplete, during a pass, to be completed in other function (or maybe not even then)
 */
inline NextToRun bfs_top_down(Graph *ginst, int &edgesFrontier, int* level, std::vector<int> edges, std::vector<int> nodes, int &current_level)
{ //Consider seeing if we can move OpenMP up a level -- Would cause an issue with the switch check
  int num_threads = omp_get_num_threads(); //4 for this system {0,1,2,3}
  int* partial_edgesExplored = new int[num_threads];  assert(partial_edgesExplored != NULL);
  //int* partial_edgesFrontier = new int[num_threads]; assert(partial_edgesFrontier != NULL);
  int* partial_nodes = new int[num_threads]; assert(partial_nodes != NULL);
  //Init
  memset(partial_edgesExplored, 0, sizeof(int)*num_threads);
  //memset(partial_edgesFrontier, 0, sizeof(int)*num_threads);
  memset(partial_nodes, 0, sizeof(int)*num_threads);
#pragma omp parallel num_threads(num_threads) default(none) firstprivate(ginst, num_threads) shared(level, partial_edgesExplored, partial_nodes, current_level, edgesFrontier, edges, nodes)
  {
    const int tid = omp_get_thread_num();
    const int section_size = ginst->vert_count/num_threads;
    const int section_beg = tid*section_size;
    const int section_end = tid==num_threads-1 ? ginst->vert_count:section_beg+section_size;
    while(true)
      {
        memset(partial_edgesExplored, 0, sizeof(int)*num_threads);
        //memset(partial_edgesFrontier, 0, sizeof(int)*num_threads);
        memset(partial_nodes, 0, sizeof(int)*num_threads);
        for (int vert_index = section_beg; vert_index < section_end; vert_index++)
          {
            if (level[vert_index]==current_level) // If unexplored
              {
                int neighbor_beg = ginst->beg_pos[vert_index];
                int neighbor_end = ginst->beg_pos[vert_index+1];
                for(int neighbor_index = neighbor_beg; neighbor_index < neighbor_end; neighbor_index++)
                  {
                    int node = ginst->csr[neighbor_index];
                    //partial_edgesFrontier[tid] += (ginst->beg_pos[node+1] - ginst->beg_pos[node]);
                    ++partial_edgesExplored[tid];
                    if(level[node] == -1)
                      {
                        level[node] = current_level+1;
                        ++partial_nodes[tid];
                      }
                  }
              }
          }
        #pragma omp barrier
        //#pragma omp parallel for reduction(+:edges) reduction(+:edgesFrontier) reduction(+:nodes)
        for(int j = 0; j < num_threads; j++)
          {
            edges.at(current_level) += partial_edgesExplored[j];
            //edgesFrontier += partial_edgesFrontier[j];
            nodes.at(current_level) += partial_nodes[j];
          }
        if(nodes.at(current_level) == 0) break;
        if(tid == 0)
          {
            printf("Level: %d\tFrontier Size: %d\tEdges: %d\n", current_level, nodes.at(current_level), edges.at(current_level));
          }
        ++current_level;
      }
  }
  delete[] partial_edgesExplored;
  //delete[] partial_edgesFrontier;
  delete[] partial_nodes;
  return Top_Down;
}
inline NextToRun bfs_bottom_up(Graph *ginst, int &edgesFrontier, unsigned num_vertices, int* level, std::vector<int> edges, std::vector<int> nodes, int &current_level)
{
  int num_threads = omp_get_num_threads();
  int* partial_edgesExplored = new int[num_threads];
  int* partial_edgesFrontier = new int[num_threads];
  //Init
  memset(partial_edgesExplored, 0, sizeof(int)*num_threads);
  memset(partial_edgesFrontier, 0, sizeof(int)*num_threads);
#pragma omp parallel num_threads(num_threads) default(none) firstprivate(ginst) shared(current_level, partial_edgesFrontier, nodes, edges)
  {
    const int tid = omp_get_thread_num();
    const int section_size = ginst->vert_count/num_threads;
    const int section_beg = tid*section_size;
    const int section_end = tid==num_threads-1 ? ginst->vert_count:section_beg+section_size;
    for (int vert_index = section_beg; vert_index < section_end; ++vert_index)
      {
        if (level[vert_index]==-1) //Unexplored (Only want to look through potential children of frontier)
          {
	    int neighbor_beg = ginst->beg_pos[vert_index];
            int neighbor_end = ginst->beg_pos[vert_index+1];
            for (int neighbor_index = neighbor_beg; neighbor_beg < neighbor_end; ++neighbor_index) //Indexes in CSR
              {
                int node = ginst->csr[neighbor_index]; //Neighbor of i
                ++partial_edgesExplored[tid];
                if (level[node] == current_level-1) //Find a neighbor that has been visited, it will be a parent (on frontier too)
                  {
                    level[vert_index] = current_level; // Now that the unexplored vertice has found a parent amonst it's neighbors/connections we can add it to the frontier (and will later mark as visited)
                    partial_edgesFrontier[tid] += (ginst->beg_pos[vert_index+1] - ginst->beg_pos[vert_index]);
                    break;
                  }
              }
          }
      }
  }
  //#pragma omp parallel for reduction(+:edgesExplored) reduction(+:edgesFrontier)
  for(int j = 0; j < num_threads; j++)
    {
      edges.at(current_level) += partial_edgesExplored[j];
      edgesFrontier += partial_edgesFrontier[j];
    }
  delete[] partial_edgesExplored;
  delete[] partial_edgesFrontier;
  return nodes.at(current_level) < ((int)num_vertices / beta) ? Top_Down:Bottom_Up;
}
