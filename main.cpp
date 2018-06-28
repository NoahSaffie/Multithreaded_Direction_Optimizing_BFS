#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include "wtime.h"

enum NextToRun { Top_Down, Bottom_Up };
typedef graph<long, long, int, long, long, char> Graph;
inline NextToRun bfs_top_down(Graph *ginst, int &edgesFrontier, int &edgesExplored, int* level, std::vector<int> &edges, std::vector<int> &nodes, int current_level, double alpha, const int section_beg, const int section_end);
inline NextToRun bfs_bottom_up(Graph *ginst, int &edgesFrontier, int &edgesExplored, unsigned num_vertices, int* level, std::vector<int> &edges, std::vector<int> &nodes, int current_level, double beta, const int section_beg, const int section_end);

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
      for(int i = 0; i < ginst->vert_count; i++)
        {   //Find first non-zero value in beg_pos, "root" is the index before that
          if(ginst->beg_pos[i+1] != 0) //We skip any "missing" vertexs (no edges) by skipping all but the last vertex that has beg_pos[i] of 0
            { level[i] = 0; root = i; break; }
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
      if(ginst->beg_pos[j] > last_num) { last_num = ginst->beg_pos[j]; ++num_vertices; }
    }
  //Stats - Init
  double start = wtime();
  NextToRun next = Top_Down;
  int current_level = 0;
  unsigned total_nodes = 0;
  //Main Loop

  #pragma omp parallel num_threads(number_threads) default(none) firstprivate(ginst, alpha, beta, number_threads, total_nodes, num_vertices) shared(level, current_level, edges, nodes, edgesFrontier, edgesExplored, next)
  {
    const int tid = omp_get_thread_num();
    const int section_size = ginst->vert_count/number_threads;
    const int section_beg = tid*section_size;
    const int section_end = (tid==number_threads-1 ? ginst->vert_count:section_beg+section_size);
    while (num_vertices > total_nodes)
      { //Problem 1: Now that the parallel was moved out, each thread is going to call each function just for itself. Which creates major issues with keeping them "together" and with the same knowledge
        //Even with criticals, one thread could output its results and then go on to the heuristic before any other thread even got to the outputting of results. Which means, that some threads
        //Will begin to do different paths since one thread may see lower value when looking at the heuristic than other will, with potential for different returns
        //"Solution": Add more barriers in this main loop. Technically the last thread out should update "next" to be correct since it should have all the info

        //Problem 2: Top Down is gonna have issues though since it doesn't have to return after each level at the moment. This means the solution to Porblem 1 doesn't help
        //with Top Down. Another major issue here is that current_level is impossible to keep accurate, without breaking out of function after every level check
        //Since different threads could each increment the current_level, making it not only way to high, but more importantly creating an issue, where a thread that is still running,
        //now just had an updated current_level and might start marking levels wrong now. (This could be solved with doing a local copy of current_level, however other issues still remain).
        //Also creates the issue of the threads branching again. This time one thread may start on another level because it saw the heuristic with a lower value than the others.
        //"Solution": Break out of function after every single level is explored.
        //--- While main doesn't have a large amount of costs, this still adds more time ---
        if(next == Top_Down)
          {
            #pragma omp barrier
            next = bfs_top_down(ginst, edgesFrontier, edgesExplored, level, edges, nodes, current_level, alpha, section_beg, section_end);
            #pragma omp barrier
            #pragma omp single
            {
              current_level++; nodes.push_back(0); edges.push_back(0);
              printf("Exited TD\n");
            }
          }
        if(next == Bottom_Up)
          {
            #pragma omp barrier
            next = bfs_bottom_up(ginst, edgesFrontier, edgesExplored, num_vertices, level, edges, nodes, current_level, beta, section_beg, section_end);
            #pragma omp barrier
            #pragma omp single
            {
              current_level++; nodes.push_back(0); edges.push_back(0);
              printf("Exited BU\n");
            }
          }
        #pragma omp single
        {
          total_nodes = 0;
          for(unsigned i = 0; i < nodes.size(); ++i) { total_nodes+= nodes.at(i); } //printf("Total nodes: %d\n", total_nodes);
        }
      }
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
  std::cout << "Timing(seconds): " << time_total << std::endl;

  // Cleanup
  free(ginst->beg_pos);
  free(ginst->csr);
  free(ginst->weight);
  delete ginst;
  delete[] level;
  return 0;
}

inline NextToRun bfs_top_down(Graph *ginst, int &edgesFrontier, int &edgesExplored, int* level, std::vector<int> &edges, std::vector<int> &nodes, int current_level, double alpha, const int section_beg, const int section_end)
{
  int partial_nodes = 0;
  int partial_edgesExplored = 0;
  int partial_edgesFrontier = 0;
  for (int vert_index = section_beg; vert_index < section_end; vert_index++) //Go through your section
    {
        if (level[vert_index]==current_level) // If unexplored
          {
            int neighbor_beg = ginst->beg_pos[vert_index];
            int neighbor_end = ginst->beg_pos[vert_index+1];
            for(int neighbor_index = neighbor_beg; neighbor_index < neighbor_end; neighbor_index++) //Go through neighbors
              {
                int node = ginst->csr[neighbor_index];
                partial_edgesFrontier += (ginst->beg_pos[node+1] - ginst->beg_pos[node]);
                ++partial_edgesExplored;
                #pragma omp critical(NeighborUpdate)
                if(level[node] == -1) //If neighbor is unexplored
                  {
                    level[node] = current_level+1; //Add to next frontier
                    ++partial_nodes;
                  }
              }
	  }
    }
  #pragma omp critical(ResultsUpdate) //Stats updating and level update
  {
    edges.at(current_level) += partial_edgesExplored;
    edgesFrontier += partial_edgesFrontier;
    nodes.at(current_level) += partial_nodes;
    edgesExplored += edges.at(current_level);
  }
  return (edgesFrontier/(ginst->edge_count-edgesExplored) > alpha) ? Bottom_Up:Top_Down;
}
inline NextToRun bfs_bottom_up(Graph *ginst, int &edgesFrontier, int& edgesExplored, unsigned num_vertices, int* level, std::vector<int> &edges, std::vector<int> &nodes, int current_level, double beta, const int section_beg, const int section_end)
{
  int partial_nodes = 0;
  int partial_edgesExplored = 0;
  int partial_edgesFrontier = 0;
  edgesFrontier = 0;
  for (int vert_index = section_beg; vert_index < section_end; ++vert_index)
    {
      if (level[vert_index]==-1) //Unexplored (Only want to look through potential children of frontier)
        {
          long int neighbor_beg = ginst->beg_pos[vert_index];
          long int neighbor_end = ginst->beg_pos[vert_index+1];
          for (long int neighbor_index = neighbor_beg; neighbor_index < neighbor_end; ++neighbor_index) //Indexes in CSR
            {
              long int node = ginst->csr[neighbor_index]; //Neighbor of vert_index
              ++partial_edgesExplored;
              if (level[node] == current_level) //Visited neighbor limited to current_level
                {
                  level[vert_index] = current_level+1;
                  ++partial_nodes; //Should be 1 for next_level and 0 for previous level (==current_level)
                  partial_edgesFrontier += (ginst->beg_pos[vert_index+1] - ginst->beg_pos[vert_index]);
                  break;
                }
            }
        }
    }
  #pragma omp critical(ResultsUpdate) //Required now that each thread is gonna run its own function
  {
    nodes.at(current_level) += partial_nodes;
    edges.at(current_level) += partial_edgesExplored;
    edgesFrontier += partial_edgesFrontier;
    edgesExplored += edges.at(current_level);
  }
  return nodes.at(current_level)/num_vertices <  beta ? Top_Down:Bottom_Up;
}
