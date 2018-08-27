#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include "wtime.h"

/*
  Author: Noah Saffie
  Purpose: Implementation of Direction-Optimizing BFS, including the use of OpenMP
  Date last modified: 8/27/2018

  Additional Notes: 
  -- The Top Down and Bottom Up BFS functions can be found at the very bottom. I moved parallelization to the highest level I could so those functions are very simple
  -- Some of the output information (node count per level) is off by a tiny bit when doing multithreaded. This is easily fixed by using a #pragma omp critical but it comes at a high performance cost
     so I was advised against including it.
  -- On latest testing with valgrind there were only two "memory leaks" (200 bytes total steak reachable), one is a OpenMP issue (192 bytes) that I can do nothing about, and the other had no tracing information available so was not able to identify the problem (8 bytes). There can exist more that show up based due to OpenMP at times as well, but it is unclear if anything is actually lost or if it is a issue between valgrind, and OpenMP. (I found on StackOverflow this is a common issue, that can not be avoided but is obviouisly very minor here)
 */
enum NextToRun { Top_Down, Bottom_Up };
typedef graph<long, long, int, long, long, char> Graph;
inline void bfs_top_down(Graph *ginst, int* level, int* results, int current_level, const int section_beg, const int section_end);
inline void bfs_bottom_up(Graph *ginst, int* level, int* results, int current_level, const int section_beg, const int section_end);

int main(int argc, char* argv[])
{
  //Required for thread cancellation
  if(!omp_get_cancellation())
    {
      char *env = (char*)"OMP_CANCELLATION=true";
      printf("Cancellations were not enabled, enabling cancellation and rerunning program\n");
      putenv(env);
      execv(argv[0], argv);
    }
  //Graph Init
  std::cout << "Input: ./exe beg csr weight number_threads root(-1 if not sure) alpha(-1 for default) beta(-1 for default) tests_run test_type(0 - Direction-Optimizing, 1 - Just Top Down, 2 - Just Bottom Up\n";
  if (argc != 10) { std::cout << "Wrong input\n"; return -1; }
  const char *beg_file = argv[1]; //Vertices -- Like RowPtr
  const char *csr_file = argv[2]; //Edges -- That lead to the Node described from the Node (index+1) in Offsets array (Beg_file)
  const char *weight_file = argv[3]; //Weight -- Unused
  const int number_threads = atoi(argv[4]);
  Graph *ginst = new Graph(beg_file, csr_file, weight_file);
  int root = atoi(argv[5]);
  double alpha = .07142; //Default values based on the source paper (Scott Beamer "Direction Optimizing BFS")
  double beta = .04166;
  //If user provided a value, set and use that
  if(atof(argv[6]) == -1)
    {
      alpha = atof(argv[6]);
    }
  if(atof(argv[7]) == -1)
    {
      beta = atof(argv[7]);
    }
  int style_to_run = atoi(argv[9]);
  if(style_to_run > 2 || style_to_run < 0)
    {
      std::cout << "Wrong input\n"; return -1;
    }
  int tests_run = atoi(argv[8]);
  int copy_tests_run = tests_run; //Used for calculating the average time at the end
  if(tests_run < 0){ std::cout << "Wrong input\n"; return -1; }
  //General Init - Not a great way of dealing with this initialization, but put this way for ease when tests_run > 1
  std::vector<int> init_nodes;
  std::vector<int> init_edges;
  init_nodes.push_back(0); init_edges.push_back(0); //First Level with only root
  unsigned init_num_vertices = 1;
  int init_edgesExplored = 0, init_edgesFrontier = 0;
  //If the user did not provide a root, we will look for the first valid node number and use that as the root
  if(root == -1)
    {
      for(int i = 0; i < ginst->vert_count; i++)
        {   //Find first non-zero value in beg_pos, "root" is the index before that
          if(ginst->beg_pos[i+1] != 0) //We skip any "missing" vertexs (no edges) by skipping all but the last vertex that has beg_pos[i] of 0
            { root = i; break; }
        }
    }
  printf("Root: %d\tThis is included in level 0\n", root);
  init_nodes.at(0) = 1;
  init_edges.at(0) = (ginst->beg_pos[root+1] - ginst->beg_pos[root]);
  init_edgesFrontier += init_edges.at(0);

  //Find true number of verticies -- Number of verticles provided by ginst is incorrect since it pads in vertices that were not declared between 0 and the highest vertex id it found
  int last_num = 0;
  for(int j = 1; j < ginst->vert_count; ++j)
    {
      if(ginst->beg_pos[j] > last_num) { last_num = ginst->beg_pos[j]; ++init_num_vertices; }
    }
  //Stats - Init
  double time_total = 0; //Time tracker - double should be plenty of space for almost all parameters (unless of course tests_run is needlessly high AND provided dataset is abnormally large)
  //While we have more tests to run
  while(tests_run > 0)
    {
      //Reset info for new test run
      int num_vertices = init_num_vertices;
      std::vector<int> nodes = init_nodes;
      std::vector<int> edges = init_edges;
      int* level = new int[ginst->vert_count];
      memset(level, -1, sizeof(int)*ginst->vert_count);
      level[root] = 0;
      int edgesFrontier = init_edgesFrontier;
      int edgesExplored = init_edgesExplored;

      double start = wtime();
      NextToRun next = Top_Down;
      int current_level = 0;
      //Main Loop
      //Problem -- Knowing what cache line size line is per user? (Mine is 64 bytes will just work with that)
      //Note: Padding seemed to provide no increase in performance, however it had no noticeable decrease either, so it will be left in
      int cache_line_size = 64; //Bytes
      int size_int = sizeof(int); // Should be 4 here
      int number_stats = 3;
      int padding = cache_line_size/size_int + 1; //Pad away one full cache size;
      int spacing = cache_line_size/(size_int*number_stats);
      int extra_padding = cache_line_size%(size_int*number_stats);
      //Necessary to prevent thread overlap when producing values - 
      int* results = new int[padding + (number_threads*(number_stats*spacing+extra_padding))];
      memset(results, 0, sizeof(int)*(padding + (number_threads*(number_stats*spacing+extra_padding)))); //Results: [nodes, edgesExplored, edgesFrontier] repeated for number_threads
      //Just bottom_Up
      if(style_to_run == 2)
        {
          #pragma omp parallel num_threads(number_threads) default(none) firstprivate(ginst, alpha, beta, number_threads, num_vertices, padding, spacing, number_stats, extra_padding) shared(level, current_level, edges, nodes, edgesFrontier, edgesExplored, next, results)
          {
            const int tid = omp_get_thread_num();
            const int section_size = ginst->vert_count/number_threads;
            const int section_beg = tid*section_size;
            const int section_end = (tid==number_threads-1 ? ginst->vert_count:section_beg+section_size);
            while (current_level == 0 || nodes.at(current_level-1) > 0)
              {
                #pragma omp barrier
                bfs_bottom_up(ginst, level, &(results[padding + (tid*(number_stats*spacing+extra_padding))]), current_level, section_beg, section_end);
                #pragma omp barrier
                #pragma omp single
                {
                  for(int i = 0; i < number_threads; i++)
                    {
                      nodes.at(current_level) += results[padding + (i*(number_stats*spacing+extra_padding))];
                      edges.at(current_level) += results[padding + (i*(number_stats*spacing+extra_padding))+1];
                      edgesFrontier += results[padding + (i*(number_stats*spacing+extra_padding))+2];
                    }
                  edgesExplored += edges.at(current_level);
                  next = Bottom_Up;
                  current_level++; nodes.push_back(0); edges.push_back(0);
                  memset(results, 0, sizeof(int)*(padding + (number_threads*(number_stats*spacing+extra_padding))));
                }
              }
          }
        }
      //Just TOP_DOWN
      if(style_to_run == 1)
        {
          #pragma omp parallel num_threads(number_threads) default(none) firstprivate(ginst, alpha, beta, number_threads, num_vertices, padding, spacing, number_stats, extra_padding) shared(level, current_level, edges, nodes, edgesFrontier, edgesExplored, next, results)
          {
            const int tid = omp_get_thread_num();
            const int section_size = ginst->vert_count/number_threads;
            const int section_beg = tid*section_size;
            const int section_end = (tid==number_threads-1 ? ginst->vert_count:section_beg+section_size);
            while(current_level == 0 || nodes.at(current_level-1) > 0)
              {
                #pragma omp barrier
                bfs_top_down(ginst, level, &(results[padding + (tid*(number_stats*spacing+extra_padding))]), current_level, section_beg, section_end);
                #pragma omp barrier
                #pragma omp single
                {
                  for(int i = 0; i < number_threads; i++)
                    {
                      nodes.at(current_level) += results[padding + (i*(number_stats*spacing+extra_padding))];
                      edges.at(current_level) += results[padding + (i*(number_stats*spacing+extra_padding))+1];
                      edgesFrontier += results[padding + (i*(number_stats*spacing+extra_padding))+2];
                    }
                  edgesExplored += edges.at(current_level);
                  next = (edgesFrontier/(ginst->edge_count-edgesExplored) > alpha) ? Bottom_Up:Top_Down;
                  current_level++; nodes.push_back(0); edges.push_back(0);
                  memset(results, 0, sizeof(int)*(padding + (number_threads*(number_stats*spacing+extra_padding))));
                }
              }
            #pragma omp cancel parallel
          }
        }
      if(style_to_run == 0)
        {
          #pragma omp parallel num_threads(number_threads) default(none) firstprivate(ginst, alpha, beta, number_threads, num_vertices, padding, spacing, number_stats, extra_padding) shared(level, current_level, edges, nodes, edgesFrontier, edgesExplored, next, results)
          {
            const int tid = omp_get_thread_num();
            const int section_size = ginst->vert_count/number_threads;
            const int section_beg = tid*section_size;
            const int section_end = (tid==number_threads-1 ? ginst->vert_count:section_beg+section_size);
            while (current_level == 0 || nodes.at(current_level-1) > 0)
              {
                if(next == Top_Down)
                  {
                    #pragma omp barrier
                    bfs_top_down(ginst, level, &(results[padding + (tid*(number_stats*spacing+extra_padding))]), current_level, section_beg, section_end);
                    #pragma omp barrier
                    #pragma omp single
                    {
                      for(int i = 0; i < number_threads; i++)
                        {
                          nodes.at(current_level) += results[padding + (i*(number_stats*spacing+extra_padding))];
                          edges.at(current_level) += results[padding + (i*(number_stats*spacing+extra_padding))+1];
                          edgesFrontier += results[padding + (i*(number_stats*spacing+extra_padding))+2];
                        }
                      edgesExplored += edges.at(current_level);
                      next = ((1.0)*edgesFrontier/(ginst->edge_count-edgesExplored) > alpha) ? Bottom_Up:Top_Down;
                      current_level++; nodes.push_back(0); edges.push_back(0);
                      memset(results, 0, sizeof(int)*(padding + (number_threads*(number_stats*spacing+extra_padding))));
                    }
                  }
                if(next == Bottom_Up)
                  {
                    #pragma omp barrier
                    bfs_bottom_up(ginst, level, &(results[padding + (tid*(number_stats*spacing+extra_padding))]), current_level, section_beg, section_end);
                    #pragma omp barrier
                    #pragma omp single
                    {
                      for(int i = 0; i < number_threads; i++)
                        {
                          nodes.at(current_level) += results[padding + (i*(number_stats*spacing+extra_padding))];
                          edges.at(current_level) += results[padding + (i*(number_stats*spacing+extra_padding))+1];
                          edgesFrontier += results[padding + (i*(number_stats*spacing+extra_padding))+2];
                        }
                      edgesExplored += edges.at(current_level);
                      next = (1.0)*nodes.at(current_level)/num_vertices <  beta ? Top_Down:Bottom_Up;
                      current_level++; nodes.push_back(0); edges.push_back(0);
                      memset(results, 0, sizeof(int)*(padding + (number_threads*(number_stats*spacing+extra_padding))));
                    }
                  }
              }
            #pragma omp cancel parallel
          }
        }
      double end = wtime();
      time_total+= (end-start);
      //Only need to output data on last run
      if(tests_run == 1)
        {
          for(unsigned i = 0; i < nodes.size(); ++i)
            {
              if(nodes.at(i) != 0)
                {
                  printf("Level: %d\tFrontier Size: %d\tEdges: %d\n", i, nodes.at(i), edges.at(i));
                }
            }
        }
      tests_run--;
      delete[] results;
      delete[] level;
    }
  std::cout << "Timing(seconds): " << time_total/((1.0)*copy_tests_run) << std::endl;

  // Cleanup
  free(ginst->beg_pos);
  free(ginst->csr);
  free(ginst->weight);
  delete ginst;
  return 0;
}

inline void bfs_top_down(Graph *ginst, int* level, int* results, int current_level, const int section_beg, const int section_end)
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
              ++results[1];
              if(level[node] == -1) //If neighbor is unexplored
                {
                  level[node] = current_level+1; //Add to next frontier
                  ++results[0];
                  results[2] += (ginst->beg_pos[node+1] - ginst->beg_pos[node]);
                }
            }
        }
    }
}
inline void bfs_bottom_up(Graph *ginst, int* level, int* results, int current_level, const int section_beg, const int section_end)
{
  for (int vert_index = section_beg; vert_index < section_end; ++vert_index)
    {
      if (level[vert_index]==-1) //Unexplored (Only want to look through potential children of frontier)
        {
          long int neighbor_beg = ginst->beg_pos[vert_index];
          long int neighbor_end = ginst->beg_pos[vert_index+1];
          for (long int neighbor_index = neighbor_beg; neighbor_index < neighbor_end; ++neighbor_index) //Indexes in CSR
            {
              long int node = ginst->csr[neighbor_index]; //Neighbor of vert_index
              ++results[1];
              if (level[node] == current_level) //Visited neighbor limited to current_level
                {
                  level[vert_index] = current_level+1;
                  ++results[0]; //Should be 1 for next_level and 0 for previous level (==current_level)
                  results[2] += (ginst->beg_pos[vert_index+1] - ginst->beg_pos[vert_index]);
                  break;
                }
            }
        }
    }
}
