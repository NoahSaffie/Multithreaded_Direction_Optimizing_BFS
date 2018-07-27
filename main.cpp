#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include "wtime.h"

#define just 0
enum NextToRun { Top_Down, Bottom_Up };
typedef graph<long, long, int, long, long, char> Graph;
inline void bfs_top_down(Graph *ginst, int* level, int* results, int current_level, const int section_beg, const int section_end);
inline void bfs_bottom_up(Graph *ginst, int* level, int* results, int current_level, const int section_beg, const int section_end);

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
  std::cout << "Input: ./exe beg csr weight number_threads root(-1 if not sure) (optional)[alpha, beta]\n";
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
  printf("Root: %d\tThis is included in level 0\n", root);
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
  //Main Loop
  //Problem -- Knowing what cache line size line is per user? (Mine is 64 bytes will just work with that)
  int cache_line_size = 64; //Bytes
  int size_int = sizeof(int); // Should be 4 here
  int number_stats = 3;
  int padding = cache_line_size/size_int + 1; //Pad away one full cache size;
  int spacing = cache_line_size/(size_int*number_stats);
  int extra_padding = cache_line_size%(size_int*number_stats);
  int* results = new int[padding + (number_threads*(number_stats*spacing+extra_padding))];
  memset(results, 0, sizeof(int)*(padding + (number_threads*(number_stats*spacing+extra_padding)))); //Results: nodes, edgesExplored, edgesFrontier *number_threads
  if(just == 2)
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
  if(just == 1)
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
  if(just == 0)
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
  delete[] results;
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
