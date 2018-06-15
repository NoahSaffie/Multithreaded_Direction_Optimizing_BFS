#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include <chrono>
#define alpha 12
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
  Investigate where checking (if we should switch) belongs  and how often it gets checked (Is this a every iteration thing? If so might want to better calculate m_f (edgesFrontier) as we go since otheriwse will be a costly calculation)
  Future Issues/Questions:
  How does the extra (empty) vertexes effect BFS (Other than larger memory footprint), if at all?
  Picking a good start (i.e 13 as a start for the CA- is basically just a very small graph, that is not conencted to any large one -- want to avoid this)
*/
enum NextToRun { Top_Down, Bottom_Up };
typedef graph<long, long, int, long, long, char> Graph;
inline NextToRun bfs_top_down(Graph *ginst, std::queue<int> &frontier, std::vector<bool> &bitmap, int &edgesExplored, int &edgesFrontier, unsigned num_vertices);
inline NextToRun bfs_bottom_up(Graph *ginst, std::queue<int> &frontier, std::vector<bool> &bitmap, int &edgesExploredm, int &edgesFrontier, unsigned num_vertices);
inline NextToRun bfs_switch_(Graph *ginst, std::queue<int> frontier, int edgesExplored, unsigned num_vertices);

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
  int edgesExplored = 0;
  int edgesFrontier = 0;
  std::vector<bool> bitmap(ginst->vert_count, 0); 
  std::queue<int> frontier;
  //Frontier - Holds Node/Vertex # which also act as an Index in beg_pos
  //Init the frontier - Find the lowest Node/Vertex that has an adjacency
  for(int i = 0; i < ginst->vert_count; i++)
    {
      //Once we have found the index in beg_pos that is no longer looking at the first index of edge list (->csr), then we are currently sitting on the true index that the adjency list
      //begining at ->csr[0] belongs to
      if(ginst->beg_pos[i+1] != 0) //We skip any "missing" vertexs (no edges) by skipping all but the last vertex that has beg_pos[i] of 0
        {
          frontier.push(i);
          bitmap.at(i) = true;
          break;
        }
    }
  //std::cout << "Node: " << frontier.front() << std::endl;
  edgesFrontier += (ginst->beg_pos[frontier.front()+1] - ginst->beg_pos[frontier.front()]);
  
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
  int edges_TD = 0;
  int edges_BU = 0;
  int edges_base = 0;
  //std::cout << "True number of Vertices: " << num_vertices << std::endl;
  double time_TD = 0.0;
  double time_BU = 0.0;
  TimeVar start = timeNow();
  int totalSwitches = 0;
  NextToRun next = Top_Down;
  while (!frontier.empty())
    {
      // Will switch to the other when the previous one returns (which it will do depending on the value it recieves from switch function)
      // Requires one check extra check somewhere if we finish program on top_down, we need to know to skip bottom_up and break
      if(next == Top_Down)
	{
	  TimeVar start_TD = timeNow();
	  next = bfs_top_down(ginst, frontier, bitmap, edgesExplored, edgesFrontier, num_vertices);
	  edges_TD += (edgesExplored-edges_base);
	  edges_base = edgesExplored;
	  TimeVar end_TD = timeNow();
	  time_TD += duration(end_TD-start_TD);
	  ++totalSwitches;
	}
      if(frontier.empty()) { std::cout << "End of BFS" << std::endl; break; }
      //std::cout << "Switching to Bottom-Up\n" << std::endl;
      if(next == Bottom_Up)
	{
	  TimeVar start_BU = timeNow();
	  next = bfs_bottom_up(ginst, frontier, bitmap, edgesExplored, edgesFrontier, num_vertices);
	  edges_BU += (edgesExplored-edges_base);
	  edges_base = edgesExplored;
	  TimeVar end_BU = timeNow();
	  time_BU += duration(end_BU-start_BU);
	  ++totalSwitches;
	}
      //std::cout << "Switching to Top-Down\n" << std::endl;
    }
  TimeVar end = timeNow();
  double time_total = duration(end-start);
  int nodesVisited = 0;
  for(unsigned i = 0; i < bitmap.size(); i++)
    {
      if(bitmap.at(i)) { ++nodesVisited; }
    }
  std::cout << "Switches: " << totalSwitches << "\tNodes Visited: " << nodesVisited << std::endl;
  std::cout << "Edges: " << "\t\t%Top_Down: " << 100*(edges_TD*1.0/edgesExplored) << "\t%Bottom_Up: " << 100*(edges_BU*1.0/edgesExplored) << "\tTotal: " << edgesExplored << std::endl;
  std::cout << "Timing(ms): " << "\t%Top_Down: " << 100*(time_TD/time_total) << "\t%Bottom_Up: " << 100*(time_BU/time_total) << "\tTotal: " << (int)time_total << std::endl;
  
  // Cleanup
  free(ginst->beg_pos);
  free(ginst->csr);
  free(ginst->weight);
  delete ginst;
  return 0;
}
inline NextToRun bfs_top_down(Graph *ginst, std::queue<int> &frontier, std::vector<bool> &bitmap, int &edgesExplored, int &edgesFrontier, unsigned num_vertices)
{ //Consider seeing if we can move OpenMP up a level -- Would cause an issue with the switch check
  while (!frontier.empty()) //Add a heuristic piece here as well?
    {
      //Check neighbors
      
      int num_threads = omp_get_num_threads(); //4 for this system {0,1,2,3}
      int* partial_edgesExplored = (int*)malloc(sizeof(int)*num_threads); assert(partial_edgesExplored != NULL);
      int* partial_edgesFrontier = (int*)malloc(sizeof(int)*num_threads); assert(partial_edgesFrontier != NULL);
      //Init
      memset(partial_edgesExplored, 0, sizeof(int)*num_threads);
      memset(partial_edgesFrontier, 0, sizeof(int)*num_threads);
      #pragma omp parallel for num_threads(num_threads) default(none) firstprivate(partial_edgesExplored, partial_edgesFrontier,ginst) shared(bitmap, frontier)
      for (int begin = ginst->beg_pos[frontier.front()]; begin < ginst->beg_pos[frontier.front()+1]; begin++)
        {
	  int tid = omp_get_thread_num();
	  int node = ginst->csr[begin]; // Part of adjacency list of node i
          if (!bitmap.at(node))
            {
              //Mark as visited
              bitmap.at(node) = true;
              //Add to queue
              frontier.push(node);
	      partial_edgesFrontier[tid] += (ginst->beg_pos[node+1] - ginst->beg_pos[node]);
	      //printf("Node: %d\n", node);
            }
          //Mark edge as explored
          ++(partial_edgesExplored[tid]);
        }
      
      #pragma omp parallel for reduction(+:edgesExplored) reduction(+:edgesFrontier)
      for(int j = 0; j < num_threads; j++)
	{
	  edgesExplored += partial_edgesExplored[j];
	  edgesFrontier += partial_edgesFrontier[j];
	}
      //Pop off queue (and update edgesFrontier)
      edgesFrontier -= (ginst->beg_pos[frontier.front()+1] - ginst->beg_pos[frontier.front()]);
      frontier.pop();
      if(frontier.size() > (num_vertices / beta)) //edgesFrontier > ((ginst->edge_count - edgesExplored) / alpha)) // Hueristic
	{
	  free(partial_edgesExplored);
	  free(partial_edgesFrontier);
	  return Bottom_Up;
	}
      free(partial_edgesExplored);
      free(partial_edgesFrontier);
    }
  return Top_Down;
}
inline NextToRun bfs_bottom_up(Graph *ginst, std::queue<int> &frontier, std::vector<bool> &bitmap, int &edgesExplored, int &edgesFrontier, unsigned num_vertices)
{
  /*
    We don't need the frontier for this function but we do need it when we transition back to Top-Down so we add the respective nodes to the frontier as we go.
    However we want to avoid telling the function that we have visited (found a connection) to the node because then it might try (and succeed) at finding conenctions to that new
    point, which can be an issue. Especially when parallelization happens since it could cause inconsistent results across different runs
  */
  while (!frontier.empty()) { frontier.pop(); } //I believe we want a cleared frontier, since we are now finding all the "children" a different way, and don't need current queue
  edgesFrontier = 0;
  int num_threads = omp_get_num_threads();
  int* partial_edgesExplored = (int*)malloc(sizeof(int)*num_threads);
  int* partial_edgesFrontier = (int*)malloc(sizeof(int)*num_threads);
  //Init
  memset(partial_edgesExplored, 0, sizeof(int)*num_threads);
  memset(partial_edgesFrontier, 0, sizeof(int)*num_threads);
  
  #pragma omp parallel for num_threads(num_threads) default(none) firstprivate(ginst) shared(partial_edgesExplored, partial_edgesFrontier, bitmap, frontier)
  for (int i = 0; i < ginst->vert_count; i++)
    {
      int tid = omp_get_thread_num(); //Ideally not have in loop at all
      if (!bitmap.at(i)) //Unexplored
        {
	  int end = ginst->beg_pos[i+1];
          for (int begin = ginst->beg_pos[i]; begin < end; begin++) //Indexes in CSR
            {
	      int node = ginst->csr[begin]; //Neighbor of i
              ++partial_edgesExplored[tid];
              if (bitmap.at(node)) //Find a neighbor that has been visited, it will be a parent (on frontier too)
                {
                  frontier.push(i);
		  //printf("Node: %d\n", i);
		  partial_edgesFrontier[tid] += (ginst->beg_pos[i+1] - ginst->beg_pos[i]);
                  break;
                }
            }
        }
    }
  #pragma omp parallel for reduction(+:edgesExplored) reduction(+:edgesFrontier)
  for(int j = 0; j < num_threads; j++)
    {
      edgesExplored += partial_edgesExplored[j];
      edgesFrontier += partial_edgesFrontier[j];
    }
  free(partial_edgesExplored);
  free(partial_edgesFrontier);

  std::queue<int> mark(frontier); //Mark frontier as visited
  while(!mark.empty()) //While not empty
    {
      bitmap.at(mark.front()) = true; //Mark the node
      mark.pop(); //Pop it off - Job is done
    }
  if(frontier.size() < (num_vertices / beta)) //Hueristic
    {
      return Top_Down;
    }
  else
    {
      return Bottom_Up;
    }
}
