#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include <bitset>
#define alpha 12
#define beta 24
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
inline void bfs_top_down(Graph *ginst, std::queue<int> frontier, std::vector<bool> bitmap, int &exploredEdges, int &edgesFrontier);
inline void bfs_bottom_up(Graph *ginst, std::queue<int> &frontier, std::vector<bool> bitmap, int &exploredEdgesm, int &edgesFrontier, unsigned num_vertices);
inline NextToRun bfs_switch_(Graph *ginst, std::queue<int> frontier, int exploredEdges, unsigned num_vertices);

int main(int argc, char* argv[])
{
  std::cout << "Input: ./exe beg csr weight\n";
  if (argc != 4) { std::cout << "Wrong input\n"; return -1; }

  const char *beg_file = argv[1]; //Vertices -- Like RowPtr
  const char *csr_file = argv[2]; //Edges -- That lead to the Node described from the Node (index+1) in Offsets array (Beg_file)
  const char *weight_file = argv[3]; //Weight -- Unused
  Graph *ginst = new Graph(beg_file, csr_file, weight_file);
  // Testing
  /*
    for(int i = 0; i < 40+1; i++)
    {
    int beg = ginst->beg_pos[i];
    int end = ginst->beg_pos[i+1];
    std::cout<<i<<"'s neighor list: ";
    for(int j = beg; j < end; j++)
    {
    std::cout<<ginst->csr[j]<<" ";
    }
    std::cout<<"\n";
    }
  */
  std::vector<bool> bitmap(ginst->vert_count+1); //
  bitmap.at(0) = true;
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
  std::cout << "Node: " << frontier.front() << std::endl;
  int exploredEdges = 0;
  int edgesFrontier = 0;
  edgesFrontier += (ginst->beg_pos[frontier.front()+1] - ginst->beg_pos[frontier.front()]);
  //Find true # of vertices
  unsigned num_vertices = 1;
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
  //std::cout << "True number of Vertices: " << num_vertices << std::endl;
  int totalSwitches = 0;
  while (!frontier.empty()) //---- TO BE CHANGED ----
    {
      // Will switch to the other when the previous one returns (which it will do depending on the value it recieves from switch function)
      // Requires one check extra check somewhere if we finish program on top_down, we need to know to skip bottom_up and break
      bfs_top_down(ginst, frontier, bitmap, exploredEdges, edgesFrontier);
      if(!frontier.empty()) { std::cout << "End of BFS" << std::endl; break; }
      std::cout << "Switching to Bottom-Up\n" << std::endl; ++totalSwitches;
      bfs_bottom_up(ginst, frontier, bitmap, exploredEdges, edgesFrontier, num_vertices);
      std::cout << "Switching to Top-Down\n" << std::endl; ++totalSwitches;
    }
  std::cout << "Total Numebr of Switches: " << totalSwitches << std::endl;
  // Cleanup
  free(ginst->beg_pos);
  free(ginst->csr);
  free(ginst->weight);
  delete ginst;
  return 0;
}
inline void bfs_top_down(Graph *ginst, std::queue<int> frontier, std::vector<bool> bitmap, int &exploredEdges, int &edgesFrontier)
{
  while (!frontier.empty()) //Add a heuristic piece here as well?
    {
      //Check neighbors
      for (int begin = ginst->beg_pos[frontier.front()]; begin < ginst->beg_pos[frontier.front()+1]; begin++)
        {
	  int node = ginst->csr[begin]; // Part of adjacency list of node i
          if (!bitmap.at(node))
            {
              //Mark as visited
              bitmap.at(node) = true;
              //Add to queue
              frontier.push(node);
	      edgesFrontier += (ginst->beg_pos[node+1] - ginst->beg_pos[node]);
	      //std::cout << "Node: " << node << std::endl;
            }
          //Mark edge as explored
          ++exploredEdges;
	  //std::cout << "edgesExplored: " << exploredEdges << std::endl;
        }
      //Pop off queue (and update edgesFrontier)
      edgesFrontier -= (ginst->beg_pos[frontier.front()+1] - ginst->beg_pos[frontier.front()]);
      frontier.pop();
      if(edgesFrontier > ((ginst->edge_count - exploredEdges) / alpha)) // Hueristic
	{
	  return;
	}
    }
}
inline void bfs_bottom_up(Graph *ginst, std::queue<int> &frontier, std::vector<bool> bitmap, int &exploredEdges, int &edgesFrontier, unsigned num_vertices)
{
  while (!frontier.empty()) { frontier.pop(); } //I believe we want a cleared frontier, since we are now finding all the "children" a different way, and don't need current queue
  for (int i = 0; i < ginst->vert_count; i++) //Node 1 on
    {
      if (!bitmap.at(i)) //Unexplored
        {
          for (int begin = ginst->beg_pos[i]; begin < ginst->beg_pos[i+1]; begin++)
            {
	      int node = ginst->csr[begin]; //Part of adjacency list of node i
              //Mark edge as explored
              ++exploredEdges;
              if (bitmap.at(node))
                {
                  frontier.push(i);
		  edgesFrontier += (ginst->beg_pos[i+1] - ginst->beg_pos[i]);
                  bitmap.at(i) = true;
                  break;
                }
            }
        }
      if(frontier.size() < (num_vertices / beta)) //Hueristic
	{
	  return;
	}
    }
}
//Obsolete -- Moved the checks directly into the functions
inline NextToRun bfs_switch(Graph *ginst, std::queue<int> frontier, int exploredEdges, int num_vertices)
{
  std::queue<int> copy_frontier(frontier); //Seems like a big downside -- Creating a copy
  int edgesFrontier = 0; //m_f -- For each Vertice in Queue, find beg_pos[vertex+1]-beg_pos[vertex]
  int edgesUnexploredNodes = ginst->edge_count - exploredEdges; //m_u --
  int verticesFrontier = copy_frontier.size(); //n_f -- Size of Queue
  while (!copy_frontier.empty())
    {
      edgesFrontier += (ginst->beg_pos[copy_frontier.front()+1] - ginst->beg_pos[copy_frontier.front()]); //Sum of degrees
      copy_frontier.pop();
    }
  if (edgesFrontier > (edgesUnexploredNodes / alpha)) // Go to Bottom-Up
    {
      return Bottom_Up;
    }
  else if (verticesFrontier < (num_vertices / beta)) //Go to Top-Down
    {
      return Top_Down;
    }
}
