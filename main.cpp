#include <cstdlib>
#include <omp.h>
#include "graph.h"
#include <queue>
#include <bitset>
#define alpha 12
#define beta 24

enum NextToRun { Top_Down, Bottom_Up };
typedef graph<int, long, int, long, long, char> Graph
inline void bfs_top_down(Graph *ginst, std::queue<int> frontier, std::bitset bitmap, int &exploredEdges);
inline void bfs_bottom_up(Graph *ginst, std::queue<int> &frontier, std::bitset bitmap, int &edgesExplored);
inline NextToRun bfs_switch(Graph *ginst, std::queue<int> frontier, int edgesExplored);
int main(int argc, char* argv[])
{ 
	/* 
	Assumes that nodes/vertices start at 1 and go from there (based on example file) 
	May need to Modify so it can work either way (Check if first node is marked as 0 or 1 and proceed based on result)
	*/
  std::cout << "Input: ./exe beg csr weight\n";
  if (args != 4) { std::cout << "Wrong input\n"; return -1; }

  const char *beg_file = argv[1]; //Vertices -- Like RowPtr
  const char *csr_file = argv[2]; //Edges -- That lead to the Node descrived from the Node (index+1) in Offsets array (Beg_file)
  const char *weight_file = argv[3]; //Weight
  //template <file_vertex_t, file_index_t, file_weight_t
  //new_vertex_t, new_index_t, new_weight_t>
  Graph *ginst = new Graph(beg_file, csr_file, weight_file);
  std::bitset<ginst->vert_count+1> bitmap; //Marking from Index 1 -> end
  bitmap.set(0);
  std::queue<int> frontier;
  int exploredEdges = 0;
  while (bitmap.count() == bitmap.size())
  {
	 // Will switch to the other when the previous one returns (which it will do depending on the value it recieves from switch function)
	 // Requires one check extra check somewhere if we finish program on top_down, we need to know to skip bottom_up
	  bfs_top_down(ginst, frontier, bitmap, exploredEdges);

	  bfs_bottom_up(ginst, frontier, bitmap, edgesExplored);
  }
    
  return 0;
}
inline void bfs_top_down(Graph *ginst, std::queue<int> frontier, std::bitset bitmap, int &exploredEdges)
{
	/* Needs Heruistic Added and potentially weight info */
	//Add intial
	frontier.push(1);
	//Mark it as visited
	bitmap.set(1);
	while (!frontier.empty()) //Add a heuristic piece here as well?
	{
		//Check neighbors
		for (int begin = ginst->beg_pos[frontier.front()-1]; begin < ginst->beg_pos[frontier.front()]; begin++)
		{
			if (!bitmap.test(begin))
			{
				//Mark as visited
				bitmap.set(begin);
				//Add to queue
				frontier.push(begin);
			}
			//Mark edge as explored
			++exploredEdges;
		}
		//Pop off queue
		frontier.pop();
	}
}
inline void bfs_bottom_up(Graph *ginst, std::queue<int> &frontier, std::bitset bitmap, int &edgesExplored)
{
	while (!frontier.empty()) { frontier.pop(); } //I believe we want a cleared frontier, since we are now finding all the "children" a different way, and don't need current queue
	for (int i = 1; i < ginst->vert_count; i++) //Node 1 on
	{
		if (!bitmap.test(i)) //Unexplored
		{
			for (int begin = ginst->beg_pos[i-1]; begin < ginst->beg_pos[i]; begin++)
			{
				//Mark edge as explored
				++exploredEdges;
				if (bitmap.test(begin))
				{
					frontier.push(i);
					bitmap.set(i);
					break;
				}
			}
		}
	}
}
inline NextToRun bfs_switch(Graph *ginst, std::queue<int> frontier, int edgesExplored)
{
	std::queue<int> copy_frontier(frontier); //Seems like a big downside -- Creating a copy
	int edgesFrontier = 0; //m_f -- For each Vertice in Queue, find beg_pos[vertex+1]-beg_pos[vertex]
	int edgesUnexploredNodes = ginst->edge_count - edgesExplored; //m_u -- 
	int verticesFrontier = copy_frontier.size(); //n_f -- Size of Queue
	while (!copy_frontier.empty())
	{
		edgesFrontier += (ginst->beg_pos[copy_frontier.front()] - ginst->beg_pos[copy_frontier.front()-1]);
		copy_frontier.pop();
	}
	if (edgesFrontier > (edgesUnexploredNodes / alpha)) // Go to Bottom-Up
	{
		return Bottom_Up;
	}
	else if (verticesFrontier < (ginst->vert_count / beta)) //Go to Top-Down
	{
		return Top_Down;
	}
}
