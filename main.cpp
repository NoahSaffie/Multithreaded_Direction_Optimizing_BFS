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
enum NextToRun { Top_Down, Bottom_Up };
typedef graph<int, long, int, long, long, char> Graph;
/* Node # starts at 1 version */
inline void bfs_top_down(Graph *ginst, std::queue<int> frontier, std::vector<bool> bitmap, int &exploredEdges);
inline void bfs_bottom_up(Graph *ginst, std::queue<int> &frontier, std::vector<bool> bitmap, int &exploredEdges);
inline NextToRun bfs_switch(Graph *ginst, std::queue<int> frontier, int exploredEdges);

/* Node # starts at 0 version */
inline void bfs_top_down_0(Graph *ginst, std::queue<int> frontier, std::vector<bool> bitmap, int &exploredEdges);
inline void bfs_bottom_up_0(Graph *ginst, std::queue<int> &frontier, std::vector<bool> bitmap, int &exploredEdges);
inline NextToRun bfs_switch_0(Graph *ginst, std::queue<int> frontier, int exploredEdges);

int main(int argc, char* argv[])
{ 
	/* 
	Assumes that nodes/vertices start at 1 and go from there (based on example file) 
	May need to Modify so it can work either way (Check if first node is marked as 0 or 1 and proceed based on result)
	*/
  std::cout << "Input: ./exe beg csr weight\n";
  if (argc != 4) { std::cout << "Wrong input\n"; return -1; }

  const char *beg_file = argv[1]; //Vertices -- Like RowPtr
  const char *csr_file = argv[2]; //Edges -- That lead to the Node descrived from the Node (index+1) in Offsets array (Beg_file)
  const char *weight_file = argv[3]; //Weight
  //template <file_vertex_t, file_index_t, file_weight_t
  //new_vertex_t, new_index_t, new_weight_t>
  Graph *ginst = new Graph(beg_file, csr_file, weight_file);
  // Testing
  for(int i = 0; i < 40+1; i++)
    {
      int beg = ginst->beg_pos[i];
      int end = ginst->beg_pos[i+1];
      std::cout<<i<<"'s neighor list: ";
      for(int j = beg; j < end; j++)
	std::cout<<ginst->csr[j]<<" ";
      std::cout<<"\n";
    } 
  //std::bitset<vertices> bitmap; //Marking from Index 1 -> end -- Can't use  due to varying size
  std::vector<bool> bitmap(ginst->vert_count+1); //
  bitmap.at(0) = true;
  std::queue<int> frontier;
  int exploredEdges = 0;
  /*
  while (true) //---- TO BE CHANGED ----
  {
	 // Will switch to the other when the previous one returns (which it will do depending on the value it recieves from switch function)
	 // Requires one check extra check somewhere if we finish program on top_down, we need to know to skip bottom_up
	  bfs_top_down(ginst, frontier, bitmap, exploredEdges);

	  bfs_bottom_up(ginst, frontier, bitmap, exploredEdges);
  }
  */
  // Cleanup
  free(ginst->beg_pos);
  free(ginst->csr);
  free(ginst->weight);
  delete ginst;
  return 0;
}
/*


 */
inline void bfs_top_down(Graph *ginst, std::queue<int> frontier, std::vector<bool> bitmap, int &exploredEdges)
{
	/* Needs Heruistic Added and potentially weight info */
	//Add intial
	frontier.push(1);
	//Mark it as visited
	bitmap.at(1) = true;
	while (!frontier.empty()) //Add a heuristic piece here as well?
	{
		//Check neighbors
		for (int begin = ginst->beg_pos[frontier.front()-1]; begin < ginst->beg_pos[frontier.front()]; begin++)
		{
			if (!bitmap.at(begin))
			{
				//Mark as visited
				bitmap.at(begin) = true;
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
inline void bfs_bottom_up(Graph *ginst, std::queue<int> &frontier, std::vector<bool> bitmap, int &exploredEdges)
{
	while (!frontier.empty()) { frontier.pop(); } //I believe we want a cleared frontier, since we are now finding all the "children" a different way, and don't need current queue
	for (int i = 1; i < ginst->vert_count; i++) //Node 1 on
	{
		if (!bitmap.at(i)) //Unexplored
		{
			for (int begin = ginst->beg_pos[i-1]; begin < ginst->beg_pos[i]; begin++)
			{
				//Mark edge as explored
				++exploredEdges;
				if (bitmap.at(begin))
				{
					frontier.push(i);
					bitmap.at(i) = true;
					break;
				}
			}
		}
	}
}
inline NextToRun bfs_switch(Graph *ginst, std::queue<int> frontier, int exploredEdges)
{
	std::queue<int> copy_frontier(frontier); //Seems like a big downside -- Creating a copy
	int edgesFrontier = 0; //m_f -- For each Vertice in Queue, find beg_pos[vertex+1]-beg_pos[vertex]
	int edgesUnexploredNodes = ginst->edge_count - exploredEdges; //m_u -- 
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


/*
  This is the starts at Vertex of 0 versions
  Changes:
      Intial push is 0 instead of 1
      Loop conditions are, i, and i+1 vs i-1, i
      Loop init is 0 not 1
 */
inline void bfs_top_down_0(Graph *ginst, std::queue<int> frontier, std::vector<bool> bitmap, int &exploredEdges)
{
	/* Needs Heruistic Added and potentially weight info */
	//Add intial
	frontier.push(0);
	//Mark it as visited
	bitmap.at(0) = true;
	while (!frontier.empty()) //Add a heuristic piece here as well?
	{
		//Check neighbors
		for (int begin = ginst->beg_pos[frontier.front()]; begin < ginst->beg_pos[frontier.front()+1]; begin++)
		{
			if (!bitmap.at(begin))
			{
				//Mark as visited
				bitmap.at(begin) = true;
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
inline void bfs_bottom_up_0(Graph *ginst, std::queue<int> &frontier, std::vector<bool> bitmap, int &exploredEdges)
{
	while (!frontier.empty()) { frontier.pop(); } //I believe we want a cleared frontier, since we are now finding all the "children" a different way, and don't need current queue
	for (int i = 0; i < ginst->vert_count; i++) //Node 1 on
	{
		if (!bitmap.at(i)) //Unexplored
		{
			for (int begin = ginst->beg_pos[i]; begin < ginst->beg_pos[i+1]; begin++)
			{
				//Mark edge as explored
				++exploredEdges;
				if (bitmap.at(begin))
				{
					frontier.push(i);
					bitmap.at(i) = true;
					break;
				}
			}
		}
	}
}

inline NextToRun bfs_switch_0(Graph *ginst, std::queue<int> frontier, int exploredEdges)
{
	std::queue<int> copy_frontier(frontier); //Seems like a big downside -- Creating a copy
	int edgesFrontier = 0; //m_f -- For each Vertice in Queue, find beg_pos[vertex+1]-beg_pos[vertex]
	int edgesUnexploredNodes = ginst->edge_count - exploredEdges; //m_u -- 
	int verticesFrontier = copy_frontier.size(); //n_f -- Size of Queue
	while (!copy_frontier.empty())
	{
		edgesFrontier += (ginst->beg_pos[copy_frontier.front()+1] - ginst->beg_pos[copy_frontier.front()]);
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
