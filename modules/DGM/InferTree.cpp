#include "InferTree.h"
#include "GraphPairwise.h"

namespace DirectGraphicalModels
{
	void CInferTree::calculateMessages(unsigned int)
	{
		const byte		nStates	= getGraph().getNumStates();
		const size_t	nNodes	= getGraph().getNumNodes();
		const size_t	nEdges	= getGraph().getNumEdges();

		// ====================================== Initialization ======================================
		vec_bool_t		isReady(nEdges, false);								// Flags indicating whether the messages were already calculated
		vec_bool_t		suspend(nEdges, false);								// Flags indicating weather the message calculation must be postponed
		
		// =================================== Computing messages ===================================
		size_t  * nFromEdges = new size_t[nNodes];							// Count number of neighbors
		std::deque<size_t> nodeQueue;
		for (ptr_node_t &node: getGraphPairwise().m_vNodes) {
			nFromEdges[node->id] = node->from.size();						// number of incoming edges
			if (nFromEdges[node->id] <= 1) nodeQueue.push_back(node->id);	// Add all leafs to the queue
		}


		double *temp = new double[nStates];
		while (!nodeQueue.empty()) {
			//for (size_t q = 0; q < nodeQueue.size(); q++) printf("%d, ", nodeQueue[q]);	printf("\n");
			
			size_t n = nodeQueue.front();									// n - node with one neighbour
			nodeQueue.pop_front();

			Node  *node = getGraphPairwise().m_vNodes[n].get();				// Node with one neighbour
			
			bool allSuspend = true;
			for (size_t e_t : node->to)
				if (!suspend[e_t]) {
					allSuspend = false;
					break;
				}

			if (allSuspend) {	// Now prepare messages for suspending edges
				for (size_t e_t : node->to) {
					if (isReady[e_t]) continue;
					
					Edge *edge_to = getGraphPairwise().m_vEdges[e_t].get();
					calculateMessage(*edge_to, temp, getMessage(e_t));
					isReady[e_t] = true;
					
					// ------
					size_t n1 = edge_to->node1;
					size_t n2 = edge_to->node2;
					auto it = std::find_if(getGraphPairwise().m_vNodes[n1]->from.begin(), getGraphPairwise().m_vNodes[n1]->from.end(), [&](size_t e) {
						return (getGraphPairwise().m_vEdges[e]->node1 == n2);
					});
					if (it != getGraphPairwise().m_vNodes[n1]->from.end())
						suspend[*it] = true;
					// ------
					
					nFromEdges[n2]--;
					if (nFromEdges[n2] <= 1) nodeQueue.push_back(n2);
				}
			} else {			// Prepare messages for all non-suspending edges
				for (size_t e_t : node->to) {
					if (suspend[e_t]) continue;
					if (isReady[e_t]) continue;
					
					Edge * edge_to = getGraphPairwise().m_vEdges[e_t].get();
					calculateMessage(*edge_to, temp, getMessage(e_t));
					isReady[e_t] = true;
					// ------
					size_t n1 = edge_to->node1;
					size_t n2 = edge_to->node2;
					auto it = std::find_if(getGraphPairwise().m_vNodes[n1]->from.begin(), getGraphPairwise().m_vNodes[n1]->from.end(), [&](size_t e) {
						return (getGraphPairwise().m_vEdges[e]->node1 == n2);
					});
					if (it != getGraphPairwise().m_vNodes[n1]->from.end())
						suspend[*it] = true;
					// ------
					
					nFromEdges[n2]--;
					if (nFromEdges[n2] <= 1) nodeQueue.push_back(n2);
				}
			}
		} // while

		delete[] temp;
		delete[] nFromEdges;
	}
}
