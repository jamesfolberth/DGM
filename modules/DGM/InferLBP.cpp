#include "InferLBP.h"
#include "GraphPairwise.h"

namespace DirectGraphicalModels
{
	void CInferLBP::calculateMessages(unsigned int nIt)
	{
		const byte		nStates = getGraph().getNumStates();				// number of states
		
		// ======================== Main loop (iterative messages calculation) ========================
#ifndef ENABLE_PPL
		double *temp = new double[nStates];
#endif
		for (unsigned int i = 0; i < nIt; i++) {								// iterations
#ifdef DEBUG_PRINT_INFO
			if (i == 0) printf("\n");
			if (i % 5 == 0) printf("--- It: %d ---\n", i);
#endif
#ifdef ENABLE_PPL
			concurrency::parallel_for_each(getGraphPairwise().m_vNodes.begin(), getGraphPairwise().m_vNodes.end(), [&, nStates](ptr_node_t &node) {		// all nodes
				double *temp = new double[nStates];
#else
			std::for_each(getGraphPairwise().m_vNodes.begin(), getGraphPairwise().m_vNodes.end(), [&](ptr_node_t &node) {
#endif
				// Calculate a message to each neighbor
				for (size_t e_t : node->to) {									// outgoing edges
					Edge *edge_to = getGraphPairwise().m_vEdges[e_t].get();		// current outgoing edge
					calculateMessage(*edge_to, temp, getMessageTemp(e_t), m_maxSum);
				} // e_t;
#ifdef ENABLE_PPL
				delete[] temp;
#endif
			}); // nodes
			swapMessages();														// Coping data from msg_temp to msg
		} // iterations
#ifndef ENABLE_PPL
		delete[] temp;
#endif
	}
}
