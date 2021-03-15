#include <cmath>
#include <limits>

#include "LogMessagePassing.h"
#include "GraphPairwise.h"
#include "macroses.h"

namespace DirectGraphicalModels
{
    void CLogMessagePassing::infer(unsigned int nIt)
	{
		const byte   nStates = getGraph().getNumStates();

		// ====================================== Initialization ======================================
        createMessages(std::log(1. / nStates));				// msg[] = 1 / nStates; msg_temp[] = log(1 / nStates);

		// =================================== Calculating messages ==================================
		calculateMessages(nIt);

		// =================================== Calculating beliefs ===================================
#ifdef ENABLE_PPL
		concurrency::parallel_for_each(getGraphPairwise().m_vNodes.begin(), getGraphPairwise().m_vNodes.end(), [&, nStates](ptr_node_t &node) {
#else
		std::for_each(getGraphPairwise().m_vNodes.begin(), getGraphPairwise().m_vNodes.end(), [&,nStates](ptr_node_t &node) {
#endif
			for (size_t e_f : node->from) {
                msg_type *msg = getMessage(e_f);				// message of current incoming edge
                msg_type epsilon = FLT_EPSILON;
				for (byte s = 0; s < nStates; s++) { 		// states
                    node->Pot.at<float>(s,0) *= std::exp(msg[s]);
                    //node->Pot.at<float>(s, 0) = (epsilon + node->Pot.at<float>(s, 0)) * (epsilon + msg[s]);		// Soft multiplication
                } //s
			} // e_f
			
			// Normalization
            accum_type SUM_pot = 0;
			for (byte s = 0; s < nStates; s++)				// states
				SUM_pot += node->Pot.at<float>(s, 0);
			for (byte s = 0; s < nStates; s++) {			// states
				node->Pot.at<float>(s, 0) /= SUM_pot;
				DGM_ASSERT_MSG(!std::isnan(node->Pot.at<float>(s, 0)), "The lower precision boundary for the potential of the node %zu is reached.\n \
						SUM_pot = %f\n", node->id, SUM_pot);
			}
		});

		deleteMessages();
	}

	// dst: usually edge msg or edge msg_temp
    void CLogMessagePassing::calculateMessage(const Edge& edge_to, msg_type* temp, msg_type* dst, bool maxSum)
	{
		Node		* node = getGraphPairwise().m_vNodes[edge_to.node1].get();		// source node
		const byte	  nStates = getGraph().getNumStates();							// number of states

        // Compute temp = sum of all incoming msgs except e_t
		for (byte s = 0; s < nStates; s++) temp[s] = node->Pot.at<float>(s, 0);		// temp = node.Pot

		for (size_t e_f : node->from) {												// incoming edges
			Edge *edge_from = getGraphPairwise().m_vEdges[e_f].get();				// current incoming edge
			if (edge_from->node1 != edge_to.node2) {
                msg_type *msg = getMessage(e_f);										// message of current incoming edge
				for (byte s = 0; s < nStates; s++)
                    temp[s] += msg[s];												// temp = temp * msg
			}
        } // e_f

        // Compute new message: new_msg = (edge_to.Pot^2)^t x temp
        // but done in log domain
        const accum_type Z = MatMul(edge_to.Pot, temp, dst, maxSum);

        // Normalization and setting new values
        for (byte s = 0; s < nStates; s++) {
            dst[s] -= Z;
        }
	}

    void CLogMessagePassing::createMessages(std::optional<msg_type> val)
	{
		const size_t nEdges = getGraph().getNumEdges();
		const byte	nStates	= getGraph().getNumStates();
		
        m_msg = new msg_type[nEdges * nStates];
		DGM_ASSERT_MSG(m_msg, "Out of Memory");
        m_msg_temp = new msg_type[nEdges * nStates];
		DGM_ASSERT_MSG(m_msg_temp, "Out of Memory");

		if (val) {
			std::fill(m_msg, m_msg + nEdges * nStates, val.value());
			std::fill(m_msg_temp, m_msg_temp + nEdges * nStates, val.value());
		}
	}

    void CLogMessagePassing::deleteMessages(void)
	{
		if (m_msg) {
			delete[] m_msg;
			m_msg = NULL;
		}
		if (m_msg_temp) {
			delete[] m_msg_temp;
			m_msg_temp = NULL;
		}
	}

    void CLogMessagePassing::swapMessages(void)
	{
        msg_type *pTemp = m_msg;
		m_msg = m_msg_temp;
		m_msg_temp = pTemp;
	}

    CLogMessagePassing::msg_type*
    CLogMessagePassing::getMessage(size_t edge)
	{ 
		return m_msg ? m_msg + edge * getGraph().getNumStates() : NULL;
	}

    CLogMessagePassing::msg_type*
    CLogMessagePassing::getMessageTemp(size_t edge)
	{ 
		return m_msg_temp ? m_msg_temp + edge * getGraph().getNumStates() : NULL;
	}

	// dst = (M * M)^T x v
    CLogMessagePassing::accum_type
    CLogMessagePassing::MatMul(const Mat& M, const msg_type* v, msg_type* dst, bool maxSum)
	{
        accum_type res = -std::numeric_limits<accum_type>::infinity();
        DGM_ASSERT(dst);
		for (int x = 0; x < M.cols; x++) {
            accum_type sum = (maxSum) ? -std::numeric_limits<accum_type>::infinity() : 0;
			for (int y = 0; y < M.rows; y++) {
                msg_type m = M.at<float>(y, x);
                //XXX: recomputing logs a whole bunch here!
                msg_type prod = v[y] + msg_type{2}*std::log(m);
				if (maxSum) { if (prod > sum) sum = prod; }
                else sum += std::exp(prod); // could use log-sum-exp, but nStates is 2 for me
            } // y
            if (!maxSum) sum = std::log(sum);
            dst[x] = sum;
            res = std::max(res, sum);
        } // x
        return res;
	}
}
