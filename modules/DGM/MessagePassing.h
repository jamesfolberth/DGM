// Base abstract class for message passing algorithms used for exact and approximate inference
// Written by Sergey G. Kosov in 2015 for Project X
#pragma once

#include "Infer.h"
#include "GraphPairwise.h"

namespace DirectGraphicalModels
{
	struct Edge;

	// ==================== Message Passing Base Abstract Class ==================
	/**
	* @ingroup moduleDecode
	* @brief Abstract base class for message passing inference algorithmes
	* @author Sergey G. Kosov, sergey.kosov@project-10.de
	*/
	class CMessagePassing : public CInfer
	{
	public:
		/**
		* @brief Constructor
		* @param graph The graph
		*/
		DllExport CMessagePassing(CGraphPairwise &graph) : CInfer(graph) {}
		DllExport virtual ~CMessagePassing(void) = default;
		
		DllExport virtual void	  infer(unsigned int nIt = 1);


	protected:
		/**
		* @brief Returns the graph
		* @return The graph
		*/
		CGraphPairwise& getGraphPairwise(void) const { return dynamic_cast<CGraphPairwise&>(getGraph()); }
		/**
		* @brief Calculates messages, associated with the edges of corresponding graphical model
		* @details > This function may modify Edge::msg and Edge::msg_temp containers of graph edges
		* @param nIt Number of iterations
		*/
		virtual void calculateMessages(unsigned int nIt) = 0;
		/**
		* @brief Calculates one message for the specified edge \b edge
		* @details > PPL-safe function.
		* @param[in] edge Graph edge
		* @param[in] temp Auxilary array of \b nStates values. Introduced for higher perfomance reasons.
		* @param[out] dst Destination array for calculated message. Usually \b edge->msg or \b edge->msg_temp.
		* @param[in] maxSum Flag indicating weather the message must be calculated according to the \a sum-product (false) or \a max-product (true) algorithm.
		*/
		void	calculateMessage(const Edge& edge, double* temp, double* dst, bool maxSum = false);
		/**
		* @brief Allocates memory for Edge::msg and Edge::msg_temp containers for all edges in the graph
		* @param val Default value to fill in the Edge::msg and Edge::msg_temp containers 
		*/
		void	createMessages(std::optional<double> val = std::nullopt);
		/**
		* @brief Deletes memory for Edge::msg and Edge::msg_temp containers for all edges in the graph
		*/
		void	deleteMessages(void);
		/**
		* @brief Swaps Edge::msg and Edge::msg_temp for all edges in the graph 
		*/
		void	swapMessages(void);
		/**
		* @brief Returns the pointer to the edge messages
		* @param edge The %Edge index
		* @return The pointer to the edge messages
		*/
		double*	getMessage(size_t edge);
		/**
		* @brief Returns the pointer to the edge temp messages
		* @param edge The %Edge index
		* @return The pointer to the edge temp messages
		*/
		double*	getMessageTemp(size_t edge);
		/**
		* @brief Specific matrix multiplication
		* @details This function calculates the result of multiplying square of matrix \b M by vector \b v as following:
		* \f$\vec{dst} = (M\cdot M)^\top\times\vec{v}\f$
		* @param[in] M Matrix: Mat(size: \a M.height x \a M.width; type: CV_32FC1)
		* @param[in] v Vector of length \a M.height
		* @param[out] dst Resulting vector of length \a M.width.
		* @param[in] maxSum Flag indicating weather the \a max-sum multiplication should be performed
		* @return The sum of all elemts in vector \b dst
		*/
		static double MatMul(const Mat& M, const double* v, double* dst, bool maxSum = false);


	private:
		double	* m_msg;			///< Message: Mat(size: nStates x 1; type: CV_32FC1)
		double	* m_msg_temp;		///< Temp Message: Mat(size: nStates x 1; type: CV_32FC1)
	};
}
