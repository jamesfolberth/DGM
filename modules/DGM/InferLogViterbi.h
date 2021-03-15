// log-Viterbi inference class interface
// Writen by James Folberth in 2021
#pragma once

#include "InferLogLBP.h"

namespace DirectGraphicalModels
{
	// ==================== Viterbi Infer Class ==================
	/**
	* @ingroup moduleDecode
	* @brief Max product Viterbi inference class
	* @author Sergey G. Kosov, sergey.kosov@project-10.de
	*/
    class CInferLogViterbi : public CInferLogLBP
	{
	public:
		/**
		* @brief Constructor
		* @param graph The graph
		*/			
        DllExport CInferLogViterbi(CGraphPairwise &graph) : CInferLogLBP(graph) { setMaxSum(true); }
        DllExport virtual ~CInferLogViterbi(void) = default;
	};

}
