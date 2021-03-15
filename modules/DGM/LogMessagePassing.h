// Base abstract class for logarithmic message passing algorithms used for exact and approximate inference
// Written by James Folberth in 2021
#pragma once

#include "Infer.h"
#include "GraphPairwise.h"
#include "MessagePassing.h"

namespace DirectGraphicalModels
{
    class CLogMessagePassing : public CInfer
    {
    public:
        using msg_type = float;
        using accum_type = float;

        DllExport CLogMessagePassing(CGraphPairwise& graph) : CInfer(graph) {}
        DllExport virtual ~CLogMessagePassing(void) override = default;

        DllExport virtual void infer(unsigned int nIt = 1) override;

    protected:
        CGraphPairwise& getGraphPairwise(void) const {
            return dynamic_cast<CGraphPairwise&>(getGraph());
        }

        virtual void calculateMessages(unsigned int nIt) = 0;

        void calculateMessage(const Edge& edge, msg_type* temp, msg_type* dst, bool maxSum = false);

        void createMessages(std::optional<msg_type> val = std::nullopt);

        void deleteMessages(void);

        void swapMessages(void);

        msg_type* getMessage(size_t edge);

        msg_type* getMessageTemp(size_t edge);

        static accum_type MatMul(const Mat& M, const msg_type* v, msg_type* dst, bool maxSum = false);

    private:
        ///< Message: Mat(size: nStates x 1; type: CV_32FC1)
        msg_type* m_msg;
        ///< Temp Message: Mat(size: nStates x 1; type: CV_32FC1)
        msg_type* m_msg_temp;
    };
}
