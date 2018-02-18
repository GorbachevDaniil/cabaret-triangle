#include "Cell.hpp"

#include <vector>
#include <algorithm>
#include <cassert>

long Cell::getNextNodeID(unsigned long nodeIDPos) {
    assert(nodeIDPos < nodeIDs.size());

    if (nodeIDPos == nodeIDs.size() - 1) {
        return nodeIDs.front();
    }
    return nodeIDs[nodeIDPos + 1];
}

long Cell::getPrevNodeID(unsigned long nodeIDPos) {
    assert(nodeIDPos < nodeIDs.size());

    if (nodeIDPos == 0) {
        return nodeIDs.back();
    }
    return nodeIDs[nodeIDPos - 1];
}

std::vector<long> Cell::getEdgeOrderedNodeIDs(std::vector<long> edgeUnorderedNodeIDs) {
    assert(edgeUnorderedNodeIDs.size() == 2);

    std::vector<long> edgeOrderedNodeIDs;
    for (unsigned long i = 0; i < nodeIDs.size(); i++) {
        if (nodeIDs[i] == edgeUnorderedNodeIDs[0]) {
            if (getNextNodeID(i) == edgeUnorderedNodeIDs[1]) {
                edgeOrderedNodeIDs = edgeUnorderedNodeIDs;
                break;
            }
            if (getPrevNodeID(i) == edgeUnorderedNodeIDs[1]) {
                std::reverse(edgeUnorderedNodeIDs.begin(), edgeUnorderedNodeIDs.end());
                edgeOrderedNodeIDs = edgeUnorderedNodeIDs;
                break;
            }
        }
    }
    
    assert(edgeOrderedNodeIDs.size() == 2);
    return edgeOrderedNodeIDs;
}