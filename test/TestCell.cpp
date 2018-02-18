#include "gtest/gtest.h"

#include "Cell.hpp"

#include <vector>

TEST(getNextNodeIDTest, ReturnTheNextValueInACircularManner) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    EXPECT_EQ(11, cell.getNextNodeID(0));
    EXPECT_EQ(12, cell.getNextNodeID(1));
    EXPECT_EQ(10, cell.getNextNodeID(2));
}

TEST(getNextNodeIDTest, DieIfIncomingNodeIDPosMoreThanSize) {
    Cell cell;

    ASSERT_DEATH(cell.getNextNodeID(1), "");
}

TEST(getPrevNodeIDTest, ReturnTheNextValueInACircularManner) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    EXPECT_EQ(12, cell.getPrevNodeID(0));
    EXPECT_EQ(10, cell.getPrevNodeID(1));
    EXPECT_EQ(11, cell.getPrevNodeID(2));
}

TEST(getPrevNodeIDTest, DieIfIncomingNodeIDPosMoreThanSize) {
    Cell cell;

    ASSERT_DEATH(cell.getNextNodeID(1), "");
}

TEST(getEdgeOrderedNodeIDs, UnorderedNodeIDs) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    std::vector<long> unorderedNodeIDs;
    unorderedNodeIDs.push_back(10);
    unorderedNodeIDs.push_back(12);

    std::vector<long> orderedNodeIDs = cell.getEdgeOrderedNodeIDs(unorderedNodeIDs);
    EXPECT_EQ(12, orderedNodeIDs[0]);
    EXPECT_EQ(10, orderedNodeIDs[1]);
}

TEST(getEdgeOrderedNodeIDs, OrderedNodeIDs) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    std::vector<long> unorderedNodeIDs;
    unorderedNodeIDs.push_back(12);
    unorderedNodeIDs.push_back(10);

    std::vector<long> orderedNodeIDs = cell.getEdgeOrderedNodeIDs(unorderedNodeIDs);
    EXPECT_EQ(12, orderedNodeIDs[0]);
    EXPECT_EQ(10, orderedNodeIDs[1]);
}

TEST(getEdgeOrderedNodeIDs, SizeMustBeTwoForIncomingVector) {
    Cell cell;

    std::vector<long> unorderedNodeIDs;
    unorderedNodeIDs.push_back(10);
    unorderedNodeIDs.push_back(11);
    unorderedNodeIDs.push_back(12);

    ASSERT_DEATH(cell.getEdgeOrderedNodeIDs(unorderedNodeIDs), "");

    unorderedNodeIDs.clear();
    unorderedNodeIDs.push_back(10);

    ASSERT_DEATH(cell.getEdgeOrderedNodeIDs(unorderedNodeIDs), "");
}

TEST(getEdgeOrderedNodeIDs, SizeMustBeTwoForOutcomingVector) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    std::vector<long> unorderedNodeIDs;
    unorderedNodeIDs.push_back(11);
    unorderedNodeIDs.push_back(11);

    ASSERT_DEATH(cell.getEdgeOrderedNodeIDs(unorderedNodeIDs), "");
}