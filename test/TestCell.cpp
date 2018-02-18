#include "gtest/gtest.h"

#include "Cell.hpp"

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