#include <gtest/gtest.h>

#include "Atom.hpp"

TEST(AtomTest, ConstructorFromIndividualCoordinates) {
    Atom atom(1, 1.0, 2.0, 3.0);

    EXPECT_EQ(atom.get_atnum(), 1);
    EXPECT_DOUBLE_EQ(
        atom.getMass(),
        1.008);  // Assuming this is the correct mass for atomic number 1
    EXPECT_EQ(atom.getSymbol(), "H");
    EXPECT_DOUBLE_EQ(atom.get_x(), 1.0);
    EXPECT_DOUBLE_EQ(atom.get_y(), 2.0);
    EXPECT_DOUBLE_EQ(atom.get_z(), 3.0);
}

TEST(AtomTest, ConstructorFromArray) {
    double coordinates[] = {1.0, 2.0, 3.0};
    Atom atom(1, coordinates);

    EXPECT_EQ(atom.get_atnum(), 1);
    EXPECT_DOUBLE_EQ(atom.getMass(), 1.008);
    EXPECT_EQ(atom.getSymbol(), "H");
    EXPECT_DOUBLE_EQ(atom.get_x(), 1.0);
    EXPECT_DOUBLE_EQ(atom.get_y(), 2.0);
    EXPECT_DOUBLE_EQ(atom.get_z(), 3.0);
}

TEST(AtomTest, ConstructorFromVector) {
    std::vector<double> coordinates = {1.0, 2.0, 3.0};
    Atom atom(1, coordinates);

    EXPECT_EQ(atom.get_atnum(), 1);
    EXPECT_DOUBLE_EQ(atom.getMass(), 1.008);
    EXPECT_EQ(atom.getSymbol(), "H");
    EXPECT_DOUBLE_EQ(atom.get_x(), 1.0);
    EXPECT_DOUBLE_EQ(atom.get_y(), 2.0);
    EXPECT_DOUBLE_EQ(atom.get_z(), 3.0);
}

TEST(AtomTest, ConstructorFromSymbolAndIndividualCoordinates) {
    Atom atom("H", 1.0, 2.0, 3.0);

    EXPECT_EQ(atom.get_atnum(), 1);
    EXPECT_DOUBLE_EQ(atom.getMass(), 1.008);
    EXPECT_EQ(atom.getSymbol(), "H");
    EXPECT_DOUBLE_EQ(atom.get_x(), 1.0);
    EXPECT_DOUBLE_EQ(atom.get_y(), 2.0);
    EXPECT_DOUBLE_EQ(atom.get_z(), 3.0);
}

TEST(AtomTest, ConstructorFromSymbolAndArray) {
    double coordinates[] = {1.0, 2.0, 3.0};
    Atom atom("H", coordinates);

    EXPECT_EQ(atom.get_atnum(), 1);
    EXPECT_DOUBLE_EQ(atom.getMass(), 1.008);
    EXPECT_EQ(atom.getSymbol(), "H");
    EXPECT_DOUBLE_EQ(atom.get_x(), 1.0);
    EXPECT_DOUBLE_EQ(atom.get_y(), 2.0);
    EXPECT_DOUBLE_EQ(atom.get_z(), 3.0);
}

TEST(AtomTest, ConstructorFromSymbolAndVector) {
    std::vector<double> coordinates = {1.0, 2.0, 3.0};
    Atom atom("H", coordinates);

    EXPECT_EQ(atom.get_atnum(), 1);
    EXPECT_DOUBLE_EQ(atom.getMass(), 1.008);
    EXPECT_EQ(atom.getSymbol(), "H");
    EXPECT_DOUBLE_EQ(atom.get_x(), 1.0);
    EXPECT_DOUBLE_EQ(atom.get_y(), 2.0);
    EXPECT_DOUBLE_EQ(atom.get_z(), 3.0);
}

TEST(AtomTest, OperatorOutput) {
    Atom atom(1, 1.0, 2.0, 3.0);
    std::stringstream ss;
    ss << atom;

    EXPECT_EQ(ss.str(),
              "  H     1    1.0080");  // Adjust format as per your ostream
                                       // operator<< definition
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
