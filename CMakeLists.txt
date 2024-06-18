cmake_minimum_required(VERSION 3.10)

# Project name
project(handleWF)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add include directory
include_directories(include)

# Source files
set(SOURCE_FILES
    src/Atom.cpp
    # src/Field.cpp
    # src/Rvector.cpp
    # src/Timer.cpp
    # src/Wavefunction.cpp
    # src/main.cpp
)

# Test files
set(TEST_FILES
    test/test_atom.cpp
)

# Find Google Test
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Create executable for tests
add_executable(runTests ${SOURCE_FILES} ${TEST_FILES})

# Link Google Test libraries
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)

# Enable testing
enable_testing()

# Add test
add_test(NAME runTests COMMAND runTests)
