cmake_minimum_required(VERSION 3.11)

include(FetchContent)

project(cabaret_triangle)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)
file(GLOB_RECURSE SOURCES "src/*.cpp")

find_package(Armadillo REQUIRED)
include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})

# cabaret target begin
add_executable(cabaret Main.cpp ${SOURCES})
target_compile_options(cabaret PRIVATE -Wall -Wextra -pedantic -Werror)
target_link_libraries(cabaret ${ARMADILLO_LIBRARIES})
# engine target end

# unit_tests target begin
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        release-1.10.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_GetProperties(googletest)
if(NOT googletest)
    FetchContent_Populate(googletest)
    add_subdirectory(
            ${googletest_SOURCE_DIR}
            ${googletest_BINARY_DIR}
            EXCLUDE_FROM_ALL
    )
endif()

enable_testing()
file(GLOB TEST_SOURCES test/*.cpp)
add_executable(unit_tests ${SOURCES} ${TEST_SOURCES})
target_link_libraries(unit_tests gtest_main ${ARMADILLO_LIBRARIES})
add_test(test_all unit_tests)

file(COPY test/resources DESTINATION ${CMAKE_BINARY_DIR}/test)
# unit_tests target end
