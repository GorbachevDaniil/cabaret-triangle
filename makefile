TARGET			= cabaret_triangle
TEST_TARGET		= cabaret_triangle_test
MAIN_CLASS		= Main.cpp

SRC_DIR			= src
INCLUDE_DIR		= include
LIB_DIR			= lib
BIN_DIR			= bin
BUILD_DIR		= build
TEST_DIR		= test
OUTPUT_DIR		= $(BIN_DIR)/output
GTEST_DIR 		= $(LIB_DIR)/googletest

CXX       		= g++
CXXFLAGS  		= -g -std=c++11 -Wall -Wextra

LINKER   		= g++ -o
LFLAGS   		= -Wall -Wextra

INCLUDES 		= $(wildcard $(INCLUDE_DIR)/*.hpp)

SOURCES  		= $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS  		= $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
SOURCES_NO_MAIN	= $(filter-out $(SRC_DIR)/$(MAIN_CLASS), $(wildcard $(SRC_DIR)/*.cpp))
OBJECTS_NO_MAIN	= $(SOURCES_NO_MAIN:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
TEST_SOURCES	= $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJECTS	= $(TEST_SOURCES:$(TEST_DIR)/%.cpp=$(BUILD_DIR)/%.o)

$(BUILD_DIR)/gtest-all.o: $(GTEST_DIR)/src/gtest-all.cc
	$(CXX) $(CXXFLAGS) -I $(GTEST_DIR)/include -I $(GTEST_DIR) \
            -c $(GTEST_DIR)/src/gtest-all.cc \
            -o $(BUILD_DIR)/gtest-all.o

$(BUILD_DIR)/gtest_main.o: $(GTEST_DIR)/src/gtest_main.cc
	$(CXX) $(CXXFLAGS) -I $(GTEST_DIR)/include -I $(GTEST_DIR) \
            -c $(GTEST_DIR)/src/gtest_main.cc \
            -o $(BUILD_DIR)/gtest_main.o

$(BUILD_DIR)/gtest_main.a: $(BUILD_DIR)/gtest-all.o $(BUILD_DIR)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

$(BIN_DIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $@ $(LFLAGS) -I $(INCLUDE_DIR) $(OBJECTS)
	@echo "Target linked successfully!"

$(OBJECTS): $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(OUTPUT_DIR)
	@$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(BIN_DIR)/$(TEST_TARGET): $(TEST_OBJECTS) $(OBJECTS_NO_MAIN) $(BUILD_DIR)/gtest_main.a
	@$(LINKER) $@ $(LFLAGS) -I $(INCLUDE_DIR) -I $(GTEST_DIR)/include $(BUILD_DIR)/gtest_main.a $(TEST_OBJECTS) $(OBJECTS_NO_MAIN)
	@echo "Test target linked successfully!"

$(TEST_OBJECTS): $(BUILD_DIR)/%.o : $(TEST_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(GTEST_DIR)/include -c $< -o $@
	@echo "Compiled "$<" successfully!"	

target: $(BIN_DIR)/$(TARGET)

target_test: $(BIN_DIR)/$(TEST_TARGET)

clean:
	@rm -f $(OBJECTS) $(TEST_OBJECTS)
	@echo "Cleaned object files successfully!"
	@rm -f $(BIN_DIR)/$(TARGET) $(BIN_DIR)/$(TEST_TARGET)
	@echo "Cleaned bin files successfully!"
	@rm -rf $(OUTPUT_DIR)/
	@echo "Cleaned result files successfully!"