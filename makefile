TARGET		= cabaret_triangle

SRC_DIR		= src
INCLUDE_DIR	= include
LIB_DIR		= lib
BIN_DIR		= bin
BUILD_DIR	= build
TEST_DIR	= test

CXX       	= g++
CXXFLAGS  	= -std=c++11 \
		 		-Wall \
		 		-I $(INCLUDE_DIR) \
		 		-I $(LIB_DIR)

LINKER   	= g++ -o
LFLAGS   	= -Wall \
		 		-I $(INCLUDE_DIR) \
		 		-I $(LIB_DIR)

SOURCES  	= $(wildcard $(SRC_DIR)/*.cpp)
INCLUDES 	= $(wildcard $(INCLUDE_DIR)/*.hpp)
OBJECTS  	= $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

$(BIN_DIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $@ $(LFLAGS) $(OBJECTS)
	@echo "Linked successfully!"

$(OBJECTS): $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

test:
	@$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJECTS)
	@echo "Cleaned object files successfully!"

remove: clean
	@rm -f $(BIN_DIR)/$(TARGET)
	@echo "Cleaned bin files successfully!"