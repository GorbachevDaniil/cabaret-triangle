#ifndef Parser_hpp
#define Parser_hpp

#include "Mesh.hpp"

class Parser {
   public:
    int LoadNodes(Mesh *mesh, std::string nodeFile);
    int LoadEdges(Mesh *mesh, std::string edgeFile);
    int LoadCells(Mesh *mesh, std::string cellFile);
};

#endif
