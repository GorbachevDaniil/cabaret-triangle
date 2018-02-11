#ifndef NODE_HPP_
#define NODE_HPP_

#include "Calc.h"
#include "Structures.h"

class Cell;
class Side;
class Calc;

class Node{

	public:

		int gid;
		int id;  // номер для вывода в текплот файл
		Coords center;
		Cell *adjCells[8];

		Node();
		~Node();
};



#endif /* NODE_HPP_ */
