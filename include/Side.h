#ifndef SIDE_HPP_
#define SIDE_HPP_

#include "Structures.h"

class Cell;
class Node;
class Calc;

class Side {
	public:
		int gid;
		Coords center;
		double S;
		Coords n,m,l; //базовые вектора
		Cell *adjCells[2];
		Side *op[2];
		Side *pbc;
		Side *jsides[2][4];
		double C[2][4];
		double C0[2];

		Var uf, ufn;
		Coords F, Fn;
		double Q, Qn;
		double *J, *Jn;
		double Cv, Cp;


		int BndNum;
		int type; // 0 - внутренняя; 1 - граничная;
		int surfNo;

		Side(int nComp);
		~Side();

		int Phase2(Calc *);
		int VPhase(Calc *);
		int TPhase(Calc *);
		int DPhase(Calc *);
		double GetInv(int, Calc *, Side *, double);

		int Bound(Calc *);
		int InOut(Calc *);
		int Slip(Calc *);
		int NoSlip(Calc *);
		int FixedQ(Calc *);

		int BoundV(Calc *);
		int InOutV(Calc *);
		int SlipV(Calc *);
		int NoSlipV(Calc *);
		int FixedQV(Calc *);

		int BoundT(Calc *);
		int I_type_T(Calc *);
		int II_type_T(Calc *);
		int III_type_T(Calc *);

		int BoundD(Calc *);
		int I_type_D(Calc *);
		int II_type_D(Calc *);
		int III_type_D(Calc *);
		int WallCondD(Calc *);

};

#endif /* SIDE_HPP_ */
