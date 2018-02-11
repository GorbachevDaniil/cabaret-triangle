#ifndef CELL_HPP_
#define CELL_HPP_

#include "Calc.h"
#include "Structures.h"

class Side;
class Calc;
class Node;

class Cell {
	public:
		int gid;
		Coords center;
		Coords h;
		double Vol;
		Side *adjSides[6];
		Node *adjNodes[8];
		int sense[6];
		int material;
		bool CutCell;

		Var uc;
		Var ucn;
		Var uc_av;

		double G;
		double *Flux;
		double *D;    // коэффициент диффузии смеси
		double mu;    // вязкость смеси
		double kapa;  // теплопроводность смеси
		double Cp, Cv;// теплоемкость смеси
		double CpN, CvN;

		double TransportCoff;

		Cell(int nComp);
		~Cell();
		int Phase1H(Calc *);
		int Phase1M(Calc *);
		int Phase3H(Calc *);
		int Phase3M(Calc *);
		int BuildFluxH(Calc *);
		int BuildFluxM(Calc *);
		int BuildFlux0H(Calc *);
		int BuildFlux0M(Calc *);
		int CalcTransportCoeff0(Calc *);
		int CalcTransportCoeff(Calc *);
		int VPhase(Calc *);
		int VFluxes(Calc *);
		int TPhase(Calc *);
		double Inv(int, Side *, Calc *);
		double GetInv(int, int, Calc *, Side *);
		double GetChar(int, Calc *, Side *);
		int Averaging(Calc *);
		int Condensation(Calc *);
};

#endif /* CELL_HPP_ */
