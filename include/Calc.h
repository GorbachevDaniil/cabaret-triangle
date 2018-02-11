#ifndef CALC_HPP_
#define CALC_HPP_

#include <stdio.h>
#include <vector>
#include <map>
#include "Structures.h"
#include "Node.h"
#include "Cell.h"
#include "Side.h"


using namespace std;

/*
 *  обращение матрицы: to = from^(-1)
 */
inline int Inverse(double from[3][3], double to[3][3]) {
	double det =
			from[0][0] * (from[1][1]*from[2][2]-from[1][2]*from[2][1]) -
			from[0][1] * (from[1][0]*from[2][2]-from[1][2]*from[2][0]) +
			from[0][2] * (from[1][0]*from[2][1]-from[1][1]*from[2][0]);
	if(det==0) {
		printf("Inverse: Null matrix.\n");
		return 0;
		//Exit();
	}
	double det1 = 1. / det;
  to[0][0] = det1 * (from[1][1] * from[2][2] - from[1][2] * from[2][1]);
  to[1][0] = det1 * (from[1][2] * from[2][0] - from[1][0] * from[2][2]);
  to[2][0] = det1 * (from[1][0] * from[2][1] - from[1][1] * from[2][0]);

  to[0][1] = det1 * (from[0][2] * from[2][1] - from[0][1] * from[2][2]);
  to[1][1] = det1 * (from[0][0] * from[2][2] - from[0][2] * from[2][0]);
  to[2][1] = det1 * (from[0][1] * from[2][0] - from[0][0] * from[2][1]);

  to[0][2] = det1 * (from[0][1] * from[1][2] - from[0][2] * from[1][1]);
  to[1][2] = det1 * (from[0][2] * from[1][0] - from[0][0] * from[1][2]);
  to[2][2] = det1 * (from[0][0] * from[1][1] - from[0][1] * from[1][0]);
  return 1;
}

//----------------------------------------------------------------------------------------------------------------------------

// умножение матрицы на вектор: result = m * r
inline void MatrixVectorProduct(double m[3][3], double r[3], double result[3]) {
  result[0] = m[0][0]*r[0] + m[0][1]*r[1] + m[0][2]*r[2];
  result[1] = m[1][0]*r[0] + m[1][1]*r[1] + m[1][2]*r[2];
  result[2] = m[2][0]*r[0] + m[2][1]*r[1] + m[2][2]*r[2];
}

//----------------------------------------------------------------------------------------------------------------------------

// давление насыщенного пара:
inline double Ps(double T){

	double T0 = 273.15;
	double P0 = 611.3;
	double T1 = 38.35;
	double cof = 17.05;

	return ( P0*exp(cof*(T-T0)/(T-T1)) ); // формула из аппрксимации экспериментальных данных

}

inline double DPs(double T){

	double T0 = 273.15;
	//double P0 = 611.3;
	double T1 = 38.35;
	double cof = 17.05;

	return ( cof*(T0 - T1)/pow((T-T1),2)*Ps(T) );

}
////////////////////////////////////////////////////////////////////////////////////////////


class Calc{
public:

	char *saveFile;

	int iter;
	int nstop, nout;
	int nsave;
	double pan;
	double tstop, tout;
	double dt, ttime, CFL;
	double sound, sound2;
	map<int, double> rho0;
	map<int, double> kapa;
	map<int, double> Cv;
	map<int, double> Qsource;
	int nComp;
	map<int, double> gammai;
	map<int, double> Mi;
	map<int, double> Ri;
	map<int, double> Cvi;
	map<int, double> Cpi;
	map<int, double> Tki;
	map<int, double> Pki;
	map<int, double> Vki;
	map<int, double> MUki;
	map<int, double> Dvi;
	double **D0;
	double Lw;
	double R0, P0, rho_0;
	double g[3];
	double tvisc;
	double OMEGA;
	int nTask;
	double lxmin, lymin, lzmin, lmin;
	bool needAv;
	double av_start, av_stop;
	double avTime;

	int nMat;
	bool needCalcT;
	bool needCalcV;
	bool needCalcD;
	bool needTransportCoeff;
	bool havePBC;

	double timeToOutput;
	int nOutput;
	bool needOutput;

	int scrout;
	bool needCond;

	double filter;

	double NaN;
	double pi;

	map<int, Node *>nRef; // все узлы
	int nNodes;
	map<int, Node *>nThermRef; // все тепловые узлы
	int nThermNodes;

	map<int, Cell *>cRef; // все ячейки
	int nCells;
	map<int, Side *>sRef; // все грани
	int nSides;
	map<int, Cell *>cHydroRef; // все гидродинамические ячейки
	int nHydroCells;
	map<int, Side *>sHydroRef; // все гидродинамические грани
	int nHydroSides;
	map<int, Cell *>cThermRef; // все тепловые ячейки
	int nThermCells;
	map<int, Side *>sThermRef; // все тепловые грани
	int nThermSides;
	map<int, Cell *>cMatRef;   // все ячейки материала
	int nMatCells;
	map<int, int> RefSides;

	int nSurf;
	map<int, vector<Side *> > Surf; // список внешних границ
	SurfType *SurfTypesV, *SurfTypesT, *SurfTypesD;

	vector<Side *> bc_cond;


	double IntTau;
	FILE *fInt, *fWCond;


	bool WallCond;
	double dM_cond, dM_cond_sum;

	Calc();
	~Calc();

	int CleanOut();
	int Parameters();
	int BndConditions();
	int InitMesh(Mesh *GRID);
	int BuildConnectivity(Mesh *GRID);
	int FixNodeAdresses(Mesh *GRID);
	int FixCellAdresses(Mesh *GRID);
	int FixSideAdresses(Mesh *GRID);
	int WebCutting();
	int SetMaterials();
	int BuildCoef();
	int SetSurfaces();
	int SetPeriodicBoundaries();
	int Init(Mesh *);
	int InitSteps();
	int Steps(Mesh *);
	int CalcDT();
	int Phase1();
	int Phase2();
	int Phase3();
	int BuildFlux();
	int VPhase();
	int TPhase();
	int DPhase();
	int Condensation();
	int Swapvar();
	int Output(Mesh *GRID);
	int TecplotOutput();
	int WriteSave();
	int ReadSave(char *filename);
	int Averaging();
	int AvOutput();
};

#endif /* CALC_HPP_ */
