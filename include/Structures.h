#ifndef STRUCTURES_HPP_
#define STRUCTURES_HPP_

#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <vector>

class Coords{
public:
	double x;
	double y;
	double z;

	inline void Clear() { x = 0.;  y = 0.;  z = 0.; }
	inline void Set(Coords &c){ x = c.x; y = c.y; z = c.z; }
	inline void SetMult(Coords &c, double d){ x = d*c.x; y = d*c.y; z = d*c.z;}
	inline void Set(double *c){ x = c[0]; y = c[1]; z = c[2];}
	inline void Set(double c0, double c1, double c2 ){x = c0; y = c1; z = c2;}
	inline double Scalar(Coords &c){ return (x*c.x + y*c.y + z*c.z); }
	static double ScalarProduct(Coords &v1, Coords &v2) {
		return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
	}

	static double ScalarProduct(double *d, Coords &v) {
		return d[0]*v.x+d[1]*v.y+d[2]*v.z;
	}

	static double ScalarProduct(double *d1, double *d2) {
		return d1[0]*d2[0]+d1[1]*d2[1]+d1[2]*d2[2];
	}
	inline void VectorProduct(Coords &a, Coords &b){
		x = a.y*b.z - a.z*b.y;
		y = a.z*b.x - a.x*b.z;
		z = a.x*b.y - a.y*b.x;}
	inline void Mult(double d){ x *= d; y *= d; z *= d; }
	inline void Plus(Coords &c){ x += c.x; y += c.y; z += c.z; }
	inline void PlusMult(Coords &c, double d){ x += d*c.x; y += d*c.y; z += d*c.z; }
	inline void Middle(Coords &a, Coords &b) {x = 0.5*(a.x + b.x); y = 0.5*(a.y + b.y); z = 0.5*(a.z + b.z);}
	inline double Length(){ return sqrt(x*x + y*y + z*z); }
	inline double Length2(){ return (x*x + y*y + z*z); }
	inline void SetMinus(Coords &c1, Coords &c2) { x=c1.x-c2.x; y=c1.y-c2.y; z=c1.z-c2.z; }
	static double Length(double x, double y, double z) { return sqrt(x*x+y*y+z*z); }
	static double GetDistance(Coords &v1, Coords &v2) { return Coords::Length(v2.x-v1.x,v2.y-v1.y,v2.z-v1.z); }

};

class Var{
public:
	double rho;
	Coords u;
	double p; // давление в слабосжимаемом приближении p = sound2*(rho-rho0)
	double T;
	double *ksi;

	inline void Clear(int nComp){ rho = 0.; u.Clear(); p = 0.; T = 0.; for(int i = 0; i<nComp; i++) ksi[i] = 0.;};
	inline void Set(int nComp,Var &c){ rho = c.rho; u.Set(c.u); p = c.p; T = c.T; for(int i = 0; i<nComp; i++) ksi[i] = c.ksi[i];};
	inline void CalcP(double sound2, double P0, double gamma, double R){ p = sound2*rho*(1. - pow(P0/(rho*R*T), 1./gamma)); }
	inline void Mult(int nComp,double d) {rho *= d; u.Mult(d); p *= d; T *= d; for(int i = 0; i<nComp; i++) ksi[i] *= d;}
	inline void PlusMult(int nComp, Var &c, double d) {rho += d*c.rho; u.PlusMult(c.u,d); p += d*c.p; T += d*c.T; for(int i = 0; i<nComp; i++) ksi[i] += d*c.ksi[i];}
	inline void CalcRho(double sound2, double P0, double gamma, double Cv) {rho = P0/(gamma-1)/Cv/T + gamma*p/sound2;}
	inline void CalcRhoNonLin(double sound2, double P0, double gamma, double Cv)
	{
		double rhos = P0/(gamma-1.)/Cv/T;
		double drho = p/sound2;
		double b = drho/rhos;
		double x = b;
		double eps = 1.e-5;
		if(abs(b) > 1.e-10){
			while(true){
				double x0 = x;
				double pw1 = pow((1.-x0), gamma);
				double pw2 = pow((1.-x0), (gamma-1.));
				x = x0 +  (b*pw1-x0)/(b*gamma*pw2 + 1);
				if(abs(x-x0)<eps*abs(x)) break;
			}

			rho = drho/x;
		}else{
			rho = rhos + gamma*drho;
		}
	}

	inline void CalcRhoQ(double sound2, double P0, double gamma, double Cv, double I, double G, double Q)
	{
		double rhos = P0/(gamma-1.)/Cv/T;
		double drho = I/(G*sound2);
		double b = drho/rhos;
		double a = Q/rhos/b/I;
		double x = b;
		double eps = 1.e-5;
		if(abs(b) > 1.e-10){
			while(true){
				double x0 = x;
				double pw1 = pow((1.-x0+a*x0*x0), gamma);
				double pw2 = pow((1.-x0+a*x0*x0), (gamma-1.));
				x = x0 + (b*pw1-x0)/( 1. - b*gamma*pw2*(2.*a*x0-1.));
				if(abs(x-x0)<eps*abs(x)) break;
			}

			rho = drho/x;
		}else{
			rho = rhos + gamma*drho;
		}
	}
};

typedef struct{
	int n1, n2, n3;
	double *x, *y, *z;
	double lx, ly, lz;
	double Volume;
}Mesh;


typedef struct{
	int type;
	double pars[20];
}SurfType;

#endif /* STRUCTURES_HPP_ */
