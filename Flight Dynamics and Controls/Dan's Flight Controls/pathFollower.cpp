// pathFollower.cpp

#include <math.h>  // Math operations, i.e. atan2
#include <iostream>
#include <cmath> // For abs function

#include "vector.hpp"

#define STRAIGHT 0
#define ORBIT 1
#define CLOCKWISE 1
#define CNTCLOCKWISE -1

#define PI M_PI
#define GRAVITY 9.81

#define K_ORBIT 0.2
#define K_PATH 0.03
#define CHI_INFINITY 50 * PI/180


#define sign(a) 2 * (a > 0) - 1

pathParams orbit(pathInputParams in);
pathParams lineFollow(pathInputParams in);

int main(int argc, char** argv){

}

pathParams orbit(pathInputParams in){

	pathParams out;

	double cOrbitN = in.cOrbit.xComponent;
	double cOrbitE = in.cOrbit.yComponent;
	double cOrbitD = in.cOrbit.zComponent;

	double rhoOrbit = in.orbitRadius;
	int lambda = in.orbitDirection;

	double pn = in.pn;
	double pe = in.pe;
	double height = in.height;

	double phiOrbit = atan2(pe - cOrbitE, pn - cOrbitN);
	double chi = in.chiRadians;

	while (std::abs(phiOrbit - chi) < PI){
		phiOrbit = phiOrbit - 2 * PI * sign(phiOrbit - chi);
	}

	double distFromCenter = std::sqrt( pow(pn - cOrbitN, 2) + pow(pe - cOrbitE, 2) ); 

	out.chiCommandedRadians = phiOrbit + lambda * (PI / 2 + atan(K_ORBIT * (distFromCenter - rhoOrbit) / rhoOrbit) );

    out.heightCommanded = -1 * cOrbitD;
    out.phiFFradians = lambda * atan(pow(in.VaCurrent, 2) / GRAVITY * rhoOrbit);
    out.VaC = in.VaD;

    return out;
}

pathParams lineFollow(pathInputParams in){

	pathParams out;

	double qPathN = in.qPath.xComponent;
	double qPathE = in.qPath.yComponent;
	double qPathD = in.qPath.zComponent;

	double rPathN = in.rPath.xComponent;
	double rPathE = in.rPath.yComponent;
	double rPathD = in.rPath.zComponent;

	double pn = in.pn;
	double pe = in.pe;
	double height = in.height;

	double chiQ = atan2(qPathE, qPathN);  // Course angle as measured from north.
	double chi = in.chiRadians;

	while (std::abs(chi - chiQ) > 2*PI){
		chiQ = chiQ + 2 * PI * sign(chi - chiQ);
	}

	Vector qVector = Vector(qPathN, qPathE, qPathD);
	Vector normalVector = Vector(0, 0, 1.0);
	Vector qCrossK = crossProduct(qVector, normalVector);
	Vector n = normalize(qCrossK);

	squareMatrix rotationMatrix = squareMatrix(cos(chiQ), sin(chiQ), 0, 
								-1*sin(chiQ), cos(chiQ), 0, 
								 0, 0, 1);
	Vector ePI = Vector(pn - rPathN, pe - rPathE, height - rPathD);

	Vector eP = rotationMatrix * ePI;

	double e_py = eP[1];

	Vector sI = ePI - dotProduct(ePI, n) * n;

	out.chiCommandedRadians = chiQ - CHI_INFINITY * 2 / PI * atan(K_PATH * e_py);
	out.heightCommanded = -rPathD + std::sqrt(pow(sI[0], 2) + pow(sI[1], 2) ) * qVector[2] / sqrt(pow(qVector[1], 2) + pow(qVector[0], 2) );
	out.phiFFradians = in.phiRadians; // atan(pow(in.VaCurrent, 2) / GRAVITY * )
	out.VaC = in.VaD;

    return out;
}