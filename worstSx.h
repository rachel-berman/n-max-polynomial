/*
 * worstSx.h
 *
 *  Created on: Jul 4, 2017
 *      Author: rberman
 */

#ifndef WORSTSX_H_
#define WORSTSX_H_

#include <fstream>
#include <iostream>
#include <new>
#include <cmath>
#include <iomanip>
#include "auxFunctions.h"

#define ADVANCED_VERSION 1 //according to theorem 2.4

using namespace std;

extern int globalDeg_sx;

struct WorstSx
{
	int* sxRegularHist; /* pointer to array of length regularHistLength */
	int* sxTotalHist; /* pointer to array of length totalHistLength */
	unsigned long int numValidUx;
	int t;
	int dt;
};

WorstSx* WorstSxCreate(DegLengthINFO* degLenInfo);
void WorstSxDestroy(WorstSx* worstSx);

STATUS compareToWorstSx(DegLengthINFO* degLenInfo, WorstSx* worstSx,
		int* currentSxTotal);

/* check if s(x) in total histogram satisfies necessary conditions:
properties 4+5 and deg(s) >= minDeg.
(deg(s) <= 2n + roperties 1+3 were already checked in the backtracking algorithm)*/
bool isSxPotentialyWorst(DegLengthINFO* degLenInfo, int* Sx, int minDeg);

// MAIN FUNCTIONS
void getWorstSx(int n, int q, int minDeg, int maxDeg);
void getWorstSxAdvanced(int n, int q, bool printPossibleWorstPoly);


void printWorstSx(DegLengthINFO* degLenInfo, WorstSx* worstSx);

//
//template<typename T> void printElement(T t, const int& width, const char& separator)
//{
//    cout << left << setw(width) << setfill(separator) << t;
//}

// Sx is given in a regular histogram form
int numValidUxDividingSx(DegLengthINFO* degLenInfo, int* Sx);

/* GIVEN: total histogram of s(x).
 * calculate the regular histogram and write it into sRegularHist,
 * RETURN: number of valid u(x) OR -1 in case of error
 * OUTPUT PARAMETER: sRegularHist.
 */
int getRegularHistAndNumValidUx(DegLengthINFO* degLenInfo, int* sTotalHist,
		 int* sRegularHist );

#endif /* WORSTSX_H_ */

