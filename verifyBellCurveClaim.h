/*
 * verifyBellCurveClaim.h
 *
 *  Created on: Oct 1, 2017
 *      Author: rberman
 */

#ifndef VERIFYBELLCURVECLAIM_H_
#define VERIFYBELLCURVECLAIM_H_

#include <fstream>
#include <iostream>
#include <sstream>      // std::istringstream
#include <string>       // std::string
#include <new>

#include <iomanip>
#include "auxFunctions.h"
using namespace std;

#define MAX_HIST_LENGTH 100

 /* PARSING RESULT.TXT FILE - DECLARATIONS */
struct SxInfo{
	int n;
	int totalLength;
	int* sxHist;
	int Lc;
	int degSx;
};
SxInfo* SxInfoCreate();
void SxInfoDestroy(SxInfo* sxInfo);
//SxInfo* parsingSxLine(istringstream &iss);

SxInfo* parsingSxLine(string line);

void printSxInfo(SxInfo* sxInfo);
/* PARSING - end of declarations*/


/* Deal only with REGULAR HISTOGRAMS in COMPLEMENTARY polynomials model */
void printSizeBkCurve(int n, int q, int *SxHist, int SxHistLength, int Lc);

void printSiCurve(SxInfo* sxInfo);

void printCurve(int *curveValues, int len);
// DegLengthINFO* degLenInfo


#endif /* VERIFYBELLCURVECLAIM_H_ */
