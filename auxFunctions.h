/*
 * auxFunctions.h
 *
 *  Created on: Jul 4, 2017
 *      Author: rberman
 */

#ifndef AUXFUNCTIONS_H_
#define AUXFUNCTIONS_H_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <new>
#include <cmath>

#define Q 2 // the size of the Galois field

using namespace std;


typedef enum status{
	SUCCESS, FAILURE, NULL_PTR, MISMATCHING_ARAAYS_LENGTHS
}STATUS;

struct DegLengthINFO
{
	int n;
	int q;
	int totalHistLength;
	int regularHistLength;
	int* numIrrPolyPerDeg; /* pointer to array of length totalHistLength */
	int* degVecForRegularHist; /* pointer to array of length regularHistLength */
};

DegLengthINFO* DegLenInfoCreate(int n, int q);
void DegLenInfoDestroy(DegLengthINFO*);

struct ComplementaryDegLengthINFO
{
	int n;
	int q;
	int regularHistLength;
	int* degVecForRegularHist; /* pointer to array of length regularHistLength */
};

ComplementaryDegLengthINFO* ComplementaryDegLenInfoCreate(int n, int q);
void ComplementaryDegLenInfoDestroy(ComplementaryDegLengthINFO*);


unsigned nChoosek( unsigned n, unsigned k );


//mu(i) = mobius(i). Ignore 0-coordinate
void getMobius(int* mu, int len);

/* return the number of irreducible polynomials over GF(q) of degree n,
   according to the formula with Mobius function. */
int numIrrPolyOfDegN(int n, int q);

/* according to formulas Lq(n), Dq(n) */
STATUS getSxHistLength(int n, int q, int* lengthTotalHist, int* lengthRegularHist);

STATUS getDegVecForRegularHist(int* numIrrPolyPerDeg, int maxDeg,
		int* degOfPolyInRegulatHist, int lengthRegular);

/* calculate the regular histogram of s(x) from the total histogram
 *  (according to property 2), write it to the output parameter "regularHist" */
STATUS getRegularHistFromTotalHist(DegLengthINFO* degLenInfo, int* totalHist,
		int* regulatHist);

// ASSUMPTION: regularArray != NULL
void distributeEvenlyMultiplicity(int* regularArray, int length, int totalMultiplicity);

// ASSUMPTION: 0<= d < totalHistLength, degLen != NULL
int getIndexOfDegRegularHist(DegLengthINFO* degLen, int d);

// ASSUMPTION: weight <= weight(vector)
unsigned long int getNumSubvectorsOfGivenWeight(int* vector, int length, int weight);
int getVectorWeight(int* vector, int length);

int getInnerProduct(int* a, int* b, int length);
int getInnerProductWithDifference(int* a1, int* a2, int* b, int length);

// ASSUMPTION: ptr != NULL, and the arrays are of the same length
void replaceArrayContent(int* array, int* newContent, int length);

int degreeTotalHist(int* totalHist, int lengthTotalHist);

//ASSUMPTIONs: array is not zero-vector, with "no holes of zeros".
// array != NULL
int getIndexOfRightmostNonZeroElement(int* array, int length);

double averageValueArray(int *array, int length);

int maxWith0(int m);

void printArray(int* arr, int len);

void printHist(int* totalHist, int lengthTotal, int* regulatHist,
		int lengthRegular);

void printDegLenStruct(DegLengthINFO*);


template<typename T> void printElement(T t, const int& width, const char& separator)
{
    cout << left << setw(width) << setfill(separator) << t;
}

#endif /* AUXFUNCTIONS_H_ */
