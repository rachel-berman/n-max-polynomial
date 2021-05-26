/*
 * verifyBellCurveClaim.cpp
 *
 *  Created on: Oct 1, 2017
 *      Author: rberman
 */

#include "verifyBellCurveClaim.h"
using namespace std;

#define ALL_IRR_LINEAR 0
#define PRINT_FOR_EXCEL 0


struct SxDetails
{
	int* SxHist; /* pointer to array of length regularHistLength */
	int Sxlength;
	int *sizeBk; /* length of array = 2*n+1 */
	int sizeBkLength;
};

SxDetails* SxDetailsCreate(int n, int *SxHist, int SxHistLength);
void SxDetailsDestroy(SxDetails* sxDetails);

// ALGORITHM: backtracking (recursive function)
// ASSUMPTION: degOfPoly, sxDetails, Ux != NULL, it was checked in "printSizeBkCurve".
void getAndCountUx(ComplementaryDegLengthINFO* degLenInfo, SxDetails* sxDetails,
					int* Ux, int indx);

SxDetails* SxDetailsCreate(int n, int *SxHist, int SxHistLength){
	SxDetails* sxDetails = new SxDetails();
	if(!sxDetails || !SxHist){
		delete sxDetails;
		return NULL;
	}
	sxDetails->Sxlength = SxHistLength;
	sxDetails->SxHist = NULL;
	sxDetails->sizeBkLength = 2*n+1;
	sxDetails->sizeBk = NULL;

	sxDetails->SxHist = new int[sxDetails->Sxlength];
	sxDetails->sizeBk = new int[sxDetails->sizeBkLength];
	if(!sxDetails->SxHist || !sxDetails->sizeBk){
		delete[] sxDetails->SxHist;
		delete[] sxDetails->sizeBk;
		delete sxDetails;
		return NULL;
	}
	for(int i=0; i < sxDetails->Sxlength; i++){
		sxDetails->SxHist[i] = SxHist[i];
	}
	for(int i=0; i < sxDetails->sizeBkLength; i++){
		sxDetails->sizeBk[i] = 0;
	}
	return sxDetails;
}

void SxDetailsDestroy(SxDetails* sxDetails){
	if(sxDetails){
		delete[] sxDetails->SxHist;
		delete[] sxDetails->sizeBk;
		delete sxDetails;
	}
}

/* PARSING 'RESULT.TXT' FILE - DEFINITIONS */
SxInfo* SxInfoCreate(){
	SxInfo* sxInfo = new SxInfo();
	if(!sxInfo){
		cout << "'SxInfoCreate' - problem in alloc sxInfo"<< endl;
		return NULL;
	}
	sxInfo->sxHist = new int[MAX_HIST_LENGTH];
	if(!sxInfo->sxHist){
		delete sxInfo;
		cout << "'SxInfoCreate' - problem in alloc sxInfo->sxHist"<< endl;
		return NULL;
	}

	// Initializing:
	sxInfo->n = -1;
	sxInfo->Lc = -1;
	sxInfo->degSx = -1;
	sxInfo->totalLength = 0;
	for(int i=0; i<MAX_HIST_LENGTH; i++){
		sxInfo->sxHist[i] = 0;
	}

	return sxInfo;
}

void SxInfoDestroy(SxInfo* sxInfo){
	delete[] sxInfo->sxHist;
	delete sxInfo;
}

SxInfo* parsingSxLine(string line){
	istringstream iss(line);

	SxInfo* sxInfo = SxInfoCreate();
	if(!sxInfo){
		return NULL;
	}

	iss >> sxInfo->n;
	iss >> sxInfo->degSx;
	iss >> sxInfo->Lc;

	sxInfo->sxHist[0] = 2*(sxInfo->n) - (sxInfo->degSx); // calculating 'r0'
	++sxInfo->totalLength;
	int i=1;
	while(iss >> sxInfo->sxHist[i]){
		++sxInfo->totalLength;
		++i;
	}

	return sxInfo;
}

void printSxInfo(SxInfo* sxInfo){
	if(!sxInfo){
		cout << "'sxInfo' - NULL pointer" << endl;
		return;
	}
	cout << "n = " << sxInfo->n << endl;
	cout << "degS = " << sxInfo->degSx << endl;
	cout << "Lc = " << sxInfo->Lc << endl;
	cout << "total length = " << sxInfo->totalLength << endl;
	printArray(sxInfo->sxHist, sxInfo->totalLength);
	cout << endl;
}
/* PARSING - END OF DEFINITIONS */

void printSizeBkCurve(int n, int q, int *SxHist, int SxHistLength, int Lc){
	ComplementaryDegLengthINFO* degLenInfo =  ComplementaryDegLenInfoCreate( n, q);
	SxDetails* sxDetails = SxDetailsCreate(n, SxHist, SxHistLength);
	if(degLenInfo == NULL || sxDetails == NULL){
		ComplementaryDegLenInfoDestroy(degLenInfo);
		SxDetailsDestroy(sxDetails);
		cout << "'printSizeBkCurve' - problem in alloc degLen and sxDetails" << endl;
		return;
	}

	if(SxHistLength > degLenInfo->regularHistLength){
		ComplementaryDegLenInfoDestroy(degLenInfo);
		SxDetailsDestroy(sxDetails);
		cout << "'printSizeBkCurve' - length of s(x) is too long!" << endl;
		return;
	}

	if(getInnerProduct(SxHist, degLenInfo->degVecForRegularHist, SxHistLength) > 2*n){
		ComplementaryDegLenInfoDestroy(degLenInfo);
		SxDetailsDestroy(sxDetails);
		cout << "'printSizeBkCurve' - deg s(x) > 2*n" << endl;
		return;
	}


	// call getAndCountUx
	int indx = 0;
	int* Ux = new int[sxDetails->Sxlength];
	if(Ux == NULL){
		cout << "'printSizeBkCurve' - problem in alloc Ux" << endl;
		ComplementaryDegLenInfoDestroy(degLenInfo);
		SxDetailsDestroy(sxDetails);
		return;
	}
	for(int i=0; i< sxDetails->Sxlength; i++){
		Ux[i] = 0;
	}
	getAndCountUx(degLenInfo, sxDetails, Ux, indx);

	// print results
	if(PRINT_FOR_EXCEL){
		cout << n << " ";
		printCurve(sxDetails->sizeBk, sxDetails->sizeBkLength);
/*
		cout << n << endl;
		printCurve(sxDetails->sizeBk, sxDetails->sizeBkLength);
		cout <<	(int)averageValueArray(sxDetails->sizeBk, sxDetails->sizeBkLength)<< endl;
*/
	}else{
		cout << "n = " << n <<" , Histogram: " << endl;
		printArray(sxDetails->SxHist, sxDetails->Sxlength);

		cout << "Lc = | Bn| --- ";
		(Lc == sxDetails->sizeBk[n]) ? (cout << "YES") : (cout << "NO");
		cout << endl;
		cout << "average size Bk = " <<
				(int)averageValueArray(sxDetails->sizeBk, sxDetails->sizeBkLength)<< endl;
		printCurve(sxDetails->sizeBk, sxDetails->sizeBkLength);
		cout << endl;
	}

	// free memory
	ComplementaryDegLenInfoDestroy(degLenInfo);
	SxDetailsDestroy(sxDetails);
	delete[] Ux;
}


void getAndCountUx(ComplementaryDegLengthINFO* degLenInfo, SxDetails* sxDetails,
					int* Ux, int indx){
	if(indx >= sxDetails->Sxlength){
		int degUx = 0;
		if(ALL_IRR_LINEAR){
			for(int i=0; i< sxDetails->Sxlength; i++){
				degUx += Ux[i];
			}
		}else{
			degUx = getInnerProduct(Ux, degLenInfo->degVecForRegularHist,
												sxDetails->Sxlength);
		}

		sxDetails->sizeBk[degUx]++;
		//printArray(Ux, length); // TODO: for DEBUG, delete this line
		return;
	}

	int Ux_IndxOriginal = Ux[indx];
	do{
		//printArray(Ux, indx+1); // TODO erase
		getAndCountUx(degLenInfo, sxDetails, Ux, indx+1);
		++Ux[indx];
	}while(Ux[indx] <= sxDetails->SxHist[indx]);

	Ux[indx] = Ux_IndxOriginal;
}

void printSiCurve(SxInfo* sxInfo){
	// void printSizeBkCurve(int n, int q, int *SxHist, int SxHistLength, int Lc){
	int degVec[24] = {1,1,1,2,3,3,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6};
	int lengthSx = sxInfo->totalLength;
/*  N=14 **************
	int lengthSx = 7;
	SxInfo* sxInfo = SxInfoCreate();
	sxInfo->n = 14;
	sxInfo->totalLength = lengthSx;
	sxInfo->degSx = 23;
	sxInfo->Lc = 382;
	int tmp[sxInfo->totalLength] = {5,5,4,2,1,1,1};
	for(int i=0; i<sxInfo->totalLength; ++i){
		sxInfo->sxHist[i]=tmp[i];
	}
*/
	//SxInfo* sxInfo_i = SxInfoCreate();

	// N=4 **********
	/*int lengthSx = 4;
	SxInfo* sxInfo = SxInfoCreate();
	sxInfo->n = 4;
	sxInfo->totalLength = lengthSx;
	sxInfo->degSx = 6;
	sxInfo->Lc = 12;
	int tmp[sxInfo->totalLength] = {2,2,2,1};
	for(int i=0; i<sxInfo->totalLength; ++i){
		sxInfo->sxHist[i]=tmp[i];
	}*/

	for(int i=0; i < lengthSx; ++i){
		//set Si(x) info
		sxInfo->totalLength = i+1;
		sxInfo->degSx = getInnerProduct(sxInfo->sxHist+1,degVec+1, sxInfo->totalLength ); // skipping r0

		//print Bk-curve of Si
		printSizeBkCurve(sxInfo->n, Q, sxInfo->sxHist, i+1, 1);
	}


	//SxInfoDestroy(sxInfo);
	//SxInfoDestroy(sxInfo_i);
}

// ASSUMPTION: THE LENGTH OF ARRAY "curveValues" IS ODD
void printCurve(int *curveValues, int len){
	int middle = len/2;

	if(PRINT_FOR_EXCEL){
		for(int i=0; i<len; i++){
			cout << curveValues[i] << " ";
		}
	}else{ // PRINT ARRAY WITH EMPHASIZING THE MIDDLE VALUE
		for(int i=0; i<len; i++){
			if(i==middle)
				cout << "[";
			cout << curveValues[i];
			if(i==middle)
				cout << "]";
			cout << " ";
		}
	}
	cout << endl;
}



