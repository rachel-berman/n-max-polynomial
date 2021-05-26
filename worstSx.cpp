/*
 * worstSx.cpp
 *
 *  Created on: Jul 4, 2017
 *      Author: rberman
 */

#include "worstSx.h"
using namespace std;

static unsigned long int validUxCounter;
static unsigned long int UxOptions;
static bool isOK;

/*
 * ALGORITHM: backtracking (recursive function), with optimizations according to
 * theory results
 * ASSUMPTION: degLenInfo, worstSx, Sx, != NULL, it was checked in "getWorstSx".
 * dealing only with TOTAL histograms
*/
void getSx(DegLengthINFO* degLenInfo, WorstSx* worstSx,int *Sx, int indx,
			int minDeg, int maxDeg);

/* ALGORITHM: backtracking (recursive function)
 * ASSUMPTION: degLenInfo, Sx, Ux != NULL, it was checked in "numValidUxDividingSx".
 * dealing only with REGULAR histograms
 * TODO: change the algorithm and use dynamic programming
*/
void getUx(DegLengthINFO* degLenInfo, int* Sx, int* Ux, int indx);

WorstSx* WorstSxCreate(DegLengthINFO* degLenInfo){
	WorstSx* worstSx = new WorstSx();
	if(worstSx == NULL || degLenInfo == NULL){
		delete worstSx;
		return NULL;
	}
	worstSx->numValidUx = 0;
	worstSx->t = 0;
	worstSx->dt = 0;

	worstSx->sxRegularHist = NULL;
	worstSx->sxRegularHist = new int[degLenInfo->regularHistLength];
	if(worstSx->sxRegularHist == NULL){
		delete worstSx;
		return NULL;
	}
	for(int i=0; i< degLenInfo->regularHistLength; i++){
		worstSx->sxRegularHist[i] = 0;
	}

	worstSx->sxTotalHist = NULL;
	worstSx->sxTotalHist = new int[degLenInfo->totalHistLength];
	if(worstSx->sxTotalHist == NULL){
		delete[] worstSx->sxRegularHist;
		delete worstSx;
		return NULL;
	}
	for(int i=0; i< degLenInfo->totalHistLength; i++){
		worstSx->sxTotalHist[i] = 0;
	}
	 worstSx->t = 0;
	 worstSx->dt = 0;

	return worstSx;
}
void WorstSxDestroy(WorstSx* worstSx){
	delete[] worstSx->sxRegularHist;
	delete[] worstSx->sxTotalHist;
	delete worstSx;
}

STATUS compareToWorstSx(DegLengthINFO* degLenInfo, WorstSx* worstSx,
		int* currentSxTotal){
	if(degLenInfo == NULL || worstSx ==NULL || currentSxTotal == NULL){
		return FAILURE;
	}

	int* currentSxRegular = NULL;
	currentSxRegular = new int[degLenInfo->regularHistLength];
	if(currentSxRegular == NULL){
		return FAILURE;
	}

	int currentNumUx =  getRegularHistAndNumValidUx(degLenInfo, currentSxTotal,
			currentSxRegular);
	if(currentNumUx == -1){
		delete[] currentSxRegular;
		return FAILURE;
	}

	if(currentNumUx > worstSx->numValidUx){
		worstSx->numValidUx = currentNumUx;
		worstSx->dt = getIndexOfRightmostNonZeroElement(currentSxTotal,
											  degLenInfo->totalHistLength)+1;
		worstSx->t = getIndexOfRightmostNonZeroElement(currentSxRegular,
											  degLenInfo->regularHistLength)+1;
		replaceArrayContent(worstSx->sxRegularHist, currentSxRegular,
				degLenInfo->regularHistLength);
	}
	delete[] currentSxRegular;
	return SUCCESS;
}

void getWorstSx(int n, int q, int minDeg, int maxDeg){
	isOK = true;

	DegLengthINFO* degLenInfo =  DegLenInfoCreate( n, q);
	WorstSx* worstSx = WorstSxCreate(degLenInfo);
	if(degLenInfo == NULL || worstSx == NULL){
		DegLenInfoDestroy(degLenInfo);
		WorstSxDestroy(worstSx);
		cout << "problem in alloc degLen and worstSx" << endl;
		return;
	}

	int* Sx = NULL;
	Sx = new int[degLenInfo->totalHistLength];
	if(Sx == NULL){
		DegLenInfoDestroy(degLenInfo);
		WorstSxDestroy(worstSx);
		cout << "problem in alloc Sx" << endl;
		return;
	}
	for(int i=0; i< degLenInfo->totalHistLength; i++){	//initialize Sx to zero
		Sx[i]=0;
	}
	int indx = 0;

	// CALL GET_AND _PRINT_SX
	getSx(degLenInfo, worstSx, Sx, indx, minDeg, maxDeg);

	if(isOK){ // print worstSx
		printWorstSx(degLenInfo, worstSx);
	}else{
		cout << "problem in getSx" << endl;
	}


	DegLenInfoDestroy(degLenInfo);
	WorstSxDestroy(worstSx);
	delete[] Sx;
}

// Sx in total histogram
bool isSxPotentialyWorst(DegLengthINFO* degLenInfo, int* Sx, int minDeg){
	bool minDegreeConstraint =
			degreeTotalHist(Sx,degLenInfo->totalHistLength) >= minDeg;

	int maxIndexNonZeroElement = getIndexOfRightmostNonZeroElement(Sx,
			degLenInfo->totalHistLength);
	if(maxIndexNonZeroElement == -1){
		return false;
	}
	bool isProperty5Satisfied = (int)floor(Sx[0]/degLenInfo->numIrrPolyPerDeg[0])
			>= maxIndexNonZeroElement;

	int minDegNotInSx = 0;
	if((int)floor(Sx[maxIndexNonZeroElement]/
			degLenInfo->numIrrPolyPerDeg[maxIndexNonZeroElement])==0 ){
		minDegNotInSx = maxIndexNonZeroElement;
	}else{
		minDegNotInSx = maxIndexNonZeroElement + 1;
	}

	bool isProperty4Satisfied = (int)ceil(Sx[0]/degLenInfo->numIrrPolyPerDeg[0])
			< 2*minDegNotInSx;

	//bool isProperty5Satisfied = true;
	//bool isProperty4Satisfied = true;
	return minDegreeConstraint && isProperty5Satisfied && isProperty4Satisfied;
}

void getSx(DegLengthINFO* degLenInfo, WorstSx* worstSx,int *Sx, int indx,
		int minDeg, int maxDeg){
	if(indx >= degLenInfo->totalHistLength){ // we have s(x) in hand

		if( isSxPotentialyWorst(degLenInfo, Sx, minDeg)
				/*degreeTotalHist(Sx,degLenInfo->totalHistLength) >= minDeg*/ ){
			STATUS status = compareToWorstSx(degLenInfo, worstSx, Sx);
			if(status != SUCCESS){
				isOK = false;
				cout << "problem in getSx!" << endl;
			}

		}
		return;
	}

	int Sx_indxOriginal = Sx[indx];

	bool condition = false;
	do{
		getSx(degLenInfo, worstSx,Sx, indx+1, minDeg, maxDeg);

		++Sx[indx];
		//set condition:
		if(indx == 0){
			condition = degreeTotalHist(Sx, indx+1) <= maxDeg;
		}else{
			condition = degreeTotalHist(Sx, indx+1) <= maxDeg &&
					Sx[indx-1] >= degLenInfo->numIrrPolyPerDeg[indx-1] && /*property 1*/
					(int)floor( Sx[indx-1]/degLenInfo->numIrrPolyPerDeg[indx-1] ) /* property 3 */
					>= (int)ceil( Sx[indx]/degLenInfo->numIrrPolyPerDeg[indx] );
		}

	}while(condition);
	Sx[indx] = Sx_indxOriginal;
}

void printWorstSx(DegLengthINFO* degLenInfo, WorstSx* worstSx){
	const char separator    = ' ';
	const int nWidth		= 6;
	const int degSxWidth	= 8;
	const int Lc1Width		= 17;
	const int tWidth 		= 6;
	const int dtWidth		= 6;
	const int r0Width 		= 6;
	const int riWidth 		= 4;

	int degSx = getInnerProduct(worstSx->sxRegularHist,
			degLenInfo->degVecForRegularHist, degLenInfo->regularHistLength);
	int r0 = 2*degLenInfo->n - degSx;
	globalDeg_sx = degSx; // TODO use global variable for purpose degSn(x) >= degS(n-1)(x)

	if(degLenInfo->n == 2){ // print the table header row
		printElement("n", nWidth, separator);
		printElement("degS(x)", degSxWidth, separator);
		printElement("Lc(1)", Lc1Width, separator);
		printElement("t", tWidth, separator);
		printElement("dt", dtWidth, separator);
		cout << "| ";
		printElement("r0", r0Width, separator);
		cout << "S(x)" << endl;
	}

	printElement(degLenInfo->n, nWidth, separator);
	printElement(degSx, degSxWidth, separator);
	printElement(worstSx->numValidUx, Lc1Width, separator);
	printElement(worstSx->t, tWidth, separator);
	printElement(worstSx->dt, dtWidth, separator);
	cout << "| ";
	printElement(r0, r0Width, separator);

	for(int i=0; i<degLenInfo->regularHistLength; i++){
		if(worstSx->sxRegularHist[i] == 0)
			break;

		printElement(worstSx->sxRegularHist[i], riWidth, separator);
	}
	cout << endl;
}


// dealing only with REGULAR histograms
int numValidUxDividingSx(DegLengthINFO* degLenInfo, int* Sx){
	if(degLenInfo==NULL || Sx==NULL){
		return -1;
	}
	validUxCounter = 0;
	UxOptions = 0; // TODO erase
	//initialize parameters for the recursive function getUx:
	int indx = 0;
	int* Ux = new int[degLenInfo->regularHistLength];
	if(Ux == NULL){
		return -1;
	}
	for(int i=0; i< degLenInfo->regularHistLength; i++){
		Ux[i] = 0;
	}

	getUx(degLenInfo, Sx, Ux, indx);

	delete[] Ux;
	Ux = NULL;

	//cout << "calls of getUx = "<< UxOptions << endl; // TODO erase
	return validUxCounter;
}


void getUx(DegLengthINFO* degLenInfo, int* Sx, int* Ux, int indx){
	if(indx >= degLenInfo->regularHistLength){
		//printArray(Ux, length); // TODO: for DEBUG, delete this line
		int degUx = getInnerProduct(Ux, degLenInfo->degVecForRegularHist,
				degLenInfo->regularHistLength);
		int degVx = getInnerProductWithDifference(Sx, Ux,
				degLenInfo->degVecForRegularHist, degLenInfo->regularHistLength);
		if(degUx <= degLenInfo->n && degVx <= degLenInfo->n){
			validUxCounter++;
			//printArray(Ux, length); // TODO: erase
		}
		return;
	}

	int Ux_IndxOriginal = Ux[indx];
	do{
		//printArray(Ux, indx+1); // TODO erase
		getUx(degLenInfo, Sx, Ux, indx+1);
		++UxOptions; // TODO erase
		++Ux[indx];
	}while(Ux[indx] <= Sx[indx] &&
			getInnerProduct(Ux, degLenInfo->degVecForRegularHist, indx+1) <= degLenInfo->n
			/*&&
			getInnerProductWithDifference(Sx, Ux, degOfPoly, indx+1) <= n */);

	Ux[indx] = Ux_IndxOriginal;
}


int getRegularHistAndNumValidUx(DegLengthINFO* degLenInfo, int* sTotalHist,
		 int* sRegularHist){

	if(degLenInfo == NULL || sTotalHist == NULL ||
			 sRegularHist == NULL){
		return -1;
	}
	STATUS status = getRegularHistFromTotalHist(degLenInfo, sTotalHist,
			sRegularHist);
	if(status != SUCCESS){
		return -1;
	}

	return numValidUxDividingSx(degLenInfo, sRegularHist);
}


/****  ADVANCED ALGORITHM, ACCORDING TO CHARACTERIZATION IN THEOREM 2.4  ****/
/* Theorem 2.4: total multiplicity of d-factors = alpha*ro/d + beta*(r0/d - 1),
 * where 0<=alpha<=L2(d) */

void calculateSxTotalHistFromAlpha(DegLengthINFO* degLenInfo,int* alpha, int* Sx, int r0);
int calculateDegFromAlpha(DegLengthINFO* degLenInfo,int* alpha, int length, int r0);
/*
 * ALGORITHM: backtracking on alpha (recursive function)
 * ASSUMPTION: degLenInfo, worstSx, Sx, alpha, != NULL, it was checked in
 * "getWorstSxAdvanced".
 * dealing only with TOTAL histograms
*/
void getSxAdvanced(DegLengthINFO* degLenInfo, WorstSx* worstSx, int *Sx,
		int* alpha, int indx, int r0, bool printPossibleWorstPoly);
STATUS compareToWorstSxAdvanced(DegLengthINFO* degLen, int r0, int* currentSxTotal,
							    WorstSx* worstSx);

unsigned long int countValidUxTotalVersion(DegLengthINFO* degLen, int r0, int* SxTotal, int* SxRegular);

/* ALGORITHM: backtracking (recursive function)
 * ASSUMPTION: degLenInfo, SxTotal, SxRegular, UxTotal != NULL,
 * it was checked in "numValidUxDividingSx".
*/
void getUxTotal(DegLengthINFO* degLen, int r0, int* SxTotal, int* SxRegular,
				int* UxTotal, int indx);

unsigned long int countConfigurationsOfUxInRegularHist(DegLengthINFO* degLen, int r0,
									int* SxTotal, int* SxRegular, int* UxTotal);


void getWorstSxAdvanced(int n, int q, bool printPossibleWorstPoly){
	isOK = true;

	DegLengthINFO* degLenInfo =  DegLenInfoCreate( n, q);
	WorstSx* worstSx = WorstSxCreate(degLenInfo);
	if(degLenInfo == NULL || worstSx == NULL){
		DegLenInfoDestroy(degLenInfo);
		WorstSxDestroy(worstSx);
		cout << "problem in alloc degLen and worstSx" << endl;
		return;
	}

	int* Sx = NULL;
	Sx = new int[degLenInfo->totalHistLength];
	if(Sx == NULL){
		DegLenInfoDestroy(degLenInfo);
		WorstSxDestroy(worstSx);
		cout << "problem in alloc Sx" << endl;
		return;
	}

	int* alpha = NULL;
	alpha = new int[degLenInfo->totalHistLength];
	if(alpha == NULL){
		DegLenInfoDestroy(degLenInfo);
		WorstSxDestroy(worstSx);
		delete[] Sx;
		cout << "problem in alloc alpha" << endl;
		return;
	}
	for(int i=0; i< degLenInfo->totalHistLength; i++){	//initialize Sx to zero
		Sx[i]=0;
		alpha[i]=0;
	}
	int indx = 0;

	// CALL GET_AND _PRINT_SX
	int r0LowerBound, r0UpperBound;
	//r0LowerBound= (int)ceil(log2(n)-log2(log2(n)));
	//r0UpperBound= (int)floor(2*log2(n)+3);
	r0LowerBound = 4; // instead of (int)ceil(log2(n)-log2(log2(n))); //TODO
	r0UpperBound = 13; // instead of (int)floor(2*log2(n)+2) //TODO
	for(int r0 = r0LowerBound; r0<= r0UpperBound; r0++){
		getSxAdvanced(degLenInfo, worstSx, Sx, alpha, indx, r0, printPossibleWorstPoly);
		if(!isOK){
			cout << "problem in getSxAdvanced" << endl;
			return;
		}
	}

	if(isOK){ // print worstSx
		printWorstSx(degLenInfo, worstSx);
	}else{
		cout << "problem in getSxAdvanced" << endl;
	}


	DegLenInfoDestroy(degLenInfo);
	WorstSxDestroy(worstSx);
	delete[] Sx;
	delete[] alpha;

}

void getSxAdvanced(DegLengthINFO* degLenInfo, WorstSx* worstSx, int *Sx,
		int* alpha, int indx, int r0, bool printPossibleWorstPoly){
	if(indx >= degLenInfo->totalHistLength){ // we have alpha in hand
		calculateSxTotalHistFromAlpha(degLenInfo,alpha, Sx, r0);
		if(r0 + degreeTotalHist(Sx, degLenInfo->totalHistLength) ==
				2*degLenInfo->n){

			for(int i=0; i<degLenInfo->totalHistLength-1; i++){
				if(Sx[i]<degLenInfo->numIrrPolyPerDeg[i] && Sx[i+1]>0){ //there is zero-hole
					return;
				}
			}
			// we have Sx (in total histogram) that satisfies theorem 2.4
			if(printPossibleWorstPoly){
				// total-to-regular
				int* SxRegular = NULL;
				SxRegular = new int[degLenInfo->regularHistLength];
				if(SxRegular == NULL){
					cout << "oy vey"<< endl;
					return;
				}
				STATUS status = getRegularHistFromTotalHist(degLenInfo, Sx,
						SxRegular);
				if(status != SUCCESS){
					cout << "oy vey2"<< endl;
					delete[] SxRegular;
					return;
				}
				// print regular
				cout << r0 << "  ";
				printArray(SxRegular, degLenInfo->regularHistLength);
				delete[] SxRegular;
				return;
			}else{
/*
				printArray(alpha, degLenInfo->totalHistLength); // DEBUG
				cout << r0 << " | ";  // DEBUG
				printArray(Sx, degLenInfo->totalHistLength); // DEBUG
*/

				/* GIVEN sx in total histogram
					   check how much valid factorizations
					   compare to current-worst Sx */
				// new algorithm - use total histograms (instead of regular histograms) for counting valid ux
				STATUS status = compareToWorstSxAdvanced(degLenInfo, r0, Sx, worstSx);

				// old algorithm - use regular histograms for counting valid ux
				//STATUS status = compareToWorstSx(degLenInfo, worstSx, Sx);


				if(status != SUCCESS){
					isOK = false;
					cout << "problem in getSxAdvanced!" << endl;
				}
			}

		}
		return;
	}

	int alpha_indxOriginal = alpha[indx];
	do{
		getSxAdvanced(degLenInfo, worstSx, Sx, alpha, indx+1, r0, printPossibleWorstPoly);
		++alpha[indx];
	}while(alpha[indx] <= degLenInfo->numIrrPolyPerDeg[indx] &&
			r0 + calculateDegFromAlpha(degLenInfo,alpha,indx+1,r0)<= 2*degLenInfo->n &&
			!(r0/(indx+1)==0 && alpha[indx]>0));

	alpha[indx]=alpha_indxOriginal;
}

STATUS compareToWorstSxAdvanced(DegLengthINFO* degLen,int r0, int* currentSxTotal,
							    WorstSx* worstSx){
	int* currentSxRegular = NULL;
	currentSxRegular = new int[degLen->regularHistLength];
	if(currentSxRegular==NULL) return FAILURE;
	STATUS status = getRegularHistFromTotalHist(degLen, currentSxTotal,currentSxRegular);
	if(status==FAILURE) return FAILURE;

	unsigned long int currentValidUxCounter = countValidUxTotalVersion(degLen, r0, currentSxTotal, currentSxRegular);
	if(!(currentValidUxCounter > 0)) return FAILURE;
	if(currentValidUxCounter > worstSx->numValidUx){
		replaceArrayContent(worstSx->sxTotalHist, currentSxTotal, degLen->totalHistLength);
		replaceArrayContent(worstSx->sxRegularHist, currentSxRegular, degLen->regularHistLength);
		worstSx->numValidUx = currentValidUxCounter;
		worstSx->dt = getIndexOfRightmostNonZeroElement(currentSxTotal,
											  degLen->totalHistLength)+1;
		worstSx->t = getIndexOfRightmostNonZeroElement(currentSxRegular,
											  degLen->regularHistLength)+1;
/*//START DEBUG
		cout << "----------------------------------------------------" << endl<< endl;
		printWorstSx(degLen, worstSx);
	}
	if(currentValidUxCounter == worstSx->numValidUx){
		replaceArrayContent(worstSx->sxTotalHist, currentSxTotal, degLen->totalHistLength);
		replaceArrayContent(worstSx->sxRegularHist, currentSxRegular, degLen->regularHistLength);
		worstSx->numValidUx = currentValidUxCounter;
		worstSx->dt = getIndexOfRightmostNonZeroElement(currentSxTotal,
				degLen->totalHistLength)+1;
		worstSx->t = getIndexOfRightmostNonZeroElement(currentSxRegular,
				degLen->regularHistLength)+1;


		printWorstSx(degLen, worstSx);
// END DEBUG */
	}

	delete[] currentSxRegular;
	return SUCCESS;
}

unsigned long int countValidUxTotalVersion(DegLengthINFO* degLen, int r0, int* SxTotal, int* SxRegular){
	int* UxTotal = NULL;
	UxTotal = new int[degLen->totalHistLength];
	if(UxTotal==NULL){
		cout << "count valid ux total - problem";
		return 0;
	}
	for(int i=0; i<degLen->totalHistLength; i++){
		UxTotal[i] = 0;
	}
	validUxCounter = 0;
	int indx = 0;
	getUxTotal(degLen, r0, SxTotal, SxRegular, UxTotal, indx);

	delete[] UxTotal;
	return validUxCounter;
}

void getUxTotal(DegLengthINFO* degLen, int r0, int* SxTotal, int* SxRegular,
				int* UxTotal, int indx){
	if(indx >= degLen->totalHistLength){
		int degUx = degreeTotalHist(UxTotal,degLen->totalHistLength);
		int degVx = 0;
		for(int i=0; i < degLen->totalHistLength; i++){
			degVx += (i+1)*(SxTotal[i]-UxTotal[i]);
		}
		if(degUx <= degLen->n && degVx <= degLen->n){
			validUxCounter += countConfigurationsOfUxInRegularHist(degLen, r0,
									SxTotal,SxRegular,UxTotal);
		}
		return;
	}
	int UxTotal_indxOriginal = UxTotal[indx];
	do{
		getUxTotal(degLen, r0, SxTotal, SxRegular, UxTotal, indx+1);
		++UxTotal[indx];
	}while(UxTotal[indx]<=SxTotal[indx] &&
			degreeTotalHist(UxTotal, indx+1) <= degLen->n);
	UxTotal[indx] = UxTotal_indxOriginal;
}



unsigned long int countConfigurationsOfUxInRegularHist(DegLengthINFO* degLen, int r0,
									int* SxTotal, int* SxRegular, int* UxTotal){
	unsigned long int configurationsCounter = 1;
	for(int dIndx = 0; dIndx<degLen->totalHistLength; dIndx++){
		if(dIndx+1 > r0/2){
			configurationsCounter *= nChoosek(SxTotal[dIndx],UxTotal[dIndx]);
		}else{
			configurationsCounter *= getNumSubvectorsOfGivenWeight(
							SxRegular+getIndexOfDegRegularHist(degLen,dIndx),
							degLen->numIrrPolyPerDeg[dIndx],
							UxTotal[dIndx]);
		}
	}
	return configurationsCounter;
}


/* assumption: pointers are not NULL */
void calculateSxTotalHistFromAlpha(DegLengthINFO* degLenInfo,int* alpha, int* Sx, int r0){
	for(int i=0; i<degLenInfo->totalHistLength; i++){
		Sx[i] = alpha[i]*(r0/(i+1)) +
				(degLenInfo->numIrrPolyPerDeg[i]-alpha[i])*maxWith0(r0/(i+1)-1);
	}

}

/* assumption: pointers are not NULL */
int calculateDegFromAlpha(DegLengthINFO* degLenInfo,int* alpha, int length, int r0){
	int sum = 0;
	for(int i=0; i< length; i++){ // degree = i+1
		int totalMultiplicity = alpha[i]*(r0/(i+1)) +
				(degLenInfo->numIrrPolyPerDeg[i]-alpha[i])*maxWith0(r0/(i+1)-1);
		sum += (i+1)*totalMultiplicity;
	}
	return sum;
}
