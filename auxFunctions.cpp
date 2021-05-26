/*
 * auxFunctions.cpp
 *
 *  Created on: Jul 4, 2017
 *      Author: rberman
 */


#include "auxFunctions.h"

using namespace std;

static unsigned long int subvectorCounter;

DegLengthINFO* DegLenInfoCreate(int n, int q){
	DegLengthINFO* degLen = new DegLengthINFO();
	if(degLen == NULL){
		return NULL;
	}

	degLen->n = n;
	degLen->q =q;
	STATUS status = getSxHistLength(n,q,
			& degLen->totalHistLength, & degLen->regularHistLength);
	if(status != SUCCESS){
		delete degLen;
		return NULL;
	}

	degLen->numIrrPolyPerDeg = NULL;
	degLen->degVecForRegularHist = NULL;
	degLen->numIrrPolyPerDeg = new int[degLen->totalHistLength];
	degLen->degVecForRegularHist = new int[degLen->regularHistLength];
	if(degLen->numIrrPolyPerDeg == NULL || degLen->degVecForRegularHist==NULL){
		delete[] degLen->numIrrPolyPerDeg;
		delete[] degLen->degVecForRegularHist;
		delete degLen;
		return NULL;
	}

	for(int d=0; d < degLen->totalHistLength; d++){
		degLen->numIrrPolyPerDeg[d] = numIrrPolyOfDegN(d+1, q);
	}

	status = getDegVecForRegularHist(degLen->numIrrPolyPerDeg, degLen->totalHistLength,
			degLen->degVecForRegularHist, degLen->regularHistLength);
	/* TODO check "status" if there are any problems. */
	return degLen;
}

void DegLenInfoDestroy(DegLengthINFO* degLen){
	delete[] degLen->numIrrPolyPerDeg;
	delete[] degLen->degVecForRegularHist;
	delete degLen;
}


ComplementaryDegLengthINFO* ComplementaryDegLenInfoCreate(int n, int q){
	ComplementaryDegLengthINFO* complementaryDegLen = new ComplementaryDegLengthINFO();
	if(!complementaryDegLen){
		return NULL;
	}

	DegLengthINFO* degLen = DegLenInfoCreate(n,q);
	if(!degLen){
		delete complementaryDegLen;
		return NULL;
	}

	complementaryDegLen->n = n;
	complementaryDegLen->q =q;
	// add coordinate for multiplicity r_0 of linear factor 'y'
	complementaryDegLen->regularHistLength = degLen->regularHistLength + 1;
	complementaryDegLen->degVecForRegularHist = NULL;
	complementaryDegLen->degVecForRegularHist =
							new int[complementaryDegLen->regularHistLength];
	if(!complementaryDegLen->degVecForRegularHist){
		delete complementaryDegLen;
		DegLenInfoDestroy(degLen);
		return NULL;
	}

	complementaryDegLen->degVecForRegularHist[0] = 1; // degree of 'y'
	for(int i=0; i< degLen->regularHistLength; i++){
		complementaryDegLen->degVecForRegularHist[i+1] =
											degLen->degVecForRegularHist[i];
	}

	DegLenInfoDestroy(degLen);
	return complementaryDegLen;
}

void ComplementaryDegLenInfoDestroy(ComplementaryDegLengthINFO* degLen){
	delete[] degLen->degVecForRegularHist;
	delete degLen;
}

unsigned nChoosek( unsigned n, unsigned k ){
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( unsigned i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}


void getMobius(int* mu, int len){
	if(mu==NULL){
		cout << "getMobius NULL ptr"<<endl;
		return;
	}
	int max = len-1;
	int sqrtMax = (int)floor((double)sqrt(max));
	for(int i = 1; i <= max; i++)
		mu[i] = 1;

	for (int i = 2; i <= sqrtMax; i++){
		if (mu[i] == 1)
		{
			for (int j = i; j <= max; j += i)
				mu[j] *= -i;
			for (int j = i * i; j <= max; j += i * i)
				mu[j] = 0;
		}
	}

	for (int i = 2; i <= max; i++)
	{
		if (mu[i] == i)
			mu[i] = 1;
		else if (mu[i] == -i)
			mu[i] = -1;
		else if (mu[i] < 0)
			mu[i] = 1;
		else if (mu[i] > 0)
			mu[i] = -1;
	}

}

int numIrrPolyOfDegN(int n, int q){
	int* mu = new int[n+1];
	getMobius(mu, n+1);

	int sum = 0, qPower = 1;

	for(int d=1; d <= n; d++){
		qPower *= q;
		if(n%d==0){
			sum += mu[n/d] * qPower;
		}
	}
	delete[] mu;
	mu = NULL;
	return sum/n;
}


STATUS getSxHistLength(int n, int q, int* lengthTotalHist, int* lengthRegularHist){
	if(lengthTotalHist==NULL || lengthRegularHist==NULL){
		return NULL_PTR;
	}
	int tmpTotalDegree = 0, tmpLengthRegularHist=0;

	for(int d=1; d <= n; d++){
		int numIrrPolyPerDeg = numIrrPolyOfDegN(d,q);
		tmpLengthRegularHist += numIrrPolyPerDeg;
		tmpTotalDegree += d*numIrrPolyPerDeg;
		if(tmpTotalDegree >= 2*n){ // d = Dq(n) = length of total Histogram
			*lengthTotalHist = d;
			*lengthRegularHist = tmpLengthRegularHist;
			return SUCCESS;
		}
	}

	return FAILURE;
}

STATUS getDegVecForRegularHist(int* numIrrPolyPerDeg, int maxDeg,
		int* degOfPolyInRegulatHist, int lengthRegular){
	if(numIrrPolyPerDeg==NULL || degOfPolyInRegulatHist==NULL){
		return NULL_PTR;
	}

	int regularIndx = 0;
	for(int d = 1; d <= maxDeg; d++){
		int lengthPerDeg = numIrrPolyPerDeg[d-1];
		for(int i=0; i < lengthPerDeg; i++){
			degOfPolyInRegulatHist[regularIndx + i] = d;
		}
		regularIndx += lengthPerDeg;
	}

	if(regularIndx > lengthRegular){
		return MISMATCHING_ARAAYS_LENGTHS;
	}

	return SUCCESS;
}


STATUS getRegularHistFromTotalHist(DegLengthINFO* degLenInfo, int* totalHist,
		int* regulatHist){

	if(degLenInfo->numIrrPolyPerDeg==NULL || totalHist==NULL || regulatHist==NULL){
		return NULL_PTR;
	}

	int regularIndx = 0;
	for(int d = 1; d <= degLenInfo->totalHistLength; d++){
		int lengthPerDeg = degLenInfo->numIrrPolyPerDeg[d-1];
		int multiplicity = totalHist[d-1] / lengthPerDeg;
		for(int i=0; i < lengthPerDeg; i++){
			regulatHist[regularIndx + i] = multiplicity;
		}
		int multResidue = totalHist[d-1] % lengthPerDeg;
		for(int i=0; i<multResidue; i++){
			regulatHist[regularIndx + i] += 1;
		}
		regularIndx += lengthPerDeg;
	}

	if(regularIndx > degLenInfo->regularHistLength){
		return MISMATCHING_ARAAYS_LENGTHS;
	}

	return SUCCESS;
}

void distributeEvenlyMultiplicity(int* regularArray, int length, int totalMultiplicity){
	int multiplicity = totalMultiplicity / length;
	for(int i=0; i < length; i++){
		regularArray[i] = multiplicity;
	}
	int multResidue = totalMultiplicity % length;
	for(int i=0; i<multResidue; i++){
		regularArray[i] += 1;
	}
}

int getIndexOfDegRegularHist(DegLengthINFO* degLen, int d){
	if(d<0 || d >= degLen->totalHistLength){
		cout << "illegal degree-index";
		return -1;
	}
	int indxD = 0;
	for(int i=0; i<d; ++i){
		indxD += degLen->numIrrPolyPerDeg[i];
	}
	return indxD;
}


int getVectorWeight(int* vector, int length){
	if(vector == NULL) return -1;
	int weight = 0;
	for(int i=0; i<length; i++){
		weight += vector[i];
	}
	return weight;
}
void getSubvector(int* vector, int* subvector, int length, int indx, int weight);
unsigned long int getNumSubvectorsOfGivenWeight(int* vector, int length, int weight){
	int indx = 0;
	subvectorCounter = 0;
	int* subvector = NULL;
	subvector = new int[length];
	if(subvector==NULL){
		return -1;
	}
	for(int i=0; i<length; i++) subvector[i]=0;
	getSubvector(vector, subvector, length, indx, weight);

	delete[] subvector;
	return subvectorCounter;
}

void getSubvector(int* vector, int* subvector, int length, int indx, int weight){
	if(indx >= length){
		if(getVectorWeight(subvector, length) == weight){
			subvectorCounter++;
		}
		return;
	}

	int subvector_indxOriginal = subvector[indx];
	do{
		getSubvector(vector, subvector, length, indx+1, weight);
		++subvector[indx];
	}while(subvector[indx] <= vector[indx] &&
			getVectorWeight(subvector,indx+1)<= weight);
	subvector[indx] = subvector_indxOriginal;
}



int getInnerProduct(int* a, int* b, int length){
	if(a==NULL || b==NULL){
		return -1;
	}
	int sum = 0;
	for(int i=0; i<length; i++){
		sum += a[i]*b[i];
	}
	return sum;
}

int getInnerProductWithDifference(int* a1, int* a2, int* b, int length){
	if(a1==NULL || a2==NULL || b==NULL){
		return -1;
	}
	int sum = 0;
	for(int i=0; i<length; i++){
		sum += (a1[i]-a2[i])*b[i];
	}
	return sum;

}

void replaceArrayContent(int* array, int* newContent, int length){
	for(int i=0; i<length; i++){
		array[i] = newContent[i];
	}
}

// TODO ASSUMPTION ptr totalHist != NULL. make sure
int degreeTotalHist(int* totalHist, int lengthTotalHist){
	int degTotal = 0;
	for(int i=0; i< lengthTotalHist; i++){
		degTotal += (i+1)*totalHist[i];
	}
	return degTotal;
}

int getIndexOfRightmostNonZeroElement(int* array, int length){
	for(int i=0; i<length; i++){
		if(array[i]==0)
			return i-1;
	}
	return length-1; // case: no zeros in array. returns the index of last element
}

double averageValueArray(int *array, int length){
	int sum = 0;
	for(int i=0; i< length; i++){
		sum+=array[i];
	}
	return (double)sum/(double)length;

}

int maxWith0(int m){
	if(m>0)
		return m;
	else
		return 0;
}
void printArray(int* arr, int len){
	for(int i=0; i<len; i++)
		cout << arr[i] << " ";

	cout << endl;
}

void printHist(int* totalHist, int lengthTotal, int* regulatHist,
		int lengthRegular){
	cout << "Total histogram:" << endl;
	printArray(totalHist, lengthTotal);

	cout << "Regular histogram:" << endl;
	printArray(regulatHist, lengthRegular);
}

void printDegLenStruct(DegLengthINFO* degLen){
	if(degLen == NULL){
		cout << "null pointer"<< endl;
	}else{
		cout <<"n = "<<degLen->n << " , total length = " << degLen->totalHistLength <<
				" , regular length = " << degLen->regularHistLength << endl;
		printHist(degLen->numIrrPolyPerDeg, degLen->totalHistLength,
				degLen->degVecForRegularHist, degLen->regularHistLength);
	}
}
