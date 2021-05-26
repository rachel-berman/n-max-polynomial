# n-max-polynomial

**Background

This code computes a combinatorial object I call an "n-maximal polynomial":
Given a finite field F=GF(q) and a positive integer n, let S(x) be a monic polynomial over F of degree at most 2n. 
An "n-factorization" of S(x) is an ordered pair (u(x), v(x)) of polynomials over F such that:
  1. S(x) = u(x) * v(x).
  2. u(x) and v(x) have degree at most n.
We say that S(x) is "n-maximal" if it has as many as possible n-factorizations (u(x), v(x)).
The research objectives: 
  1. Compute or give upper and lower bounds on the number of n-factorizations of an n-maximal polynomial, in terms of n and q. 
  2. Characterize the n-maximal polynomials in F.

This work is motivated by an application in the field of Coding Theory, specifically, the list decoding of a certain type of Rank-Metric codes. More details can be found in:
*R.M. Roth, ``On decoding rank-metric codes over large fields``, IEEE Trans. Inf. Theory 64 (2018)
944–951 (Section 4). 

The results of this research can be found in:
*Berman R.N. & Roth R.M., ``On the Number of Factorizations of Polynomials over
Finite Fields``, Journal of Combinatorial Theory, Series A. Volume 182 (2021, August), 105462,
https://doi.org/10.1016/j.jcta.2021.105462

**About the code:

An n-maximal polynomial is called `worstSx`.

The module `auxFunctions` includes auxiliary functions; some of them stem from the mathematical model in which I represent a polynomial.

The main module is `worstSx`. The main functions in this module are:
  1. void getWorstSx(int n, int q, int minDeg, int maxDeg);
  2. void getWorstSxAdvanced(int n, int q, bool printPossibleWorstPoly);
Given n and q, these functions compute an n-maximal polynomial and its number of n-factorizations. The algorithm in the first function is backtracking, with optimizations (pruning branches) due to mathematical properties of n-maximal polynomials that I have proved. The second function is more efficient compared to the first function. The algorithm in this function is based on a theorem that I have proved, which generalized the mathematical properties mentioned above.

Module `verifyBellCurveClaim` calculates certain quantities where the goal was to test the feasibility of a mathematical hypothesis developed during my research. The main function in this module is:
void printSizeBkCurve(int n, int q, int *SxHist, int SxHistLength, int Lc);
