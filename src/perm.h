#ifndef __PERM_H__
#define __PERM_H__

#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "random.h"
using namespace std;     

const double almostfull = 0.95;              // for hash table
const int scratch = 100;                      // uint32[scratch][n]

class Perm{
public:
	Perm();
	~Perm();

	void init(int n, int members);  
	int get_n();
	int get_members();
	uint32 hash(uint32 in, uint32 max);
	void display(uint32* A);
	void displayall();
	void multiply(uint32* A, uint32* B, uint32* result);
	uint32* newperm(uint32 k);
	void representbase(uint32 k, int result[]);               // always represent in base "members"
			// Example: members == 3, k == 531, then the rep. of 531 base 3 = 201200
			// result[0] = 2
			// ...
			// result[5] = 0
			// result[6] = -1 as terminating value
			// note: user is responsible for allocating enough storage to "result"
			// this leads to the permutation "2" * "0" * "1" * "2" * "0" * "0"

private:
	int n;												// n = n-permutation (i.e. perm on n elements)
	int members;											// number of distinct random perms we use as base for future perms.
	int representation[scratch];										// representation of <number> base <k>
	uint32** perms;
	uint32** interm;											// temporary storge when multiplying permutations together
	//Random r;

};
extern Perm per;
#endif
