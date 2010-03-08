#include "perm.h"

static int initperm = 0;

//************************************

Perm::Perm()
{
	perms = NULL;
	interm = NULL;
}

//************************************

Perm::~Perm()
{
	int i;
	if(perms != NULL){
		for(i = 0; i < members; i++){
			if(perms[i] != NULL){
				delete [] perms[i];
			}
		}
		delete [] perms;
	}
	if(interm != NULL){
		for(i = 0; i < scratch; i++){
			if(interm[i] != NULL){
				delete [] interm[i];
			}
		}
		delete [] interm;
	}
}

//************************************

void Perm::init(int a, int b)
{
  // int a: permutation on "a" objects
  // int b: number of random "base" permutation we will set up
	// a) assign memmory to "uint32** perms" -> perms[members][N]
	// b) populate perms[members] with random valid permutations used as base for future permutations
	int i,k,t;
	int* tally;
	uint32 rand;
	int filled;
	int done;

	if(initperm == 0){
		this->n = a;
		this->members = b;
		tally = new int[n];
//		representation = new int[scratch];
		interm = new uint32*[scratch];
		for(i = 0; i < scratch; i++){
			interm[i] = new uint32[n];
		}
		perms = new uint32*[members];
		for(i = 0; i < members; i++){
			perms[i] = new uint32[n];
			for(t = 0; t < 2; t++){
				for(k = 0; k < n; k++){
					tally[k] = 0;
				}
				filled = 0;
				while(filled < n){
					rand = (uint32)r.ranint(n);
					if(tally[rand] == 0){             // cell still unoccupied
						tally[rand]++;	
					}
					else if((tally[rand] == 1) && (almostfull > (double)(filled)/(double)(n))){
						done = 0;
						while(done == 0){
							rand = (uint32)r.ranint(n);
							if(tally[rand] == 0){
								done = 1;
								tally[rand]++;
							}
							else{
								rand = hash(rand, n);
								if(tally[rand] == 0){
									done = 1;
									tally[rand]++;
								} // end if(tally[rand] == 0)
							} // end if(tally[rand] == 0)
						} // end while(done == 0)

					}
					else{							// the "tally[]" hashtable is "almostfull" (like >= 95%), now fill it up all the way
						rand = (uint32)r.ranint(n);
						while(tally[rand] == 1){
							rand++;
							rand %= n;
						}
						tally[rand]++;
					} // end if(tally[rand] == 0)
					//perms[i][filled] = rand;
					interm[t][filled] = rand;
					filled++;
				} // end while(filled < n)
			} // end for(t = 0; t <2; t++)
			multiply(interm[0], interm[1], perms[i]);
		} // end for(i = 0; i < members; i++)
		delete [] tally;
		initperm = 1;
	}
//	else{
//		cout <<"Perm::init. Permutations have already been initialized previously."<<endl;
//	}
}

//************************************

int Perm::get_members()
{
	return(members);
}

//************************************

int Perm::get_n()
{
	return(n);
}

//************************************

uint32 Perm::hash(uint32 in, uint32 max)
{
	uint32 result;
	uint32 temp;
	uint32 sum = 0x00000000;
	int i;
	
	for(i = 0; i < 5; i++){
		temp = in >> i;
		sum += temp;
	}
	result = sum ^ in;
	result = result % max;
	return(result);
}

//************************************

void Perm::display(uint32* A)
{
	int i;
	for(i = 0; i < n; i++){
		cout <<i+1<<" ";
	}
	cout <<endl;
	for(i = 0; i < n; i++){
		cout <<A[i]+1<<" ";
	}
	cout <<endl;
}

//************************************

void Perm::displayall()
{
	int i,j;
	cout << "Displaying all base permutations, N = "<<n<<":"<<endl;
	for(i = 0; i < members; i++){
		cout <<"Base permutation "<<i<<":"<<endl; 
		display(perms[i]);
	}
	uint32* res = new uint32[n];
	for(i = 0; i < members; i++){
		for(j = 0; j < members; j++){
			multiply(perms[i], perms[j], res);
			cout <<"#########"<<endl;
			display(perms[i]);
			cout <<" times "<<endl;
			display(perms[j]);
			cout <<"equals"<<endl;
			display(res);
		}
	}
	delete [] res;
}

//************************************

void Perm::multiply(uint32* A, uint32* B, uint32* result)
{
	int i;
	uint32* temp = B;
	for(i = 0; i < n; i++, temp++){
		result[i] = A[*temp];
	}
}

//************************************

void Perm::representbase(uint32 k, int result[])
{
	uint32 divisor = (uint32)(this->members);
	int temp[20];

	int i = 0;
	if(k == 0){
		temp[i] = 0;
		result[i] = 0;
		i++;
	}
	while(k > 0){
		result[i] = k % divisor;
		temp[i] = k % divisor;
		k /= divisor;
		i++;
	}
	temp[i] = -1;
	result[i] = -1;
}

//************************************

uint32* Perm::newperm(uint32 k)
{
	int i;
	uint32* result = new uint32[n];
	representbase(k, representation);									// now representation holds k base "members", [0] = most sig. "digit", "-1" to terminate
	if((representation[0] < 0) || (representation[0] >= members)){
		cout <<"Perm::newperm. representation[0] = "<<representation[0]<<" ...bailing out."<<endl;
		exit(0);
	}
	if(representation[1] == -1){
		for(i = 0; i < n; i++){
			result[i] = perms[representation[0]][i];
		}
//		cout <<"in Perm::newperm, checking base permutations"<<endl;
//		for(i = 0; i < this->members; i++){
//			display(perms[i]);
//			cout <<"###"<<endl;
//		}
	}
	else if(representation[2] == -1){
		multiply(perms[representation[0]], perms[representation[1]], result);
	}
	else if(representation[3] == -1){
		multiply(perms[representation[0]], perms[representation[1]], interm[0]);
		multiply(interm[0], perms[representation[2]], result);
	}
	else{
		multiply(perms[representation[0]], perms[representation[1]], interm[0]);
//		cout <<"interm[0] = "<<endl;
//		display(interm[0]);
		i = 2;
		while(representation[i+1] != -1){ // will execute ONE or more times! so iterm[] are valid permutations
			multiply(interm[i-2], perms[representation[i]], interm[i-1]);
			i++;
//			cout <<"interm["<<i-1<<"] = "<<endl;
//			display(interm[i-1]);
		}
		multiply(interm[i-2], perms[representation[i]], result);
//		result = interm[i-2];
	}
//	cout <<"Perm::newperm. result = "<<endl;
//	display(result);
// start debug ////
//int dbg[1200];
//for(i=0;i<n;i++){
//	dbg[i] = result[i];
//}
	/// end debug /////
	return(result);
}

//************************************

