#include <dune/MIOCP/Dmax.hh>

#include<iostream>
#include<vector>
#include<math.h>

using namespace std;

extern "C" {
  double ddot_(const int*,const double*,const int*,const double*,const int*);
}

int Dmax::separate(int n, const double*u, double*cut, double&b) const
{
	vector<int> index; // indices of local maxima and minima; start with maximum
	double sum=0.0; 
	double exval=0.0; // value of last extrema

	// correct value if lower bound is exceeded
	// Checking whether the first point is a local maxima or not 
	if (u[0] > u[1]){
		index.push_back(0); // update indices of local maxima 
		exval=u[0]; // update value of last extrema
	}

	// Iteration over all elements to check local minima and maxima 
	for (int i=1; i <n-1; i++) 
	{
		if( (u[i-1]>= u[i]) && (u[i]< u[i+1]) &&  u[i]<exval){ // u[i] local minimum
			index.push_back(i); // update indices with local minimum
			exval=u[i]; // update value of last extrema
		}
		else if ( (u[i-1]<= u[i]) && (u[i]> u[i+1]) && u[i]>exval){ // u[i] local maximum
			index.push_back(i); // update indices with local maxima 
			exval=u[i]; // update value of last extrema
		}
	}

	// Last point is a local maximum or minimum
	index.push_back(n-1);
	
	// clear vector to save cutting plane 
	for(int i=0; i < n; i++)
		cut[i] = 0;

	// case (number of extrema - sigma) even 
	if ((sigma & 1)==(index.size() & 1)){  
		for(int i=0; i < index.size()-1; i++) 
			cut[index[i]]=(i&1)? -1: 1; 
	}
	// cases (number of extrema - sigma) odd
	else{
		for(int i=0; i < index.size(); i++) 
			cut[index[i]]=(i&1)? -1: 1;
	}
		
	// prove feasibility 
	int inc=1;
	sum = ddot_(&n,u,&inc,cut,&inc);
	b=floor(sigma*0.5); 
	if(sum <= b*1.01){
		cout << "\n control is feasible" << endl;
		return 0; 
	}
	else{
		cout << "\n control is infeasible" << endl;
		cout << "violation: " << sum-b << endl;
		return 1;
	}	
}

int Dmax::optimize(int n,const double* c, double* u) const	
{ 
	double M=1.0; 
	for(int i=0; i < n; i++) 
		M+=abs(c[i]); 
			
	// dynamic programming over timesteps i, values b in {0,1} and number l of switchings
	// P[i][b][l] encodes \sum_{k=0}{i}c_i u^*_i of best solution u^* with 
	// - value b 
	// - at time point i and 
	// - at most l switchings before time point i
	double*** P = new double**[n];

	// S[i][b][l]=1, 1<=i<=n, if best solution u^* for sub problem P[i][b][l] has
	// - switched from value 1-b at time point i-1 
	// - to b at time point i 
	// - and has at most l switchings
	double*** S = new double**[n]; 

	for(int i=0; i < n; i++){
		P[i]=new double*[2];
		S[i]=new double*[2]; 
	}
	for (int i=0; i < n; i++){
		for(int b=0; b < 2; b++){
			P[i][b]=new double[sigma+1];
			S[i][b]=new double[sigma+1];
		}
	}

	// initial values
	for(int l=0; l <= sigma; l++){
		
		// P[0][0][l]=0 for 0<=l<=sigma, S[0][0][l]=1 
		P[0][0][l]=0; 
		S[0][0][l]=0;
			
		//  P[0][1][0]=M and P[0][1][l]=c_0, S[0][1][] for 1<=l<=sigma
		if(l==0){
			P[0][1][l]= M; 
			S[0][1][l]=0;
		}
		else{
			P[0][1][l]= c[0];  
			S[0][1][l]=1;
		}
	}
		
	// recursion formula: 
	// P[i][0][l]=min( P[i-1][0][l], P[i-1][1][l-1] )	
	// P[i][1][l]=min( P[i-1][1][l], P[i-1][0][l-1] ) + c_i	
	for(int i=1; i < n; i++){
		for(int l=0; l <= sigma; l++){
			if(l==0){ 
				P[i][0][l]=P[i-1][0][l]; 
				S[i][0][l]=0;
				P[i][1][l]=P[i-1][1][l]+c[i];
				S[i][1][l]=0;
			}
			else{
				if(P[i-1][0][l]<= P[i-1][1][l-1]){
					P[i][0][l]=P[i-1][0][l]; 
					S[i][0][l]=0;
				}
				else{
					P[i][0][l]=P[i-1][1][l-1];
					S[i][0][l]=1;
				}
				if(P[i-1][0][l-1]<=P[i-1][1][l]){
					P[i][1][l]=P[i-1][0][l-1]+c[i]; 
					S[i][1][l]=1; 
				}
				else{
					P[i][1][l]=P[i-1][1][l]+c[i];
					S[i][1][l]=0;
				}
			}
		}
	}
	// backtracking of best solution
	int k=sigma;
	int b; // value at current time point
	// start: min( P[n-1][1][sigma], P[n-1][0][sigma]) 
	if(P[n-1][0][sigma]<=P[n-1][1][sigma])
		b=0;
	else 
		b=1;
	u[n-1]=b;
	// recursion
	for(int i=n-1; i>0; i--){
		if(S[i][b][k]==1){ // shifted from value 1-b at time point i-1 to b at time point i
			u[i-1]=1-b; 
			b=1-b; 
			k-=1; // one swichting less left
		}
		else
			u[i-1]=b;	
	}
	// free memory	
	for(int i=0; i<n; i++){
		for(b=0; b < 2; b++){
			delete[] P[i][b];
			delete[] S[i][b];
		}
		delete[] P[i];
		delete[] S[i];
	}
	delete[] P; 
	delete[] S;

	return 0;
} 
