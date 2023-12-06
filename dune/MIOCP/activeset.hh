// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later

#ifndef ACTIVE_SET_HH
#define ACTIVE_SET_HH

#include "heatdriver.hh" 
#include "adjointdriver.hh"
#include "linsolver.hh"
#include <dune/istl/matrix.hh>

using namespace Dune; 
using namespace std; 

/* Structure for active set method result
 */
struct ActiveSetResult
{
	int status; 
	/* status overview:
	 * -2: max iteration count reached
	 * -1: time limit reached
	 *  0: problem solved
	 */
	double defect; // active set defect
	int iter_activeset; // active set iterations
	int iter_lin; // linear solver iterations
	double time_solve; // linear solver time in active set method 	
	double time_assemble; // linear assemble time in active set method 

	ActiveSetResult(){
		status=0;
		iter_activeset=0;
		iter_lin=0; 
		time_solve=0;
		time_assemble=0;	
	} 
};


/* Primal dual active set method to solve optimal control problems of the form 
 * 
 *    min J(y,u)=1/2|y-y_d|_{L^2(Q)}^2 + \apha/2 |u-u_d|_{L^2(0,T)}^2
 *    s.t.  \partial_t y - \Delta y = Gu + f = \sum_j u_j(t)*Psi_j(x) + f	in Q=\Omega x (0,T)
 *                       		   y = g   				on \Sigma_D=\Gamma_D x (0,T)
 *         		    \nabla y \cdot n = j   				on \Sigma_N=\Gamma_N x (0,T)
 *                        		y(0) = y0 		
 * 					u_a <= u <= u_b 			a.e. in (0,T)
 *
 * [with switching constraints Au\leq b]
 * 
 * 
 *
 */ 
template<class PStruct> 
class ActiveSet
{
	public:
	typedef typename PStruct::GFS GFS;
	typedef typename PStruct::GV GV;
	typedef typename PStruct::RangeType RangeType; // range type discrete grid function
	typedef typename PStruct::vtype vtype; // type of Dof vector 

	typedef typename PStruct::ControlVector ControlVector; // type of control vector
	typedef typename ControlVector::field_type field_type; // field type of control vector
	typedef typename PStruct::Matrix Matrix; // matrix type

	typedef typename PStruct::AdjProblem AdjProblem;
	typedef typename PStruct::CGProblem CGProblem;

	using Vector = BlockVector<Dune::FieldVector<double,1>>; 
	using Result= ActiveSetResult; 

	enum{n=PStruct::n}; // number of switches
	
	private: 
	const GFS& gfs; //  spatial grid function space
	const vector<double>& dt; // temporal grid
	HeatProblemInterface<GV,1>& heat; // heat problem (with boundary data)
	const ControlVector& Gstar_eta; //  encodes G*S*(Sf-yd) for residual 

	int iter_activeset;
	int iter_lin;
	double reduction; 
	double alpha;
	double rho;  
	double tol;

	public: 
	ActiveSet(const GFS& gfs_, 
	const vector<double>& dt_, 
	HeatProblemInterface<GV,1>& heat_,
	const ControlVector& Gstar_eta_, 
	int iter_active, 
	int iter, 
	double factor, 
	double alpha_, 
	double rho_, 
	double tol_) : gfs(gfs_), dt(dt_), heat(heat_),  Gstar_eta(Gstar_eta_), iter_activeset(iter_active), iter_lin(iter), reduction(factor), alpha(alpha_), tol(tol_), rho(rho_)
	{}
	
	// active set method (without switching constraints)
	void apply(ControlVector& u, vector<vtype>& y, ControlVector& Gstar_zeta, Result &result) const;

	// active set method (with switching constraints)
	void apply(ControlVector& u, vector<vtype>& y,  ControlVector& Gstar_zeta,Vector& mu,const Matrix& A,const Vector& b, Result &result,bool reopt, double timelimit=86400) const;

	// update active switching constraints
	void update(const Matrix& A,const Vector& b,const ControlVector& u, const Vector& mu,vector<int>& B) const; 

	// update active box constraint (if switching constraints existent)
	void update(const Matrix& A,const ControlVector& Gstar_zeta, const Vector& mu, vector<pair<int,int>>& A_minus, vector<pair<int,int>>& A_plus, vector<pair<int,int>>& I) const;
	// update active box constraint (if switching constraints not existent)
	void update(const ControlVector& Gstar_zeta, vector<pair<int,int>>& A_minus, vector<pair<int,int>>& A_plus, vector<pair<int,int>>& I) const; 

	// cg method (if switching constraints not existent)
	void solve(const vector<pair<int,int>>& A_minus, const vector<pair<int,int>>& A_plus, 
	const vector<pair<int,int>>& I, ControlVector& u,Result &result) const;

	// minimum residual solver (if swixchting constraints existent)
	void solve(const vector<pair<int,int>>& A_minus, const vector<pair<int,int>>& A_plus, 
	const vector<pair<int,int>>& I,const vector<int>& B, ControlVector& u,Vector& mu, const Matrix& A,const Vector& b, Result &result) const;

};



/* Primal dual active set method (without swichting constraints) 
*/
template<class PStruct> 
void ActiveSet<PStruct>::apply(
ControlVector& u, // start control
vector<vtype>& y, // encodes SGu
ControlVector& Gstar_zeta, // encodes G*S*SGu 
Result & result ) const //storage of results
{	


	cout << "\n === Active Set Algorithm \n" << endl;
	
	bool flag = false; 
	int step = 0; 

	// previous active box constraints
	vector<pair<int,int>> A_minus_old; 
	vector<pair<int,int>> A_plus_old; 
	

	while (step <= iter_activeset && !flag) 
	{
	// active and inactive box constraints (at a time and point) 
	vector<pair<int,int>> A_minus; 
	vector<pair<int,int>> A_plus; 
	vector<pair<int,int>> I; 
	
	
	if(step > 0){
		// update active box constraints
		update(Gstar_zeta,A_minus,A_plus,I);
	}
	else{
		for(int i=0; i < dt.size(); i++)
			for(int j=0; j < n; j++)
				I.push_back(make_pair(i,j)); 
	}
	
	if ( (A_minus_old == A_minus) && (A_plus == A_plus_old) && step >0){
		flag=true;
		continue;
	}

	// calculate residual of old solution under new active sets
	double ini_defect=0; 
	// A- part
	for (auto const& p : A_minus){
		double time=accumulate(dt.begin(), dt.begin()+p.first+1,0.0);
		double val=alpha*(u[p.first][p.second]-heat.ua(time)[p.second]);
		ini_defect+=dt[p.first]*pow(val,2);
	}
	// A+ part
	for (auto const& p : A_plus){
		double time=accumulate(dt.begin(), dt.begin()+p.first+1,0.0);
		double val=alpha*(u[p.first][p.second]-heat.ub(time)[p.second]);
		ini_defect+=dt[p.first]*pow(val,2);
	}
	// I part
	for (int i=0; i < I.size(); i++){
		double time=accumulate(dt.begin(), dt.begin()+I[i].first+1,0.0);
		double val=-Gstar_eta[I[i].first][I[i].second]-Gstar_zeta[I[i].first][I[i].second]-alpha*(u[I[i].first][I[i].second]-heat.ud(time)[I[i].second]); 
		val*=dt[I[i].first];
		ini_defect+=pow(val,2);
	}

	if (sqrt(ini_defect) < tol && step>0){
		flag=true;
	}
	else 
	{
		
		cout << "\nActive Set Iteration " << step << "\n" << endl; 
#ifdef NEWTONDEBUG
		cout << "Indices in A_minus: "<< endl; 
		for(auto const& p: A_minus) 
			cout << "("<< p.first << "," << p.second << ")" << endl;

		cout << "Indices in A_plus: "<< endl; 
		for(auto const& p: A_plus) 
			cout << "("<< p.first << "," << p.second << ")" << endl;
	
		cout << "Indices in I: "<< endl; 
		for(auto const& p: I) 
			cout << "("<< p.first << "," << p.second << ")" << endl;
#endif		
		
		// CG method 
		solve(A_minus,A_plus,I,u,result);

		
		// accept step
#ifdef NEWTONDEBUG
		cout << "new control: " << endl;
		for(int i=0; i < u.size(); i++) 
			cout << u[i] << endl;
#endif
		A_minus_old.clear();
		A_minus_old = A_minus; 
		A_plus_old.clear();
		A_plus_old = A_plus; 
 		
		// prepare box constraints update
		CGProblem cgproblem(u,dt); 	
		heatdriver(gfs,cgproblem,dt,y);
		AdjProblem adjProblem(dt,y,gfs);
		vector<vtype> zeta;
		adjointdriver(gfs,adjProblem,zeta);
		Gstar(gfs,heat,zeta,Gstar_zeta);

		step++;
	}
	}
	cout << "\n=== End Active Set Algorithm " << endl;
	if(!flag){
		result.status=-2;
		cout << "Max iteration count of Active Set method reached" << endl;
	}
	
	// add number of active set iterations
	result.iter_activeset+=step; 
	
}



/* Primal dual active set method (with switching constraints)
 * update of active sets all at once
 */
template<class PStruct> 
void ActiveSet<PStruct>::apply(
ControlVector& u, // start control
vector<vtype>& y, // encodes SGu
ControlVector& Gstar_zeta, // encodes G*S*SGu 
Vector& mu, // start lagrange multiplicator of switching constraints 
const Matrix& A, // matrix of switching constraints
const Vector& b, // right hand side of switching constraints
Result &result, //storage of results
bool reopt, // reoptimization flag
double timelimit) const // rest running time
{	
	cout << "\n ==== Active Set Algorithm \n" << endl;
	

	bool flag = false; 
	int step = 0; 

	// previous active box constraints and switching constraints
	vector<pair<int,int>> A_minus_old; 
	vector<pair<int,int>> A_plus_old; 
	vector<int> B_old; 

	double time=0.0;
	while (step < iter_activeset && !flag) 
	{
	clock_t cstart=clock(); 
	// active and inactive box constraints (at a time and point) 
	vector<pair<int,int>> A_minus; 
	vector<pair<int,int>> A_plus; 
	vector<pair<int,int>> I; 
	vector<int> B;

	if((step==0) && (A.N()==1 || !reopt)){
		for(int i=0; i < dt.size(); i++)
			for(int j=0; j < n; j++)
				I.push_back(make_pair(i,j)); 
		B.push_back(A.N()-1);
	}
	else{
		// update of active switching constraints
		update(A,b,u,mu,B); 
		// update of active box constraints
		update(A,Gstar_zeta,mu,A_minus,A_plus,I); 
	}
	
	if ( (A_minus_old == A_minus) && (A_plus == A_plus_old) && (B_old == B) && step>0){ // active sets unchanged
		flag=true;
		continue;
	}
	
	// calculate residual of old solution under new active sets
	double ini_defect=0; 
	// A- part
	for (auto const& p : A_minus){
		double time=accumulate(dt.begin(), dt.begin()+p.first+1,0.0);
		double val=alpha*(u[p.first][p.second]-heat.ua(time)[p.second]);
		ini_defect+=dt[p.first]*pow(val,2);
	}
	// A+ part
	for (auto const& p : A_plus){
		double time=accumulate(dt.begin(), dt.begin()+p.first+1,0.0);
		double val=alpha*(u[p.first][p.second]-heat.ub(time)[p.second]);
		ini_defect+=dt[p.first]*pow(val,2);
	}
	// I part
	Vector h(A.M()); // auxiliary variable
	A.mtv(mu,h); // h=A^Tmu
	for (int i=0; i < I.size(); i++){
		double time=accumulate(dt.begin(), dt.begin()+I[i].first+1,0.0);
		double val=-Gstar_eta[I[i].first][I[i].second]-Gstar_zeta[I[i].first][I[i].second]-alpha*(u[I[i].first][I[i].second]-heat.ud(time)[I[i].second]); 
		val*=dt[I[i].first];
		val-=h[I[i].first*n+I[i].second];
		ini_defect+=pow(val,2);
	}
	for(int i=0; i < dt.size(); i++)
		for(int j=0; j < n; j++)
				h[i*n+j]=u[i][j];

	Vector z(A.N()); // auxiliary variable
	A.mv(h,z); // calculate Au
	for(int i=0; i < B.size(); i++)
		ini_defect+=pow(b[B[i]]-z[B[i]],2);
	// N part
	for(int i=0; i< A.N(); i++){
		if(find(B.begin(),B.end(),i)==B.end())
			ini_defect+=pow(rho*mu[i],2); 
	}	
	if(sqrt(ini_defect)<tol && step>0) // initial defect small enough
		flag=true;
	else 
	{
		if(time>=timelimit) // time limit reached -> stop active set 
			break;

		cout << "\n Active Set Iteration " << step   << "\n" << endl; 


#ifdef NEWTONDEBUG
		cout << "Indices in A_minus: "<< endl; 
		for(auto const& p: A_minus) 
			cout << "("<< p.first << "," << p.second << ")" << endl;

		cout << "Indices in A_plus: "<< endl; 
		for(auto const& p: A_plus) 
			cout << "("<< p.first << "," << p.second << ")" << endl;
	
		cout << "Indices in I: "<< endl; 
		for(auto const& p: I) 
			cout << "("<< p.first << "," << p.second << ")" << endl;

		cout << "Active switching constraints: "<< endl; 
		for(auto& i: B) 
			cout << i << endl;
#endif		
				

		if (B.size() > 0 && I.size()>0){ // active switching constraints

			// matrix restricted to the active indices in B			
			Matrix A_B;
			A_B.setBuildMode(Matrix::row_wise);
			A_B.setSize(B.size(),A.M());
			for(auto row=A_B.createbegin(); row!=A_B.createend(); ++row){
				for( auto it=A[B[row.index()]].begin(); it!=A[B[row.index()]].end(); ++it) 
					row.insert(it.index()); 
			}
			for(int i=0; i < B.size(); i++) 
				A_B[i]=A[B[i]]; 

			// MINRESSolver
			solve(A_minus,A_plus,I,B,u,mu,A_B,b,result);
			

		}
		else if(B.size()==0 && I.size()>0){ // no active switching constraints
			
			// CG method 
			solve(A_minus,A_plus,I,u,result);

			// new lagrange multiplicator
			mu=0.0;
		}
		else{
			// auxiliary variable to save boundary values of active indices
			ControlVector h(dt.size());
			for (auto const& p : A_minus){
				double time=accumulate(dt.begin(),dt.begin()+p.first+1,0.0);
				h[p.first][p.second]=heat.ua(time)[p.second]; 
			}
    			for (auto const& p : A_plus){
				double time=accumulate(dt.begin(),dt.begin()+p.first+1,0.0);
				h[p.first][p.second]=heat.ub(time)[p.second]; 
			}

			// new control 
			u=h;

			// new lagrange multiplicator
			mu=0.0;
		}
		
	
		// accept step
#ifdef NEWTONDEBUG
		cout << "new control: " << endl;
		for(int i=0; i< u.size(); i++) 
			cout << u[i] << endl; 
	
		cout << "new lagrange multiplicator: " << endl;
		for(int i=0; i< mu.size(); i++) 
			cout << mu[i] << endl; 
#endif

		A_minus_old.clear();
		A_minus_old = A_minus; 
		A_plus_old.clear();
		A_plus_old = A_plus; 
		B_old.clear();
		B_old = B;

		// prepare box constraints update
		CGProblem cgproblem(u,dt); 	
		heatdriver(gfs,cgproblem,dt,y);
		AdjProblem adjProblem(dt,y,gfs);
		vector<vtype> zeta;
		adjointdriver(gfs,adjProblem,zeta);
		Gstar(gfs,heat,zeta,Gstar_zeta);

		step++;
			
	}
	time+=(clock()-cstart)/CLOCKS_PER_SEC;
	}
	cout << "\n===" << "End Active Set Algorithm \n" << endl;
	if(!flag && step >= iter_activeset){
		result.status=-2;
		cout << "Max iteration count of Active Set method reached" << endl;
	}
	else if(!flag && time > timelimit) 
		result.status=-1;
	

	// add number of active set iterations
	result.iter_activeset+=step; 
}


/* update active switching constraints
 */
template<class PStruct> 
void ActiveSet<PStruct>::update(
const Matrix& A, // matrix of switching constraints
const Vector& b, // right hand side of switching constraints
const ControlVector& u, // current control u
const Vector& mu, // current lagrange multiplicator of switching constraints
vector<int>& B) const // storage for active switching constraints
{
	// clear old active set 
	B.clear();

	// calculate h=Au+rho*mu
	Vector h = mu;
	h*=rho;  
	// store u appropiately
	Vector x(A.M());
	for(int i=0; i<dt.size(); i++){
		for(int j=0; j < n; j++)
			x[i*n+j]=u[i][j];
	}
	// add Au to h
	A.umv(x,h);
 
	for(int i=0; i < h.size(); i++){
		if( h[i]>b[i]) // index in B
			B.push_back(i);
	}
}


/* update active box constraints (if switching constraints not existent)
 */
template<class PStruct> 
void ActiveSet<PStruct>::update(
const ControlVector &Gstar_zeta, // encodes G*S*SGu 
vector<pair<int,int>>& A_minus, // storage for active lower box constraints  A⁻
vector<pair<int,int>>& A_plus, // storage for active upper box constraints  A⁺
vector<pair<int,int>>& I) const // storage for inactive box constraints
{

	// clear old sets
	A_minus.clear();
	A_plus.clear();
	I.clear();
	
	// update of active box constraints
	for(int i=0; i < dt.size(); i++){
		double time=accumulate(dt.begin(), dt.begin()+i+1,0.0);
		FieldVector<double,n> ua= heat.ua(time);
		FieldVector<double,n> ud= heat.ud(time);
		FieldVector<double,n> ub= heat.ub(time); 
		for(int j=0; j < n; j++){
			double value=-Gstar_zeta[i][j]-Gstar_eta[i][j]+alpha*ud[j]; 
			value*=dt[i];
			if ( ( value-alpha*dt[i]*ua[j]) < 0.0) // index in A⁻
				A_minus.push_back(make_pair(i,j)); 
			else if ((  value-alpha*dt[i]*ub[j] )> 0.0) // index in A⁺
				A_plus.push_back(make_pair(i,j)); 
			else 
				I.push_back(make_pair(i,j)); 
		}
	}
}


/* update active box constraints (if switching constraints existent) 
 */
template<class PStruct> 
void ActiveSet<PStruct>::update(const Matrix& A, // matrix of switching constraints
const ControlVector &Gstar_zeta, // encodes S*SGu 
const Vector& mu,  // current lagrange multiplicator of switching constraints
vector<pair<int,int>>& A_minus, // storage for active lower box constraints  A⁻
vector<pair<int,int>>& A_plus, // storage for active upper box constraints  A⁺
vector<pair<int,int>>& I) const // storage for inactive box constraints
{

	// clear old sets
	A_minus.clear();
	A_plus.clear();
	I.clear();
	
	// update of box constraints
	Vector h(A.M());		
	A.mtv(mu,h); // h=A^Tmu
	for(int i=0; i < dt.size(); i++){
		double time=accumulate(dt.begin(), dt.begin()+i+1,0.0);
		FieldVector<double,n> ua= heat.ua(time);
		FieldVector<double,n> ud= heat.ud(time);
		FieldVector<double,n> ub= heat.ub(time); 
		for(int j=0; j < n; j++){
			double value=-Gstar_zeta[i][j]-Gstar_eta[i][j]+alpha*ud[j]; 
			value*=dt[i];
			value-=h[i*n+j];
			if ( ( value-alpha*dt[i]*ua[j] ) < 0.0) // index in A⁻
					A_minus.push_back(make_pair(i,j)); 
				else if (( value-alpha*dt[i]*ub[j] ) > 0.0 ) // index in A⁺
					A_plus.push_back(make_pair(i,j)); 
				else 
					I.push_back(make_pair(i,j)); 
		}
	}
}



/* cg method (if switching constraints not existent) 
 */
template<class PStruct> 
void ActiveSet<PStruct>::solve(
const vector<pair<int,int>>& A_minus, // active lower box constraints  A⁻
const vector<pair<int,int>>& A_plus, // active upper box constraints A⁺
const vector<pair<int,int>>&I, // inactive box constraints
ControlVector& u, // current control
Result &result) const // storage for results
{
	clock_t cstart=clock();	
	ControlVector h(dt.size());  // auxiliary variable to save boundary values of active indices
	bool null=true; // indicate whether boundary values are all zero
	for (auto const& p : A_minus){
		double time=accumulate(dt.begin(), dt.begin()+p.first+1,0.0);
		h[p.first][p.second]=heat.ua(time)[p.second];  
		if(heat.ua(time)[p.second]!=0)
			null=false;
	}
    	for (auto const& p : A_plus){
		double time=accumulate(dt.begin(),dt.begin()+p.first+1,0.0);
		h[p.first][p.second]=heat.ub(time)[p.second]; 
		if(heat.ub(time)[p.second]!=0)
			null=false;
	}

	// u restricted to inactive box constraints
	Vector u_I(I.size()); 
	for (int i=0;i<I.size(); i++ )
		u_I[i]=u[I[i].first][I[i].second]; 

	// setup right hand side
	Vector r(I.size());
	// add \Chi_I(alpha*ud+G*S*yd-G*S*Sf)-\Chi_I(G*S*SGua+G*S*SGub)
	{
		if(!null){
			cout << "\n Assemble G*S*SG(ua+ub)" << endl;
			// calculate (S*SGua+S*SGub)
			CGProblem cgproblem(h,dt); 
			vector<vtype> w; 
			heatdriver(gfs,cgproblem,dt,w);
			AdjProblem adjointProblem=AdjProblem(dt,w,gfs); 
    			vector<vtype> v; 
			adjointdriver(gfs,adjointProblem,v);
			// z= G*S*SG(ua+ub)
			ControlVector z(dt.size());
			Gstar(gfs,cgproblem,v,I,z); 

			for (int i=0;i<I.size(); i++ ){
				double time=accumulate(dt.begin(), dt.begin()+I[i].first+1,0.0);
				r[i]=alpha*heat.ud(time)[I[i].second]-Gstar_eta[I[i].first][I[i].second]-z[I[i].first][I[i].second];
				r[i]*=dt[I[i].first]; 
			}
		}
		else{
			for (int i=0;i<I.size(); i++ ){
				double time=accumulate(dt.begin(), dt.begin()+I[i].first+1,0.0);
				r[i]=alpha*heat.ud(time)[I[i].second]-Gstar_eta[I[i].first][I[i].second];
				r[i]*=dt[I[i].first]; 
			}
		}
	}
	result.time_assemble+=(double)(clock()-cstart)/CLOCKS_PER_SEC;

	// linear operator
	typedef LinOperator<PStruct> LinearOperator;
	LinearOperator op(gfs,dt,alpha,I);

	ScalarProduct<Vector> sp;

	// CG solver
	Dune::Richardson<Vector,Vector> prec(1.0); // no preconditioning
	OwnCGSolver<Vector> solver(op,sp,prec,reduction, iter_lin, 2,0);

	// storage for results
 	Dune::InverseOperatorResult res;
		
	// call CG solver
	cstart=clock();	
	solver.apply(u_I,r,tol,res);
	result.time_solve+=(double)(clock()-cstart)/CLOCKS_PER_SEC;
	
	// accept step
	// add boundary values on active box constraints
	u=h; 
	// add values on inactive box constraints
	for(int i=0; i< I.size(); i++) 
		u[I[i].first][I[i].second]=u_I[i]; 

	result.iter_lin+=res.iterations;
}



/* minressolver (if switching constraints existent) 
 */
template<class PStruct> 
void ActiveSet<PStruct>::solve(
const vector<pair<int,int>>& A_minus, // active lower box constraints  A⁻
const vector<pair<int,int>>& A_plus, // active upper box constraints A⁺
const vector<pair<int,int>>&I, // inactive box constraints
const vector<int>& B, // active switching constraints
ControlVector& u, // current control
Vector& mu, // current lagrange multiplicator of switching constraints
const Matrix& A_B, // matrix of active switching constraints
const Vector& b,  // right hand side of active switching constraints
Result &result) const // storage for results
{
	clock_t cstart=clock();	
	ControlVector h(dt.size());  // auxiliary variable to save boundary values of active indices
	bool null=true; // indicate whether boundary values are all zero
	for (auto const& p : A_minus){
		double time=accumulate(dt.begin(), dt.begin()+p.first+1,0.0);
		h[p.first][p.second]=heat.ua(time)[p.second];  
		if(heat.ua(time)[p.second]!=0)
			null=false;
	}
    	for (auto const& p : A_plus){
		double time=accumulate(dt.begin(),dt.begin()+p.first+1,0.0);
		h[p.first][p.second]=heat.ub(time)[p.second]; 
		if(heat.ub(time)[p.second]!=0)
			null=false;
	}
	// initialize right hand side
	Vector rhs(I.size()+B.size());
	// initialize start vector
	Vector x(I.size()+B.size());	

	// right hand side restricted to the active indices in B
	for(int i=0; i < B.size(); i++) 
		rhs[i+I.size()]=b[B[i]];
	// subtract A_Bh (h - extended vector of boundary values)
	{
		Vector z(B.size());
		// store h appropiately
		Vector nu(A_B.M());
		for(int i=0; i< dt.size(); i++){
			for(int j=0; j < n; j++)
				nu[i*n+j]=h[i][j];
		}
		// calculate A_Bh
		A_B.mv(nu,z);
		for(int i=0; i < B.size(); i++) 
			rhs[i+I.size()]-=z[i];
	}

	// lagrange multiplicator restricted to the active sets
	for(int i=0; i < B.size(); i++) 
		x[i+I.size()]=mu[B[i]];


	// setup right hand side
	// add \Chi_I(alpha*ud+G*S*yd-G*S*Sf)-\Chi_I(G*S*SGua+G*S*SGub)
	{
		if(!null){
			cout << "\n Assemble G*S*SG(ua+ub)" << endl;
			// calculate (S*SGua+S*SGub)
			CGProblem cgproblem(h,dt); 
			vector<vtype> w; 
			heatdriver(gfs,cgproblem,dt,w);
			AdjProblem adjointProblem=AdjProblem(dt,w,gfs); 
    			vector<vtype> v; 
			adjointdriver(gfs,adjointProblem,v); 
			// z= G*S*SG(ua+ub)
			ControlVector z(dt.size());
			Gstar(gfs,cgproblem,v,I,z);

			for (int i=0;i<I.size(); i++ ){
				double time=accumulate(dt.begin(), dt.begin()+I[i].first+1,0.0);
				rhs[i]=alpha*heat.ud(time)[I[i].second]-Gstar_eta[I[i].first][I[i].second]-z[I[i].first][I[i].second];
				rhs[i]*=dt[I[i].first]; 
			}
		}
		else{
			for (int i=0;i<I.size(); i++ ){
				double time=accumulate(dt.begin(),dt.begin()+I[i].first+1,0.0);
				rhs[i]=alpha*heat.ud(time)[I[i].second]-Gstar_eta[I[i].first][I[i].second];
				rhs[i]*=dt[I[i].first]; 
			}
		}
	}
	result.time_assemble+=(double)(clock()-cstart)/CLOCKS_PER_SEC;

			
	// u restricted to inactive box constraints
	for (int i=0;i<I.size(); i++ )
		x[i]=u[I[i].first][I[i].second]; 

	// linear operator
	typedef LinCutOperator<PStruct> LinearOperator;
	LinearOperator op(gfs,dt,alpha,I,A_B);

	ScalarProduct<Vector> sp;

	// MINRESSolver
	Dune::Matrix<double> AAT(B.size(),B.size()); 
	for(int i=0; i < B.size(); i++){
		for(int k=0; k < B.size(); k++){
			AAT[i][k]=0; 
			for(int l=0; l < I.size(); l++){
				int h=I[l].first*n+I[l].second;
				if(A_B.exists(i,h) && A_B.exists(k,h))
					AAT[i][k]+=A_B[i][h]*A_B[k][h];
			}
		}
	}

	Matrix P=Matrix(I.size()+B.size(),I.size()+B.size(),Matrix::row_wise);
	for(auto row=P.createbegin(); row!=P.createend(); ++row){
                if (row.index()< I.size())
			row.insert(row.index());
		else{
			for(int j=I.size(); j < P.N(); j++) 
				row.insert(j); 
		}
	}
	for(int i=0; i < I.size(); i++)
		P[i][i]=alpha; 
	for(int i=0; i < B.size(); i++)
		for(int j=0; j < B.size(); j++) 
			P[I.size()+i][I.size()+j]=1.0/alpha*AAT[i][j]; 
	InverseP<Matrix,Vector,Vector> invOperator(P);
	InverseOperator2Preconditioner<InverseP<Matrix,Vector,Vector>> prec(invOperator); // preconditioning
	OwnMINRESSolver<Vector> solver(op,sp,prec,reduction, iter_lin, 2);

	// storage for results
 	Dune::InverseOperatorResult res;
	
	// call MINRESSolver
	cstart=clock();	
	solver.apply(x,rhs,tol,res);
	result.time_solve+=(double)(clock()-cstart)/CLOCKS_PER_SEC;

	// accept step 
	// add boundary values on active box constraints
	u=h; 
	// add values on inactive box constraints
	for(int i=0; i< I.size(); i++) 
		u[I[i].first][I[i].second]=x[i];
	// new lagrange multiplicator
	mu=0.0;
	for(int i=0; i < B.size(); i++) 
		mu[B[i]]=x[i+I.size()]; 

	result.iter_lin+=res.iterations;
}




#endif //ifned ACTIVE_SET_HH
