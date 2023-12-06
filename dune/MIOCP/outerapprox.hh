// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef OUTER_APPROX_HH
#define OUTER_APPROX_HH

#include <dune/MIOCP/activeset.hh> 

using namespace Dune; 
using namespace std; 

template<class PStruct> 
class OuterApprox
{
public:
	// type of dof vector
    	typedef typename PStruct::vtype vtype;
	// discrete grid function
	typedef typename PStruct::DGF DGF;
	// range type of  discrete grid funtions
    	typedef typename PStruct::RangeType RangeType; 
	// grid type 
	typedef typename PStruct::GV GV; 
	
	// type of control vector
	typedef typename PStruct::ControlVector ControlVector;
	// matrix type
	typedef typename PStruct::Matrix Matrix;

	enum{n=PStruct::n}; // number of switches

	using Vector = BlockVector<FieldVector<double,1>>; 

	typedef ActiveSet<PStruct> ActiveSetSolver; 
	typedef SwitchPoly BaseSwitchSet; 

	typedef typename PStruct::AdjProblem AdjProblem;
	typedef typename PStruct::CGProblem CGProblem;
 
private: 
	PStruct& pstruct; // parameter structure
	const vector<vtype>& Sf; // encodes Sf
	const BaseSwitchSet& polytope; // switching polytope
	
public: 
	OutputData& data; // storage for iteration results

	OuterApprox(PStruct& pstruct_, const vector<vtype>& Sf_, //encodes Sf
		const BaseSwitchSet& set, // switching polytope 
		OutputData& data_) // storage for iteration result
	: pstruct(pstruct_),  
	Sf(Sf_), 
	polytope(set),
	data(data_)
	{}
	 
	double objective(ControlVector& u, vector<vtype>& y) const;
	vector<vtype> apply(const ControlVector& Gstar_eta, ControlVector& u, string dir) const;

};


/* objective value calculation 
*/
template<class PStruct> 
double OuterApprox<PStruct>::objective(ControlVector& u, vector<vtype>& y) const
{
	
	vector<vtype> g(y.size(),vtype(pstruct.gfs));
	g=y; 
	for(int i=0; i < pstruct.dt.size(); i++) 
		g[i]+=this->Sf[i];
    
    	double bound=0; 
    	// distribution 1/2 int_\Omega \int_[0,T] (y+Sf-y_d)^2 dt dx
	L2Deviation(pstruct.gfs,pstruct.dt,g,pstruct.YdProblem,bound);

    	// distribution alpha/2 \int_[0,T] ||u-u_d||^2 dt 
    	for (int i=0; i < pstruct.dt.size(); i++){
		double time=accumulate(pstruct.dt.begin(),pstruct.dt.begin()+i+1,0.0);
		auto diff=u[i]-pstruct.heat.ud(time); 
       		bound+=diff.two_norm2()*pstruct.alpha*pstruct.dt[i]*0.5; 
	}

    	return bound;
}




/* outer approximation algorithm 
*/
template<class PStruct> 
vector<typename PStruct::vtype> OuterApprox<PStruct>::apply(
const ControlVector& Gstar_eta,  // encodes G*S*(Sf-yd)
ControlVector& u, // start control
string dir) const // dir to ouput solutions const 
{

	ActiveSetSolver assolver=ActiveSetSolver(pstruct.gfs,pstruct.dt, pstruct.heat, Gstar_eta, pstruct.iter_activeset, pstruct.iter_lin, pstruct.reduction, pstruct.alpha, pstruct.rho, pstruct.tol);

	ControlVector u_=u;
	vector<vtype> y_; 
	// initialize SGu
	vector<vtype> y; 
	// initialize S*SGu
	ControlVector Gstar_zeta(pstruct.dt.size(),0);
		
	int count =0; 

	// initialize data for cutting planes
   	Matrix A; 
    	Vector b;
    	Vector mu;

   	double time=0.0; // overall runtime in seconds
	double obj; // objective of previous iteration
	
	cout << "\n---------------------- Outer Approximation Algorithm ----------------------\n" << endl; 
	
	// create ActiveSetResult;
	ActiveSetResult res;
    	clock_t cstart=clock();
	// apply solver 
	assolver.apply(u_, y_, Gstar_zeta, res); 
    	data.time_outer_solve.push_back((double)(clock()-cstart)/CLOCKS_PER_SEC); 
	data.iter_activeset+=res.iter_activeset;
	data.iter_lin+=res.iter_lin;
	data.time_lin_solve.push_back(res.time_solve);
	data.time_lin_assemble.push_back(res.time_assemble);
 
    	// accept iteration and output new control
	u=u_;
	y=y_;
	cout << "control" << ": " << endl;
	for(int i=0; i< u.size(); i++)
		cout << u[i] << endl;
	


	// calculate objective function value of the outer approximation iteration
    	double bound=objective(u,y);
    	data.lower_bound.push_back(bound); 

    	bool stop_outer = false; 
        
	time+=(double)(clock()-cstart)/CLOCKS_PER_SEC;

	int nt=pstruct.dt.size()*n;
	double* cut=new double[nt]; 
	double* v=new double[nt]; 
    	while(!stop_outer && count< pstruct.iter_outer){
		// write current control to file
		stringstream fcut;
		fcut << MIOCP_OUTPUT_PATH << dir << "cut_" << count << ".txt"; 
		fstream rescut; 
		rescut.precision(10);
		rescut.open(fcut.str(), ios::out | ios::trunc);
		rescut <<"[Grid]" << endl; 
		rescut << "nodes:" << pstruct.dt.size() << endl;
		for(int i=0; i < pstruct.dt.size(); i++)
			rescut << "dt_" << i << ":" << pstruct.dt[i] << endl;
		rescut << "[Control]" << endl;
		for(int i=0; i< u.size(); i++)
			rescut << "u_" << i << ":" << u[i] << endl;
		rescut << "[Objective]" << endl;
		rescut << bound << endl;
		rescut.close();
		// calculate heuristic solution (if desired)
#ifdef HEURISTIC
			cout << "\n=== Calculation of heuristic solution" << endl;
			double* c=new double[nt];
			double* w=new double[nt];
			for(int i=0; i<pstruct.dt.size(); i++){
				for(int j=0; j < n; j++){
					w[i*n+j]=u[i][j];
					c[i*n+j]=0.5-w[i*n+j];
				}
			}	
			polytope.optimize(nt,c,w); 

			ControlVector hsol(pstruct.dt.size()); 
			for(int i=0; i<pstruct.dt.size(); i++){
				for(int j=0; j < n; j++)
					hsol[i][j]=w[i*n+j];
			}
			CGProblem cgproblem(hsol,pstruct.dt);
			vector<vtype> Sh;
			heatdriver(pstruct.gfs,cgproblem,pstruct.dt,Sh);
			double hval=objective(hsol,Sh);
			stringstream hfile; 
			hfile << MIOCP_OUTPUT_PATH << dir << "hsol_"<<count << ".txt";
			fstream file; 
			file.precision(10);
			file.open(hfile.str(),  ios::out | ios::trunc); 
			file <<"[Grid]" << endl; 
			file << "nodes:" << pstruct.dt.size()<<endl;
			for(int i=0; i < pstruct.dt.size(); i++)
				file << "dt_" << i << ":" << pstruct.dt[i] << endl;
			file << "[Control]" << endl;
			for(int i=0; i < pstruct.dt.size(); i++)
				file << "u_" << i << ":" << hsol[i] << endl;
			file << "[Objective]" << endl;
			file << hval << endl;
			file.close();
			delete[] w;
			delete[] c;
#endif

		// store current control
		for(int i=0; i<pstruct.dt.size(); i++){
			for(int j=0; j < n; j++)
				v[i*n+j]=u[i][j];
		}
		// store current objective 
		obj=bound; 
		double rhs;

		// separation of control
		cstart=clock();
		int status = polytope.separate(nt,v, cut, rhs); 
        	data.time_outer_separation.push_back((double)(clock()-cstart)/CLOCKS_PER_SEC);
		time+=(double)(clock()-cstart)/CLOCKS_PER_SEC;
        	if(status==0){ // control feasible: stop outer approximation
            		stop_outer = true; 
			data.status = 0; 
        	}
        	else{ // control infeasible: add cutting plane and start next outer approximation iteration
	       	 	count++;
			cout << "\n=== Outer Approximation Iteration " << count << endl; 
 
			// append cutting plane to matrix
            		Matrix A_new=Matrix(A.N()+1,nt,Matrix::row_wise); 
            		for(auto row=A_new.createbegin(); row!=A_new.createend(); ++row){
                		if (row.index()< A.N()){
			        	for( auto it=A[row.index()].begin(); it!=A[row.index()].end(); ++it) 
						row.insert(it.index());
                		}
               			else{
                    			for(int k=0; k<nt; k++){
                        			if(cut[k]!=0) 
                            				row.insert(k);  
                			}
            			} 
			}

            		for(int i=0; i< A.N(); i++) 
                		A_new[i]=A[i];
            		int end = A.N();
            		for(int k=0; k<nt; k++){
				if(cut[k]!=0)
                    			A_new[end][k]=cut[k]; 
			}
			A_new[end]*=pstruct.alpha;

			// extend right hand side
			Vector b_new(end+1);
			for(int i=0; i < end; i++) 
				b_new[i]=b[i]; 
			b_new[end]=rhs;
			b_new[end]*=pstruct.alpha;
			
			// extend lagrange multiplicator 
			Vector mu_new(end+1);
			for(int i=0; i < end; i++)
            			mu_new[i]=mu[i];
		
			// start outer approximation iteration
			res = ActiveSetResult(); // reset data
            		clock_t ostart=clock();
			assolver.apply(u_, y_, Gstar_zeta, mu_new, A_new, b_new, res, pstruct.reopt,pstruct.timeouter-time);
			time+=(double)(clock()-ostart)/CLOCKS_PER_SEC;

			if(res.status==0){ // accept iteration and store output data
				data.time_outer_solve.push_back((double)(clock()-ostart)/CLOCKS_PER_SEC); 
				data.iter_activeset+=res.iter_activeset;
				data.iter_lin+=res.iter_lin;
				data.time_lin_solve.push_back(res.time_solve);
				data.time_lin_assemble.push_back(res.time_assemble);
		

            			// accept iteration
				u=u_;
				y=y_;
            			A=A_new;
				mu=mu_new; 
				b=b_new;
            
            			// output new control
            			cout << "control u_" << count << ": " << endl;
				for(int i=0; i< u.size(); i++) 
					cout << u[i] << endl; 

            			// calculate objective function value of the outer approximation iteration
				bound=objective(u,y);
            			data.lower_bound.push_back(bound); 
			}
			else if(res.status==-1){ // time limit reached
				data.status=-2;
				stop_outer=true;
			}
			else{ // maximum iterations of active set reached
				data.status=-3;
				stop_outer=true;
			}
       		}
    	}
    	cout << "\n---------------------- End Outer Approximation Algorithm ----------------------\n" << endl;

    	if(count>= pstruct.iter_outer) 
		data.status=-1;

	// save number of outer approximation iterations
    	data.iter_cut = count; 

	delete[] v;
	delete[] cut;

	return y;
}

#endif //ifned OUTER_APPROX_HH
