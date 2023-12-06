// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
// dune istl includes
#include <dune/istl/operators.hh> 
#include<dune/istl/vbvector.hh>
#include <dune/istl/solvers.hh> 
#include <dune/istl/umfpack.hh>


#include <dune/MIOCP/hfunctions.hh>

using namespace Dune; 
using namespace std; 

/* Linear Operator without cutting planes 
*/
template<class PStruct>
class LinOperator : public LinearOperator< BlockVector<FieldVector<double,1>>, BlockVector<FieldVector<double,1>> >
{
protected:
	typedef BlockVector<FieldVector<double,1>> X;
	
	typedef typename PStruct::GFS GFS; 
    	typedef typename PStruct::vtype vtype; // type of dof vector
	typedef typename PStruct::ControlVector ControlVector; // type of control vector

	typedef typename PStruct::AdjProblem AdjProblem;
	typedef typename PStruct::CGProblem CGProblem;

public:
	LinOperator(const GFS& gfs_,const vector<double>& dt_,double alpha_, const vector<pair<int,int>>& I_) : gfs(gfs_), dt(dt_),alpha(alpha_), I(I_) {};

	/* y=Cx */
	void apply(const X& x, X& y) const 
	{
		// clear y
		y=0.0;
		
		// fill x with zeros
		ControlVector z(dt.size(),0);
		for (int i=0;i<I.size(); i++ )
			z[I[i].first][I[i].second]=x[i];
 
		// calculation of w=SGz;
		CGProblem cgproblem(z,dt); 
		vector<vtype> w; 
		heatdriver(gfs,cgproblem,dt,w);
    
    	
		// calculation of v=S*SGz;
    		AdjProblem adjointProblem=AdjProblem(dt,w,gfs); 
    		vector<vtype> v; 
		adjointdriver(gfs, adjointProblem,v);  

    		
    		// y= \alpha x + \Chi_I(G*S*SGz) 
		y=x;
        	y*=alpha;
    		// contibution of \Chi_I(G*S*SGz) = \Chi_I(\int Psi(x)v(t,x) dx )
		ControlVector Gstar_z(z.size()); 
    		Gstar(gfs,cgproblem,v,I,Gstar_z);
		for (int i=0;i<I.size(); i++ ){
			y[i]+=Gstar_z[I[i].first][I[i].second];	
			y[i]*=dt[I[i].first];
		}	

	}

	/* y= y+alpha*Cx */
	void applyscaleadd(field_type alpha, const X& x, X& y) const
	{
		X z=y; // auxiliary variable to save result Cx
		apply(x,z); 
		
		z*=alpha; 
		y+=z;	
	}

	SolverCategory::Category category() const{
		return SolverCategory::sequential;
	}


	private: 
	const GFS& gfs; // spatial grid function space
	const vector<double>& dt; // temporal grid
	double alpha; // Thikonov term
	const vector<pair<int,int>>& I; // indices of inactive box constraints
	
	
};


/* Linear operator with cutting planes 
*/
template<class PStruct>
class LinCutOperator : public LinearOperator<BlockVector<FieldVector<double,1>>,BlockVector<FieldVector<double,1>>>
{
protected:
	typedef typename PStruct::ControlVector ControlVector; // type of control vector
	typedef typename PStruct::Matrix Matrix; // matrix type

	typedef typename PStruct::GFS GFS;
    	typedef typename PStruct::vtype vtype; // type of dof vector

	typedef typename PStruct::AdjProblem AdjProblem;
	typedef typename PStruct::CGProblem CGProblem;

	typedef BlockVector<FieldVector<double,1>> X;

	enum{n=PStruct::n}; // number of switches


public:
	LinCutOperator(const GFS& gfs_,const vector<double>& dt_,double alpha_, const vector<pair<int,int>>& I_,const Matrix& A_) : gfs(gfs_), dt(dt_),alpha(alpha_), I(I_), A(A_) {};

	/* y=Cx */
	void apply(const X& x, X& y) const 
	{
		// clear y
		y=0.0;
	
		// fill x[0] with zeros
		ControlVector z(dt.size());
		for (int i=0;i<I.size(); i++ )
			z[I[i].first][I[i].second]=x[i];

		// calculation of w=SGx[0]
		CGProblem cgproblem(z,dt); 
		vector<vtype> w; 
		heatdriver(gfs,cgproblem,dt,w);
    
    		
		// calculation of v=S*SGx[0];
    		AdjProblem adjointProblem=AdjProblem(dt,w,gfs); 
    		vector<vtype> v; 
		adjointdriver(gfs,adjointProblem,v);  

		// y[0]= \alpha x[0] + \Chi_I(G*S*SGx[0]+A^Tx[1])
		for (int i=0;i<I.size(); i++){
			y[i]=x[i];
        		y[i]*=alpha;
		} 
    		// contibution of \Chi_I(G*S*SGx[0]) = \Chi_I(\int Psi(x)v(t,x) dx )
		ControlVector Gstar_z(dt.size()); 
    		Gstar(gfs,cgproblem,v,Gstar_z);
		for (int i=0;i<I.size(); i++ ){
			y[i]+=Gstar_z[I[i].first][I[i].second];
			y[i]*=dt[I[i].first];
		}
		// calculate A^Tx[1] and then add only \Chi_I(A^Tx[1])
		X h1(A.M()); // auxiliary variable
		X h2(A.N());
		for (int i=0;i<A.N(); i++)
			h2[i]=x[i+I.size()];
		// h1=A^Th2
		A.mtv(h2,h1);
		for (int i=0;i<I.size(); i++ )
			y[i]+=h1[I[i].first*n+I[i].second];
		
		// y[1]=Az (z- x[0] vector filled with zeros)
		h2.resize(A.N());
		h2=0; 
		// store z appropiately
		for(int i=0; i < dt.size(); i++){
			for(int j=0; j < n; j++)
				h1[i*n+j]=z[i][j];
		}
		// calculate Az
		A.mv(h1,h2);
		for (int i=0;i<A.N(); i++ )
			y[i+I.size()]=h2[i];

	}
	/* y=y+alpha*Cx */
	void applyscaleadd(field_type alpha, const X& x, X& y) const
	{
		X z=y; // auxiliary variable to save result Cx 
		apply(x,z); 
		
		z*=alpha; 
		y+=z;	
	}

	SolverCategory::Category category() const{
		return SolverCategory::sequential;
	}

	private: 
	const GFS& gfs; // spatial grid function space
	const vector<double>& dt; // temporal grid
	double alpha; // Thikonov term
	const vector<pair<int,int>>& I; // indices of inactive box constraints
	const Matrix& A; // matrix of active cutting planes
	
};

/* Inverse Operator for Preconditioning */ 
template<class M,class X,class Y> 
class InverseP : public InverseOperator<X,Y>{
  public:	
    	typedef M matrix_type;
 	typedef X domain_type;
    	typedef Y range_type;

	// Constructor 
	InverseP(const M& P_): P(P_) {}; 

	// apply inverse operator
	void apply(X& x,Y& b, InverseOperatorResult& res){
		// call UMFPack
		UMFPack<M> Solver(P); 
		Solver.apply(x,b,res);
	}

	// apply inverse operator, with given convergence criteria
	void apply(X& x,Y& b,double reduction, InverseOperatorResult& res){
		// call UMFPack
		UMFPack<M> umfpackSolver(P); 
		umfpackSolver.apply(x,b,reduction,res);
	}

	// Solver category 
	SolverCategory::Category category() const{
		return SolverCategory::sequential;
	}
  private: 
	const M& P; // Preconditioning matrix
};



/* modified MINRRESSolver with additional termination criteria */
template<class X> 
class OwnMINRESSolver : public MINRESSolver<X> 
{
private: 
	typedef MINRESSolver<X> BaseType;
	using typename BaseType::IterativeSolver::scalar_real_type; 
	using typename BaseType::real_type; 
	using typename BaseType::field_type; 
	using Iteration = typename BaseType::Iteration;

	using BaseType::_op; 
	using BaseType::_prec; 
	using BaseType::_sp;
	using BaseType::_maxit;  
public:

	 // copy base class constructors
    	using MINRESSolver<X>::MINRESSolver;

	// copy apply and only add new stopping criteria
	void apply (X& x, X& b,double tol, InverseOperatorResult& res)
    	{
      		using std::sqrt;
      		using std::abs;
      		Iteration iteration(*this, res);
      		// prepare preconditioner
      		_prec->pre(x,b);

		real_type nb=_sp->norm(b);
      		// overwrite rhs with defect
     		_op->applyscaleadd(-1,x,b);

      		// compute residual norm
      		real_type def = _sp->norm(b);
      		if(iteration.step(0, def)|| Simd::max(def)<tol){
        		_prec->post(x);
        		return;
      		}

      		// recurrence coefficients as computed in Lanczos algorithm
      		field_type alpha, beta;
       	 	// diagonal entries of givens rotation
      		std::array<real_type,2> c{{0.0,0.0}};
        	// off-diagonal entries of givens rotation
      		std::array<field_type,2> s{{0.0,0.0}};

      		// recurrence coefficients (column k of tridiag matrix T_k)
      		std::array<field_type,3> T{{0.0,0.0,0.0}};

      		// the rhs vector of the min problem
      		std::array<field_type,2> xi{{1.0,0.0}};

      		// some temporary vectors
      		X z(b), dummy(b);

      		// initialize and clear correction
      		z = 0.0;
      		_prec->apply(z,b);

      		// beta is real and positive in exact arithmetic
      		// since it is the norm of the basis vectors (in unpreconditioned case)
      		beta = sqrt(_sp->dot(b,z));
      		field_type beta0 = beta;

      		// the search directions
      		std::array<X,3> p{{b,b,b}};
      		p[0] = 0.0;
      		p[1] = 0.0;
      		p[2] = 0.0;

      		// orthonormal basis vectors (in unpreconditioned case)
      		std::array<X,3> q{{b,b,b}};
      		q[0] = 0.0;
      		q[1] *= real_type(1.0)/beta;
      		q[2] = 0.0;

      		z *= real_type(1.0)/beta;

      		// the loop
      		int i = 1;
      		for( ; i<=_maxit; i++) {

        		dummy = z;
        		int i1 = i%3,
          		i0 = (i1+2)%3,
          		i2 = (i1+1)%3;

        	// symmetrically preconditioned Lanczos algorithm (see Greenbaum p.121)
        	_op->apply(z,q[i2]); // q[i2] = Az
        	q[i2].axpy(-beta,q[i0]);
        	// alpha is real since it is the diagonal entry of the hermitian tridiagonal matrix
        	// from the Lanczos Algorithm
        	// so the order in the scalar product doesn't matter even for the complex case
        	alpha = _sp->dot(z,q[i2]);
        	q[i2].axpy(-alpha,q[i1]);

        	z = 0.0;
        	_prec->apply(z,q[i2]);

        	// beta is real and positive in exact arithmetic
        	// since it is the norm of the basis vectors (in unpreconditioned case)
       		beta = sqrt(_sp->dot(q[i2],z));

        	q[i2] *= real_type(1.0)/beta;
        	z *= real_type(1.0)/beta;

        	// QR Factorization of recurrence coefficient matrix
        	// apply previous givens rotations to last column of T
        	T[1] = T[2];
        	if(i>2) {
          		T[0] = s[i%2]*T[1];
          		T[1] = c[i%2]*T[1];
        	}
        	if(i>1) {
         		T[2] = c[(i+1)%2]*alpha - s[(i+1)%2]*T[1];
          	T[1] = c[(i+1)%2]*T[1] + s[(i+1)%2]*alpha;
        	}
        	else
          	T[2] = alpha;

        	// update QR factorization
        	generateGivensRotation(T[2],beta,c[i%2],s[i%2]);
        	// to last column of T_k
        	T[2] = c[i%2]*T[2] + s[i%2]*beta;
        	// and to the rhs xi of the min problem
        	xi[i%2] = -s[i%2]*xi[(i+1)%2];
        	xi[(i+1)%2] *= c[i%2];

        	// compute correction direction
        	p[i2] = dummy;
        	p[i2].axpy(-T[1],p[i1]);
        	p[i2].axpy(-T[0],p[i0]);
        	p[i2] *= real_type(1.0)/T[2];

        	// apply correction/update solution
        	x.axpy(beta0*xi[(i+1)%2],p[i2]);

        	// remember beta_old
        	T[2] = beta;

        	// check for convergence
        	// the last entry in the rhs of the min-problem is the residual
        	def = abs(beta0*xi[i%2]);
        	if(iteration.step(i, def)|| Simd::max(def)<tol){ // additional termination criteria
         	 	break;
        	}
      	} 

        // postprocess preconditioner
        _prec->post(x);
    }

private: 
	void generateGivensRotation(field_type &dx, field_type &dy, real_type &cs, field_type &sn)
    	{
      		using std::sqrt;
      		using std::abs;
      		using std::max;
      		using std::min;
      		const real_type eps = 1e-15;
      		real_type norm_dx = abs(dx);
      		real_type norm_dy = abs(dy);
      		real_type norm_max = max(norm_dx, norm_dy);
     		real_type norm_min = min(norm_dx, norm_dy);
      		real_type temp = norm_min/norm_max;
      		// we rewrite the code in a vectorizable fashion
      		cs = Simd::cond(norm_dy < eps,
        		real_type(1.0),
        		Simd::cond(norm_dx < eps,
         	 real_type(0.0),
         	 Simd::cond(norm_dy > norm_dx,
            	real_type(1.0)/sqrt(real_type(1.0) + temp*temp)*temp,
            	real_type(1.0)/sqrt(real_type(1.0) + temp*temp)
          	)));
      		sn = Simd::cond(norm_dy < eps,
       	 	field_type(0.0),
        	Simd::cond(norm_dx < eps,
          	field_type(1.0),
          	Simd::cond(norm_dy > norm_dx,
           	 // dy and dx are real in exact arithmetic
            	// thus dx*dy is real so we can explicitly enforce it
            	field_type(1.0)/sqrt(real_type(1.0) + temp*temp)*dx*dy/norm_dx/norm_dy,
            	// dy and dx is real in exact arithmetic
            	// so we don't have to conjugate both of them
            	field_type(1.0)/sqrt(real_type(1.0) + temp*temp)*dy/dx
          	)));
    	}
};


/* modified CG solver with additional termination criteria */
template<class X>
class OwnCGSolver : public CGSolver<X>
{
private: 
	typedef CGSolver<X> BaseType; 

	using typename BaseType::IterativeSolver::scalar_real_type; 
	using typename BaseType::real_type; 
	using typename BaseType::field_type; 
	using Iteration = typename BaseType::Iteration;

	using BaseType::_op; 
	using BaseType::_prec; 
	using BaseType::_sp;
	using BaseType::_maxit;  
public: 
	 // copy base class constructors
	using CGSolver<X>::CGSolver;

	// copy apply method and only add stopping criteria
	void apply (X& x, X& b,double tol, InverseOperatorResult& res)
    	{
      		Iteration iteration(*this,res);
      		_prec->pre(x,b);             // prepare preconditioner
		
		real_type nb=_sp->norm(b);
      		_op->applyscaleadd(-1,x,b);  // overwrite b with defect

      		real_type def = _sp->norm(b); // compute norm
      		if(iteration.step(0, def)|| Simd::max(def)<tol){
       			_prec->post(x);
        		return;
      		}

      		X p(x);              // the search direction
      		X q(x);              // a temporary vector

      		// Remember lambda and beta values for condition estimate
      		std::vector<real_type> lambdas(0);
      		std::vector<real_type> betas(0);

      		// some local variables
      		field_type rho,rholast,lambda,alpha,beta;

      		// determine initial search direction
      		p = 0;                          // clear correction
      		_prec->apply(p,b);               // apply preconditioner
      		rholast = _sp->dot(p,b);         // orthogonalization

      		// the loop
     		int i=1;
     		for ( ; i<=_maxit; i++ ){
       			// minimize in given search direction p
        		_op->apply(p,q);             // q=Ap
        		alpha = _sp->dot(p,q);       // scalar product
        		lambda = rholast/alpha;     // minimization
				
        		x.axpy(lambda,p);           // update solution
        		b.axpy(-lambda,q);          // update defect
    
        		// convergence test
       	 		def=_sp->norm(b); // comp defect norm
        		if(iteration.step(i, def) || Simd::max(def)<tol) // additional termination criteria
          			break;

        		// determine new search direction
        		q = 0;                      // clear correction
        		_prec->apply(q,b);           // apply preconditioner
        		rho = _sp->dot(q,b);         // orthogonalization
        		beta = rho/rholast;         // scaling factor
        	
       	 		p *= beta;                  // scale old search direction
        		p += q;                     // orthogonalization with correction
       			rholast = rho;              // remember rho for recurrence
      		}

      			_prec->post(x);                  // postprocess preconditioner

	}
	
};

