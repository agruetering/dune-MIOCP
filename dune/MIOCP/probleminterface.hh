// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef PROBLEMINTERFACE_HH
#define PROBLEMINTERFACE_HH

using namespace Dune; 
using namespace std; 


/* Problem interfaces for optimal control problems of the form
 *
 *    min J(y,u)=1/2|y-y_d|_{L^2(Q)}^2 + \apha/2 |u-u_d|_{L^2(0,T;R^n)}^2
 *    s.t.  \partial_t y - \Delta y = Gu + f = \sum_j u_j(t)*Psi_j(x) + f	in Q=\Omega x (0,T)
 *                       		   y = g   				on \Sigma_D=\Gamma_D x (0,T)
 *         		    \nabla y \cdot n = j   				on \Sigma_N=\Gamma_N x (0,T)
 *                        		y(0) = y0 		
 * 					u_a <= u <= u_b 			a.e. in (0,T)
 * 
 * 
 *
 */ 

/* ProblemInterface (boundary data, initial data and right hand side f) */
template<class GridView>
class ProblemInterface
{
protected:
	
	typedef GridView GV; // spatial grid view
	enum {dim=GV::dimension}; // dimension of spatial grid
	enum {interdim=GV::dimension-1}; // facet dimension

	typedef typename GV::Traits::template Codim<0>::Entity E;  // elements type
	typedef typename GV::Traits::Intersection I; // facet type
	typedef typename GV::ctype ctype; // coordinate type of the spatial grid
	typedef FieldVector<ctype,dim> X; // vector type of element coordinates
	typedef FieldVector<ctype,interdim> Y; // vector type of facet coordinates

	double t; //current time

public:

  	// Constructor
  	ProblemInterface () :  t(0.0) {}

	// right hand side contribution q(=f) 
	virtual double q(const E& e, const X& x) const
	{
		return 0.0; 
	}


  	// boundary condition type function (true = only Dirichlet, i.e., \Sigma_D=\partial\Omega x (0,T) )
  	virtual bool b (const I& i, const Y& x) const
  	{
   		 return true;
  	}

  	// Dirichlet boundary function g (default 0)
  	virtual double g (const E& e, const X& x) const
  	{
   		return 0.0;
  	}

  	// Neumann boundary function j (default 0)
	virtual double j (const I& i, const Y& x) const
  	{
    		return 0.0;
  	}

	// initial data y0
  	virtual double ini(const E& e, const X& x) const
	{
		return 0.0;
	}

	// Set time 
	void setTime(double t_)
  	{
		 t=t_;
  	}
}; 


/* Problem to calculate right hand side contribution 
 * q(=f) (with boundary and initial data)
 */ 

template<class GridView, int number>
class HeatProblemInterface : public ProblemInterface<GridView> 
{
protected: 
	enum {n=number}; // number of switches

	typedef ProblemInterface<GridView> BaseType;

	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::t;


public:
	HeatProblemInterface() : BaseType() {}

	// form functions (need to be specified)
  	virtual FieldVector<double,n> Psi (const E& e, const X& x) const =0;


	// lower bounds on controls
	virtual FieldVector<double,n> ua(double t_) const 
	{
		return FieldVector<double,n>(0.0);
	}

	// upper bounds on controls
	virtual FieldVector<double,n> ub(double t_) const
	{
		return FieldVector<double,n>(1.0);
	}

	// desired control
	virtual FieldVector<double,n> ud(double t_) const
	{
		return FieldVector<double,n>(0.5);
	}

	// check if problem is homogeneous 
	virtual bool isHomogeneous() {return true;}

	
};

/* problem interface for calculating SGu */
template<class GridView, int number> 
class CGProblemInterface : public ProblemInterface<GridView>
{
protected: 
	enum {n=number}; // number of switches

	typedef BlockVector<FieldVector<double,n>> ControlVector; // type of control vector

	const ControlVector& u; // control vector
	vector<double> dt; // step wides

	typedef ProblemInterface<GridView> BaseType;

	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::t;

public:
	// constructor equidistant time stepping
	CGProblemInterface(const ControlVector& w, double T,int timesteps) : u(w), dt(vector<double>(timesteps,T/timesteps)), BaseType() {}
	// constructor
	CGProblemInterface(const ControlVector& w, vector<double> dt_) : u(w), dt(dt_), BaseType() {}
	
	// form functions (need to be specified)
  	virtual FieldVector<double,n> Psi (const E& e, const X& x) const =0;

	// right hand side q(=Gu)
  	double q(const E& e, const X& x) const
	{
		// determine index i (corresponding to time t)
		int i=0;
		double time=0; 
		while(time+dt[i] < t-1e-10){
			time+=dt[i];
			i++; 
		}
		FieldVector<double,n> psi=Psi(e,x);
		return u[i].dot(psi);
	}

};


/* Problem interfaces for adjoint equations of the form
 * -\partial_t p - \Delta p = q	 	in Q=\Omega x (0,T)
 *                        p = g  	on \Sigma_D=\Gamma_D x (0,T)
 *         \nabla p \cdot n = j 	on \Sigma_N=\Gamma_N x (0,T)
 *                     p(T) = pT 
 */ 

/* Adjoint problem interface (boundary data, end data and right hand side) */
template<class GridView>
class AdjointProblemInterface : public ProblemInterface<GridView>
{
protected:
	double T; //end time
	vector<double> dt; // step wides
	int i; // index for current step wide

	typedef ProblemInterface<GridView> BaseType;

	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::t;

public:

	// constructor equidistant time stepping
	AdjointProblemInterface (double endTime_, int timesteps) : T(endTime_),  dt(vector<double>(timesteps,T/timesteps)), i(0), BaseType() {}
  	// constructor
  	AdjointProblemInterface (vector<double> dt_) : dt(dt_), i(0), BaseType() {T=accumulate(dt.begin(),dt.end(), 0.0);}

	// default constructor
	AdjointProblemInterface () : T(1), dt(vector<double>(100,1e-2)), BaseType() {}

	// right hand side at time point t_l (need to be specified)
	virtual double rhs (const E& e, const X& x, int l) const = 0;

	// right hand side contribution (given by the temporal discretization)
	double q (const E& e, const X& x) const 
	{	
		double value=0.0;
		int end=dt.size()-1;
    		if(i==end){ // 1/6 (2*yd(t)+ yd(t+dt(end)))
			value+=rhs(e,x,i); 
			value*=2.0;
			value+=rhs(e,x,i+1);
			value*=1.0/6.0; 
				 
		}  
		else{ // 1/6 (2(dt(i+1)/dt(i)+1)*yd(t+dt(i))+ yd(t)+ dt(i+1)/dt(i)*yd(t+dt(i)+dt(i+1)))
			value+=2.0*(dt[i+1]/dt[i]+1)*rhs(e,x,i+1); 
			value+=rhs(e,x,i);
			value+=dt[i+1]/dt[i]*rhs(e,x,i+2);
			value*=1.0/6.0;
		}
		return value; 
	}

	// get end time
	double endTime()
	{
		return T;
	}

	// change step wides 
	void setdeltaT(vector<double> dt_)
	{
		dt=dt_;
	}

	
	// get step wide at specific index
	double deltaT(int count)
	{
		return dt[count];
	}


	// get number of time steps 
	double timesteps(){
		return dt.size();
	}

	// Set time and index for current step wide
	void setTime (double t_, int count)
  	{
		 t=t_;
		 i=count;
  	}

	using BaseType::setTime;
	
};

/* problem interface to calculate S*SGu */
template<class GFS>
class AdjointProblem : public AdjointProblemInterface<typename GFS::Traits::GridView>
{
protected:
	typedef AdjointProblemInterface<typename GFS::Traits::GridView> BaseType; 

	typedef typename PDELab::Backend::Vector<GFS,double> Y; // type of Dof vectors

	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::T;
	using BaseType::dt;
	using BaseType::i;
	using BaseType::t;

	vector<Y>& y; // Dof vector y(t)=SGv(t) at every time point t
	const GFS& gfs; // grid function space 

public: 
 	AdjointProblem(double endTime_, int timesteps, vector<Y>& y_,const GFS& gfs_)  : y(y_), gfs(gfs_), BaseType(endTime_, timesteps)  {}
	AdjointProblem(vector<double> dt_,vector<Y>& y_,const GFS& gfs_)  : y(y_), gfs(gfs_), BaseType(dt_)  {}

	
	double rhs (const E& e, const X& x, int l) const
  	{
		typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF;
    		typedef typename DGF::Traits::RangeType RangeType;

		// return value y(t_l,x)
		DGF ydgf(gfs,y[l]);
		RangeType value; 
    		ydgf.evaluate(e,x,value); 
		return value;  
  	}

};


#endif // ifndef PROBLEMINTERFACE_HH
