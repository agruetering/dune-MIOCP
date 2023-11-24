#include "adjointheatfem.hh"
#include<sys/stat.h>

using namespace Dune; 
using namespace std; 

/* scheme to solve adjoint equation 
 *
 *  -\partial_t p - \Delta p = q   in Q=\Omega x (0,T)
 *                         p = g   on \Sigma_D=\Gamma_D x (0,T)
 *          \nabla y \cdot n = j   on \Sigma_N=\Gamma_N x (0,T)
 *                      p(T) = pT
 *
 */

template<class GFS>
void adjointdriver (const GFS& gfs, // grid function space
AdjointProblemInterface<typename GFS::Traits::GridView>& problem, // problem class
vector<typename PDELab::Backend::Vector<GFS,double>>& sol)  // vector to solve p(t)
{  
#ifdef NEWTONDEBUG
	cout << "-- solve backward equation --" << endl;
#endif
	typedef typename GFS::Traits::GridView GV; 
	typedef typename GFS::Traits::FiniteElementMap FEM; 
	typedef typename PDELab::Backend::Vector<GFS,double> Y;
	typedef typename PDELab::DiscreteGridFunction<GFS,Y> YDGF; 
	typedef AdjointProblemInterface<GV> Problem;

  	// dimension of the spatial grid 
  	const int dim = GV::dimension;  

 	 //////////////////  get problem data and set time //////////////////
	// end time 
  	double time = problem.endTime();
	problem.setTime(time); 
	
	// lambda function for dirichlet boundary function
 	auto glambda = [&](const auto& e, const auto& x)
    	{
		return problem.g(e,x);
	};
	// make grid function for dirichlet boundary values
 	auto g = PDELab::makeInstationaryGridFunctionFromCallable(gfs.gridView(),glambda,problem);

	// lambda function for dirichlet boundary condition 
 	 auto blambda = [&](const auto& i, const auto& x)
   	{
		return problem.b(i,x);
	};
	// make boundary conditions on grid
 	auto b = PDELab::makeBoundaryConditionFromCallable(gfs.gridView(),blambda);

  	
 	// Assemble dirichlet constraints (determine constrained dofs)
  	typedef typename GFS::template ConstraintsContainer<double>::Type CC;
  	CC cc;
 	PDELab::constraints(b,gfs,cc); // determine constrained dofs
  	

	// lambda function for end function 
 	auto inilambda = [&](const auto& e, const auto& x)
    	{
		return problem.ini(e,x);
	}; 
	auto initial = PDELab::makeInstationaryGridFunctionFromCallable(gfs.gridView(),inilambda,problem);
  	//////// Dof vector ///////////
	// initialization
  	Y p(gfs); 
  	// Make a grid function out of it
  	YDGF pdgf(gfs,p);
	// fill the Dof vector 
	PDELab::interpolate(initial,gfs,p);



  	//////////////////  Make instationary grid operator //////////////////

	// spatial local operator (volume, right hand side and neuman boundary term)
  	typedef LinearHeatFEM<Problem,FEM> LOP;
  	LOP lop(problem);
	// matrix backend
  	typedef PDELab::ISTL::BCRSMatrixBackend<> MBE;
	// estimation of nonzero entries
  	MBE mbe((int)pow(-1+2*gfs.finiteElementMap().maxLocalSize(),dim));
	// grid operator
 	typedef PDELab::GridOperator<
    		GFS,GFS,   // ansatz and test space 
    		LOP,      // local operator
    		MBE,      // matrix backend 
   		double,double,double, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space
		> GO; 
  	GO go(gfs,cc,gfs,cc,lop,mbe);	


	// temporal local operator
  	typedef NL2<FEM> TLOP;
  	TLOP tlop;
	// temporal grid operator
  	typedef PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		TLOP,      // temporal local operator
    		MBE,      // matrix backend 
   		double,double,double, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space
		> TGO; 
  	TGO tgo(gfs,cc,gfs,cc,tlop,mbe);

	// grid operator for implicit euler method (one step method)
  	typedef PDELab::OneStepGridOperator<GO,TGO> IEGO;
  	IEGO iego(go,tgo);

	// Select a linear solver backend
  	typedef PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<IEGO> LS;
  	LS ls(5000,0);
	
  	// solve linear instationary problem
  	typedef Dune::PDELab::StationaryLinearProblemSolver<IEGO,LS,typename IEGO::Traits::Domain> PDESOLVER;
  	PDESOLVER pdesolver(iego,ls,1e-15,1e-99,0);

  	// prepare time-stepping scheme -> implicit euler method
 	PDELab::OneStepThetaParameter<double> method(1.0);
  	PDELab::TimeSteppingParameterInterface<double>* pmethod=&method;
  	
	// set up implicit euler method
  	typedef PDELab::OneStepMethod<double,IEGO,PDESOLVER,Y,Y> OSM;
  	OSM  osm(*pmethod,iego,pdesolver);
  	osm.setVerbosityLevel(0);

	////////////////// Time loop //////////////////

	sol.clear();
	sol.push_back(p); // save p(T)

	clock_t cstart=clock();
	int step=problem.timesteps();
 	for(int i=step;i>0;i--)
    	{
      		double dt=problem.deltaT(i-1);
      		// assemble constraints for new time step
     		problem.setTime(time-dt,i-1);
      		PDELab::constraints(b,gfs,cc);
	
     		 // do time step
      		Y pnew(p);
      		osm.apply(time,-dt,p,g,pnew);


      		// accept time step
      		p = pnew;
		sol.push_back(p); // save p(time-dt)
      		time-=dt;
    	}
#ifdef NEWTONDEBUG
	cout << " elapsed time: " << (double)(clock()-cstart)/CLOCKS_PER_SEC  << "\n" << endl; 
#endif
	 // save result in correct order (p(0),p(dt),...,p(T))
	reverse(sol.begin(),sol.end());

	
}
