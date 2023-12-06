// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// C, C++ includes
#include<math.h>
#include<iostream>
#include<vector>

// dune-common includes
#include<dune/common/parametertreeparser.hh>
#include<dune/common/fvector.hh> 
#include<dune/common/fmatrix.hh> 

// dune-grid includes
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/uggrid.hh> 
// dune-istl includes
#include<dune/istl/bvector.hh> 
#include<dune/istl/bcrsmatrix.hh>
// dune-pdelab includes
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include <dune/MIOCP/probleminterface.hh>
#include <dune/MIOCP/gurobi.hh> 

using namespace Dune; 
using namespace std; 

/* Problem to calculate right hand side contribution 
 * q(=f) (with boundary and initial data)
 */
template<class GridView>
class HeatProblem : public HeatProblemInterface<GridView,1>
{
protected:
	typedef HeatProblemInterface<GridView,1> BaseType;

	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::t; 

	
public: 
 	HeatProblem() : BaseType() {}

  	// form function Psi
  	FieldVector<double,1> Psi (const E& e, const X& x) const
  	{
		int dimDomain = x.size(); 
		auto global = e.geometry().global(x);
		double value =12*pow(M_PI,2);
	 	for(int i=0; i< global.size(); i++) 
			value*=exp(global[i])*sin(M_PI*global[i]);
		return value; 
  	}

 
};

/* problem to calculate desired temperature contribution yd (with boundary data) */ 
template<class GridView>
class AdjointDesired : public AdjointProblemInterface<GridView>
{
protected:
	typedef AdjointProblemInterface<GridView> BaseType;

	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::dt;

public:
	
 	AdjointDesired(double endTime_, int timesteps)  : BaseType(endTime_,timesteps) {}

	AdjointDesired(vector<double> dt_)  : BaseType(dt_) {}

	double rhs (const E& e, const X& x, int l) const
  	{
		double t=accumulate(dt.begin(), dt.begin()+l, 0.0);
		auto global = e.geometry().global(x); 
		double y_val = 2*pow(M_PI,2)* max(cos(2*M_PI*t),0.0);
		for(int i=0; i< global.size(); ++i)
	 		y_val*=sin(M_PI*global[i]); 
    		return y_val;  
  	}
 
};



/* problem for calculating SGu in CG algorithm */
template<class GridView>
class CGProblem : public CGProblemInterface<GridView,1>
{
protected:

	protected:
	typedef CGProblemInterface<GridView,1> BaseType;

	using typename BaseType::ControlVector;
	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

public:

	CGProblem(const ControlVector& w, double T,int timesteps) : BaseType(w,T,timesteps) {}

	CGProblem(const ControlVector& w, vector<double> dt_) : BaseType(w,dt_) {}

	FieldVector<double,1> Psi (const E& e, const X& x) const
  	{
		int dimDomain = x.size(); 
		auto global = e.geometry().global(x);   		
		double value = 12*pow(M_PI,2);
	 	for(int i=0; i< global.size(); i++) 
			value*=exp(global[i])*sin(M_PI*global[i]);
		return value; 
  	}
 
};


int main(int argc, char** argv)
{
    	
       //===============================================================
       // Make parameter structure 
       //===============================================================
    	// open ini file (contains grid information, parameter and folder name for output)
    	ParameterTree ptree;
    	ParameterTreeParser ptreeparser;
        stringstream file;
        file << MIOCP_SOURCE_PATH << "master-gurobi.ini"; 
    	ptreeparser.readINITree(file.str(),ptree);
    	ptreeparser.readOptions(argc,argv,ptree);

	
        //===============================================================
        // make 2D grid
        //===============================================================
        const int dim=2;
        typedef UGGrid<dim> Grid;	
        FieldVector<double,dim> L;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.L" << i; 
		L[i] = ptree.get<double>(name.str(),(double)0.0);
	}
	FieldVector<double,dim> U;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.U" << i; 
		U[i] = ptree.get<double>(name.str(),(double)1.0);
	}
        array<unsigned int,dim> N;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.N" << i; 
		N[i] = ptree.get<unsigned int>(name.str(),(unsigned int) 100);
	} 
	shared_ptr<Grid> gridp = StructuredGridFactory<Grid>::createSimplexGrid(L,U,N);

	typedef Grid::LeafGridView GridView;
        GridView gv=gridp->leafGridView();


	double T=ptree.get<double>("problem.T",(double)1.0); 
	int timesteps= ptree.get<int>("problem.Nt",(int)100);
	vector<double> dt(timesteps,T/timesteps);

	HeatProblem<GridView> heat; 
	AdjointDesired<GridView> AdjDesired(dt);

	//===============================================================
       	// Initialization
         //===============================================================


	int mod=0;
	/* mod overview 
	 * 0: integer problem with additional variables for tv norm
	 * 1: naive linear relaxation of tv norm
	 * 2: tailored convexification
	 * 3: no tv constraint
	 */
	  
	int smax=2;

	// write data to file
	stringstream outputfile;
	outputfile << MIOCP_OUTPUT_PATH << ptree.get("output.dir","output") << "NX" << N[0] << "x" << N[1] << "_step" <<timesteps << "_mod" << mod <<".txt";
	fstream stream;
	stream.precision(10);
	stream.open(outputfile.str(), ios::out | ios::app); 

	stream << "[sgrid]" << endl; 
	stream << "LX=" << L[0] << endl; 
	stream << "LY=" << L[1] << endl;
	stream << "UX=" << U[0] << endl; 
	stream << "UY=" << U[1] << endl;
	stream << "NX=" << N[0] << endl; 
	stream << "NY=" << N[1] << "\n" << endl;

	stream << "[tgrid]" << endl; 
	stream << "nodes=" << dt.size() << endl;
	for(int i=0; i < dt.size(); i++)
		stream << "dt_" << i << "=" << dt[i] << endl;

	stream << "[Parameter]" << endl;
	stream << "degree=" << 1 << endl;
	stream << "smax=" << 2 << endl; 


	stream.close();
	try{
		stringstream logfile; 
		logfile << MIOCP_OUTPUT_PATH << ptree.get("output.dir","output")  << "NX" << N[0] << "x" << N[1] << "_step" <<timesteps << "_mod" << mod << "_logfile.txt";
		
		// recalculation of objective within the gurobi optimizer
		// spatial triangulation with 100x100 nodes
		// temporal discretization with Nt=200
		array<unsigned int,dim> Nh;
		Nh.fill(100);
		shared_ptr<Grid> gpoint = StructuredGridFactory<Grid>::createSimplexGrid(L,U,Nh);
        	GridView grid=gpoint->leafGridView();
		int Nt=200;
		
		cout << "\n---------------------- Call Gurobi Optimizer ----------------------\n" << endl;
		typedef PDELab::ConformingDirichletConstraints CON; // conforming dirichlet constraints
		Gurobi<GridView,CON,CGProblem<GridView>,1> gsolver(gv,dt,heat,AdjDesired,grid,Nt); 	
		gsolver.optimize(mod, logfile.str(),outputfile.str(),smax);
		
	}
	catch (Exception &e){
		cerr << "Dune reported error: " << e << endl;
  	}
  	catch(exception & e){
    		cerr << "Error: " << e.what() << endl;
 	}
  	catch (...){
    		cerr << "Unknown exception thrown!" << endl;
  	}
	return 0;
}













