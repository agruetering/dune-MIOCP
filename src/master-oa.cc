// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// C, C++ includes
#include<ctime> 
#include<math.h>
#include<iostream>
#include<vector>
#include <stdlib.h> 

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
// dune pdelab includes
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
#include <dune/MIOCP/Parameter.hh> 
#include <dune/MIOCP/Dmax.hh> 
#include <dune/MIOCP/outerapprox.hh> 

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

	FieldVector<double,1> Psi (const E& e, const X& x) const
  	{
		int dimDomain = x.size(); 
		auto global = e.geometry().global(x);   		
		double value = 1.5-2.0*pow(global[0]-0.5,2)-2.0*pow(global[1]-0.5,2);
		return value; 
  	}

	
};


/* problem for calculating SGu in CG algorithm */
template<class GridView>
class CGProblem : public CGProblemInterface<GridView,1>
{
protected:

	typedef CGProblemInterface<GridView,1> BaseType;

	using typename BaseType::ControlVector;
	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 
	
	using BaseType::t;

public:

	CGProblem(const ControlVector& w, double T,int timesteps) : BaseType(w,T,timesteps) {}

	CGProblem(const ControlVector& w, vector<double> dt_) : BaseType(w,dt_) {}

	FieldVector<double,1> Psi (const E& e, const X& x) const
  	{
		int dimDomain = x.size(); 
		auto global = e.geometry().global(x);   		
		double value = 1.5-2.0*pow(global[0]-0.5,2)-2.0*pow(global[1]-0.5,2);
		return value; 
  	}

	
};

/* problem to calculate desired temperature contribution yd (with boundary data) */ 
template<class GFS>
class YDProblem : public AdjointProblem<GFS>
{
protected:
	typedef   AdjointProblem<GFS> BaseType;
	typedef BlockVector<FieldVector<double,1>> Vector;

	using typename  BaseType::Y;
	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::t;
	using BaseType::T;
	using BaseType::dt;
	using BaseType::gfs;
	using BaseType::y;
	
	vector<double> tp;
	double alpha;
	int m;

public: 
	YDProblem(vector<double> jp,double alpha_,double endTime_, int timesteps, vector<Y>& y_,const GFS& gfs_) : tp(jp), alpha(alpha_), m(tp.size()-1),BaseType(endTime_,timesteps,y_,gfs_) {}

	YDProblem(vector<double> jp,double alpha_,vector<double> dt_, vector<Y>& y_,const GFS& gfs_) : tp(jp), alpha(alpha_), m(tp.size()-1),BaseType(dt_,y_,gfs_) {}
	

	double rhs(const E& e, const X& x, int l) const
	{
		/* p(t,x)=-\alpha*c*(ud(t)-0.5)*sin(pi*x1)*sin(pi*x2)
		* with 	c= 1/(\int_\Omega \Psi*sin(pi*x1)*sin(pi*x2))=pi^4/(32+2*pi^2)
		* 	\nabla_t p(t,x)=-\alpha*c*ud'(t)*sin(pi*x1)*sin(pi*x2)
		* 	\Delta p(t,x) = \alpha*c*(ud(t)-0.5)*2*pi^2*sin(pi*x1)*sin(pi*x2)
		*
		* yd(t,x) = S(ud\Psi)+\nabla_t p(t,x)+\Delta p(t,x) 
		*/
		typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF;
    		typedef typename DGF::Traits::RangeType RangeType;

		double t=accumulate(dt.begin(), dt.begin()+l, 0.0);
		auto global = e.geometry().global(x);
		double val =-0.5; 
		double dval=0; // ud'(t)
		if(t>=tp[m]){
			// ud(t)
			val+=0.5*pow(-1,m)*pow(t-tp[m],2)/pow(T-tp[m],2);
			val+=pow(-1,m+1)*pow(t-tp[m],2)*(t-T)/pow(T-tp[m],3);
			val+=m%2;

			// ud'(t)
			dval+=1.0*pow(-1,m)*(t-tp[m])/pow(T-tp[m],2);
			dval+=1.0*pow(-1,m+1)*pow(t-tp[m],2)/pow(T-tp[m],3);
			dval+=2.0*pow(-1,m+1)*(t-tp[m])*(t-T)/pow(T-tp[m],3);
		}
		else{
			int j=1; 
			while(tp[j] < t)
				j++;
			
			// ud(t)
			val+=pow(-1,j+1)*pow(t-tp[j-1],2)/pow(tp[j]-tp[j-1],2);
			val+=2*pow(-1,j)*pow(t-tp[j-1],2)*(t-tp[j])/pow(tp[j]-tp[j-1],3);
			val+=(j+1)%2;

			// ud'(t)
			dval+=2.0*pow(-1,j+1)*(t-tp[j-1])/pow(tp[j]-tp[j-1],2);
			dval+=2.0*pow(-1,j)*pow(t-tp[j-1],2)/pow(tp[j]-tp[j-1],3);
			dval+=4.0*pow(-1,j)*(t-tp[j-1])*(t-tp[j])/pow(tp[j]-tp[j-1],3);
		}
		// val = 2*pi^2(ud(t)-0.5)
		val*=2; 
		val*=pow(M_PI,2);
		// val -= ud'(t)
		val-=dval;
		// val*=alpha*c
		val*=alpha; 
		val*=pow(M_PI,4); 
		val/=(2*(16+pow(M_PI,2)));
		// val*=sin(pi*x1)*sin(pi*x2)
		val*=sin(M_PI*global[0]);
		val*=sin(M_PI*global[1]); 

		int i=round(t/0.005);
		DGF ydgf(gfs,y[i]);
		RangeType value; 
    		ydgf.evaluate(e,x,value); 
		value+=val;
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
        file << MIOCP_SOURCE_PATH << "master-oa.ini"; 
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
	
   	typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();
	
	//===============================================================
	// create instance 
	//===============================================================
	
	cout << "=== create instance" << endl;
        typedef PDELab::ConformingDirichletConstraints CON; // conforming dirichlet constraints
	typedef PDELab::ISTL::VectorBackend<> VBE; // ISTL backend to storage Dofs
	typedef typename GV::ctype DF; // coordinate type of spatial grid
	typedef typename PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM; // finite element map
	HeatProblem<GV> heat; // heat problem with boundary data
	typedef typename PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS; // grid function space
	typedef typename PDELab::Backend::Vector<GFS,double> vtype; // type of Dof vector
	typedef BlockVector<FieldVector<double,1>> ControlVector; // type of control vector

	double T=ptree.get<double>("problem.T",(double)1.0); 
	int timesteps= ptree.get<int>("problem.Nt",(int)100);
	vector<double> dt(timesteps,T/timesteps);
	double alpha= ptree.get<double>("problem.alpha",(double) 0.0);
	
	
        int jumps=11; 
	// random generated jump points
	vector<int> jpoints(jumps); 
	bool stop;
	// set up seed
	unsigned int seed=time(0);
        srand(seed);
	
	// generate random jump points
	for (int i = 0; i < jumps; i++){
        	bool same;
        	do
       	 	{
           	 	same = false;
            		jpoints[i] = rand() % (timesteps-1);
            		// Check if the newly generated number is a duplicate
            		for (int j = 0; j < i;j++){
               			if (jpoints[i] == jpoints[j]){
                    			same = true;
                    			break;
                		}
            		}
        	} while (same);
    	}
	sort(jpoints.begin(),jpoints.end()); 

	// calcuate ud
	ControlVector ud(400,0);
	double h=T/400;
	vector<double> tp(jumps+1);
	tp[0]=0;
	for(int i=1; i < jumps+1; i++)
		tp[i]=(jpoints[i-1]+1)*(T/timesteps);
		
	// output jump points
        cout << "Jump times" << endl;
        for(int i=1; i < jumps+1; i++)
		cout << "t_" << i << ":" << tp[i]<< endl;

	// cubic spline interpolation between in [t_{i-1},t_i] 
	// with ud(0)=0, ud(T)=0.5 and 
	// ud(t_i)=1 if i odd, else ud(t_i)=0
	int j=0;
	for(int i=1; i < jumps+1; i++){
		while( (j+1)*h < tp[i]){
			double t=(j+1)*h;
			ud[j]+=pow(-1,i+1)*pow(t-tp[i-1],2)/pow(tp[i]-tp[i-1],2);
			ud[j]+=2*pow(-1,i)*pow(t-tp[i-1],2)*(t-tp[i])/pow(tp[i]-tp[i-1],3);
			ud[j]+=(i+1)%2;
			j++;
		}
	}
	for(int k=j; k < 400; k++){
		double t = (k+1)*h;
		ud[k]=0.5*pow(-1,jumps)*pow(t-tp[jumps],2)/pow(T-tp[jumps],2);
		ud[k]+=pow(-1,jumps+1)*pow(t-tp[jumps],2)*(t-T)/pow(T-tp[jumps],3);
		ud[k]+=jumps%2;
	}
	
	// calculation of desired temperature yd(t,x)
	FEM fem(gv); 
	GFS gfs(gv,fem);
	vector<vtype> yd; 
	vector<double> tgrid(400,h);
	CGProblem<GV> problem(ud,tgrid);
	heatdriver(gfs,problem,tgrid,yd);
	YDProblem<GFS> AdjDesired(tp,alpha,dt,yd,gfs);
	
		
        typedef Param<GV,CON,CGProblem<GV>,1,1> PStruct; 
        PStruct pstruct(gv,     /* grid */
			heat, 	/* heat problem */
			AdjDesired, /* desired temperature problem */
			alpha, /* Tikhonov parameter */
			1e-5, /* parameter for active cutting planes */
			T,  /* end time */
			timesteps, /* number of time intervals */
                        50,     /* maximum iterations active set method */
                        50,    	/* maximum linear solver iterations */ 
			1e-8,  /* absolute tolerance*/
			1e-8,  /* reduction factor */
                        6,   	/* maximum iterations outer approximation */
			86400); /* time limit */

	//===============================================================
       	// Initialization
       	//===============================================================
	cout << "\n---------------------- Initialization ----------------------\n" << endl; 
	
	// initialize u^0
        ControlVector u(timesteps);
           
	// allowed switchings 
	int smax=2; 		

	// store result 
	OutputData data;
	clock_t cstart=clock();

	
	vector<vtype> Sf(timesteps+1,vtype(gfs));
	
	cout << "=== calculation of desired temperature contribution" << endl; 
	ControlVector Gstar_eta(timesteps);
	{
		// calculation eta=S*yd
    		vector<vtype> eta; 
		adjointdriver(gfs,AdjDesired,eta);
        	
		// contibution of G*S*yd=\int Psi(x)eta(t,x) dx 
		Gstar(gfs,heat,eta,Gstar_eta);
	}
	Gstar_eta*=-1; // encodes G*S*(Sf-yd)=-G*S*yd
	data.time_preprocess=(clock()-cstart)/CLOCKS_PER_SEC;
 
	// write data to file
	stringstream outputfile;
	outputfile << MIOCP_OUTPUT_PATH << ptree.get("output.dir","output") << "NX" << N[0] <<"x"<< N[1] << "_step" << timesteps << ".txt";
	fstream stream;
	stream.precision(10);
	stream.open(outputfile.str(), ios::out | ios::trunc); 
	
	stream << "\n[Jump times]" << endl;
	for(int i=1; i < jumps+1; i++)
		stream << "t_" << i << ":" << tp[i]<< endl;
		
	stream << "\n[sgrid]" << endl; 
	stream << "LX=" << L[0] << endl; 
	stream << "LY=" << L[1] << endl;
	stream << "UX=" << U[0] << endl; 
	stream << "UY=" << U[1] << endl;
	stream << "NX=" << N[0] << endl; 
	stream << "NY=" << N[1] << endl;

	stream << "\n[tgrid]" << endl;
	stream << "nodes=" << pstruct.dt.size() << endl;
	for(int i=0; i < pstruct.dt.size(); i++)
		stream << "dt_" << i << "=" << pstruct.dt[i] << endl; 

	stream << "\n[Parameter]" << endl;
	stream << "Seed=" << seed << endl;
	stream << "degree=" << pstruct.degree << endl;  
	stream << "Tikhonov=" << pstruct.alpha << endl; 
	stream << "Rho=" << pstruct.rho << endl;
	stream << "Iter_ActiveSet=" << pstruct.iter_activeset << endl;
	stream << "Iter_Lin=" << pstruct.iter_lin << endl;
	stream << "Tol=" << pstruct.tol << endl;
	stream << "Reduction=" << pstruct.reduction << endl;
	stream << "Iter_Outer=" << pstruct.iter_outer << endl;
	stream << "TimeLimit=" << pstruct.timeouter << endl;
	stream << "smax=" << smax << endl;
	
	stream.close();
	
	vector<vtype> y;
	try{
		Dmax D(smax);
		// outer approximation of relaxed problem
       		OuterApprox<PStruct> oasolver(pstruct, Sf, D, data);
		y = oasolver.apply(Gstar_eta,u,ptree.get("output.dir","output")); 
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

	stream.open(outputfile.str(), ios::out | ios::app); 
	
	// write output data to file
	stream << "\n[Output]" << endl; 
	stream << "Iter_ActiveSet=" << data.iter_activeset << endl;
	stream << "Iter_Lin=" << data.iter_lin << endl;
	stream << "Iter_OuterApprox="<< data.iter_cut << endl;
	stream << "Time_Preprocessing=" << data.time_preprocess << endl;
	stream << "\n LowerBounds" << endl;
	for(int i=0; i < data.lower_bound.size(); i++)
		stream <<  "l_" << i << ":" << data.lower_bound[i] << endl;

	stream << "\nTimeOuter" << endl;
	for(int i=0; i < data.time_outer_solve.size(); i++)
		stream << "o_" << i << ":" << data.time_outer_solve[i] << endl;
	stream << "o_total=" << accumulate(data.time_outer_solve.begin(),data.time_outer_solve.end(),(double)0.0) << endl;

	stream << "\nTimeLinSolve" << endl;
	for(int i=0; i < data.time_lin_solve.size(); i++)
		stream << "r_" << i << ":" << data.time_lin_solve[i] << endl;
	stream << "r_total=" << accumulate(data.time_lin_solve.begin(),data.time_lin_solve.end(),(double)0.0)<< endl;

	stream << "\nTimeLinAssemble: " << endl;
	for(int i=0; i < data.time_lin_solve.size(); i++)
		stream << "a_" << i << ":" << data.time_lin_assemble[i] << endl;
	stream << "a_total=" << accumulate(data.time_lin_assemble.begin(),data.time_lin_assemble.end(),(double)0.0) << endl;

	stream << "\nTimeSeparation" << endl;
	for(int i=0; i < data.time_outer_separation.size(); i++)
		stream << "s_" << i << ":" << data.time_outer_separation[i] << endl;
	stream << "s_total=" << accumulate(data.time_outer_separation.begin(),data.time_outer_separation.end(),(double)0.0)<<endl;

	stream << "\n[End_control]" << endl; 
	for(int i=0; i < u.size(); i++)
		stream << "u_"<< i << ":" << u[i] << endl;

	stream << "\nOptimization stopped with status " << data.status << endl;
	
	stream.close();
	
	return 0;
}
