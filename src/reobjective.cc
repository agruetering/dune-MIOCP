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
#include <dune/MIOCP/heatdriver.hh>
#include <dune/MIOCP/hfunctions.hh>


using namespace Dune; 
using namespace std;

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


/* problem for calculating Sd in CG algorithm */
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

	CGProblem(const ControlVector& w, double dt_,int timesteps) : BaseType(w,dt_,timesteps) {}

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
template<class GV>
class YDProblem : public AdjointProblemInterface<GV>
{
protected:
	typedef AdjointProblemInterface<GV> BaseType;
	typedef BlockVector<FieldVector<double,1>> Vector;

	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::t;
	using BaseType::T;
	using BaseType::dt;
	
	vector<double> tp;
	double alpha;
	int m;

public: 
	YDProblem(vector<double> jp,double alpha_,double endTime_, int timesteps) : tp(jp), alpha(alpha_), m(tp.size()-1),BaseType(endTime_,timesteps) {}

	YDProblem(vector<double> jp,double alpha_,vector<double> dt_) : tp(jp), alpha(alpha_), m(tp.size()-1),BaseType(dt_) {}
	

	double rhs(const E& e, const X& x, int l) const
	{
		/* p(t,x)=-\alpha*c*(ud(t)-0.5)*sin(pi*x1)*sin(pi*x2)
		* with 	c= 1/(\int_\Omega \Psi*sin(pi*x1)*sin(pi*x2))=pi^4/(32+2*pi^2)
		* 	\nabla_t p(t,x)=-\alpha*c*ud'(t)*sin(pi*x1)*sin(pi*x2)
		* 	\Delta p(t,x) = \alpha*c*(ud(t)-0.5)*2*pi^2*sin(pi*x1)*sin(pi*x2)
		*
		* yd(t,x) = S(ud\Psi)+\nabla_t p(t,x)+\Delta p(t,x) 
		*/
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
		val*=2*pow(M_PI,2);
		// val -= ud'(t)
		val-=dval;
		// val*=alpha*c
		val*=alpha; 
		val*=pow(M_PI,4)/(2*(16+pow(M_PI,2)));
		// val*=sin(pi*x1)*sin(pi*x2)
		val*=sin(M_PI*global[0]);
		val*=sin(M_PI*global[1]); 

		return val;	
	}
};


int main(int argc, char** argv){

       //===============================================================
       // Make parameter structure 
       //===============================================================
    	// open ini file (contains grid information, filenames for output)
    	Dune::ParameterTree ptree;
    	Dune::ParameterTreeParser ptreeparser;
        std::stringstream file;
        file << MIOCP_SOURCE_PATH << "reobjective.ini"; 
    	ptreeparser.readINITree(file.str(),ptree);

	//===============================================================
        // make 2D grid
        //===============================================================
        const int dim=2;
        typedef Dune::UGGrid<dim> Grid;	
        Dune::FieldVector<double,dim> L;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.L" << i; 
		L[i] = ptree.get<double>(name.str(),(double)0.0);
	}
	Dune::FieldVector<double,dim> U;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.U" << i; 
		U[i] = ptree.get<double>(name.str(),(double)1.0);
	}
        std::array<unsigned int,dim> N;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.N" << i; 
		N[i] = ptree.get<unsigned int>(name.str(),(unsigned int) 100);
	} 
	std::shared_ptr<Grid> gridp = StructuredGridFactory<Grid>::createSimplexGrid(L,U,N);
	
   	typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();

        typedef Dune::PDELab::ConformingDirichletConstraints CON; // conforming dirichlet constraints
	typedef PDELab::ISTL::VectorBackend<> VBE; // ISTL backend to storage Dofs
	typedef typename GV::ctype DF;
	typedef typename PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM; 
	HeatProblem<GV> heat; // heat problem with boundary data
	typedef typename PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS; // grid function space
	typedef typename PDELab::Backend::Vector<GFS,double> vtype;
	typedef BlockVector<FieldVector<double,1>> Vector; // type of control vector

	//===============================================================
	// load instance instance 
	//=============================================================== 
	
	cout << "=== load instance" << endl;
	double T=ptree.get<double>("problem.T",(double)1.0); 
	int timesteps= ptree.get<int>("problem.Nt",(int)100);
	double h=T/timesteps;
	vector<double> tgrid(timesteps,T/timesteps);
	double alpha= ptree.get<double>("problem.alpha",(double) 0.0);
	int jumps=ptree.get<double>("jump.s",(int)10); 
	vector<double> tp(jumps+1);
	tp[0]=0;
	for(int i=1; i < jumps+1; i++){
		stringstream name; 
		name << "jump.t_" << i; 
		tp[i] = ptree.get<double>(name.str());
	}
	
	Vector ud(timesteps,0);
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
	for(int k=j; k < timesteps; k++){
		double t = (k+1)*h;
		ud[k]=0.5*pow(-1,jumps)*pow(t-tp[jumps],2)/pow(T-tp[jumps],2);
		ud[k]+=pow(-1,jumps+1)*pow(t-tp[jumps],2)*(t-T)/pow(T-tp[jumps],3);
		ud[k]+=jumps%2;
	}

	int nt=argc-2;
	vector<double> dt(nt,T/nt);
	BlockVector<FieldVector<double,1>> v(nt); 
	for(int i=0; i < nt; i++) 
		v[i]=atof(argv[i+2]);
	BlockVector<FieldVector<double,1>> u(timesteps);

	/* time points of finer temporal grid */
	vector<double> t_new(tgrid.size());
	t_new[0]=tgrid[0]; 
	for(int i=1; i < tgrid.size(); i++)
		t_new[i]=t_new[i-1]+tgrid[i];

	/* time points of coarser temporal grid */
	vector<double> t(dt.size());
	t[0]=dt[0]; 
	for(int i=1; i < dt.size(); i++)
		t[i]=t[i-1]+dt[i];

	for(int i=0; i < tgrid.size(); i++){
		// first time point of coarser temporal grid that is not less than t_new[i]
		auto start=lower_bound(t.begin(),t.end(), t_new[i],[](const double& a, double value){return a < value-1e-10;});
		int index=distance(t.begin(),start);
		// set value of u on (t_new[i-1],t_new[i]]
		if(index==0 || t_new[i-1]+1e-10>=t[index-1]) // (t_new[i-1],t_new[i]] \subseteq (t[index-1],t[index]]
			u[i]=v[index]; 
		else{ // coarser grid jumps in (t_new[i-1],t_new[i]]
			double factor=(t_new[i]-t[index-1])/dt[i]; 
			u[i]=(1-factor)*v[index]+factor*v[index-1];

		}
	}

	BlockVector<FieldVector<double,1>> w(timesteps);
	w=u;
	w-=ud;

	// calculation of S(u-ud)
	FEM fem(gv); 
	GFS gfs(gv,fem);
	vector<vtype> y; 
	CGProblem<GV> problem(w,tgrid);
	heatdriver(gfs,problem,tgrid,y);

	YDProblem<GV> desired(tp,alpha,tgrid);
	
	//===============================================================
	// recalculation of objective 
	//===============================================================
	cout<< "=== calculate objective "<< endl;
    	double bound=0; 
    	// distribution 1/2 int_\Omega \int_[0,T] (y-y_d)^2 dt dx
    	// loop over all elements of leaf grid
    	L2Deviation(gfs,tgrid,y,desired,bound);

    	// distribution alpha/2 \int_[0,T] ||u-u_d||^2 dt 
    	for (int i=0; i < u.size(); i++){
		auto diff=u[i]-heat.ud((i+1)*h); 
       		bound+=diff.two_norm2()*alpha*h*0.5; 
	}

	stringstream outputfile;
	outputfile << MIOCP_OUTPUT_PATH << "objective.txt";
	fstream stream;
	stream.precision(10);
	stream.open(outputfile.str(), std::ios::out | std::ios::app);

	stream  << "l_" << argv[1] << ":" << bound << endl;

	stream.close();

	return 0;
}

