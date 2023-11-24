#ifndef GUROBI_HH
#define GUROBI_HH

// dune includes
#include <dune/istl/matrixindexset.hh>

#include <gurobi_c++.h>

#include <dune/MIOCP/probleminterface.hh>
#include <dune/MIOCP/callback.hh>
#include <dune/MIOCP/heatdriver.hh>
#include <dune/MIOCP/hfunctions.hh>

using namespace Dune; 
using namespace std; 


/* discretization of the optimization problem
 * min J(y,u)=1/2|y-y_d|_{L^2(Q)}^2 
 *    s.t.  \partial_t y - \Delta y = Gu + f = \sum_j u_j(t)*Psi_j(x) + f	in Q=\Omega x (0,T)
 *                       		   y = 0   		on \partial \Omega x (0,T)
 *                        		y(0) = y0		
 */
template<class GFS, int k> 
class PDEElements
{
public: 
	typedef BCRSMatrix<FieldMatrix<double,1,1>>Matrix; 

	typedef typename GFS::Traits::GridView GV; // grid view
	typedef typename GFS::Traits::FiniteElementMap FEM; // finite element map 
	typedef typename GV::ctype ctype; // coordinate type of the grid
	typedef typename GV::IndexSet LeafIndexSet; // leaf index set
  	typedef typename GV::IntersectionIterator IntersectionIterator; // intersection iterator

	enum {n=k}; // number of switches
	enum {dim=GV::dimension}; // dimension of spatial ansatz

  	typedef typename GV::template Codim<0>::Geometry::JacobianInverseTransposed JacobianInverseTransposed;
	typedef typename PDELab::Backend::Vector<GFS,double> Y; // type of Dofs vectors
	typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
	typedef PDELab::LFSIndexCache<LFS> LFSCache; // lfs cache
	typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
	typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache
	typedef typename LocalBasis::Traits::JacobianType JacobianType; // jacobian type
	typedef typename LocalBasis::Traits::RangeType RangeType; // range type

	typedef HeatProblemInterface<GV,n> HeatProblem; 
	typedef AdjointProblemInterface<GV> AdjProblem;
	
public: 
	Matrix massMatrix; // mass matrix
	Matrix mMatrix;	// mass matrix with unit row for dirichlet nodes
	Matrix stiffMatrix; // stiffness matrix
	vector< BlockVector<double>> Myd; // mass matrix multiplied with desired temperature Dof vector
	vector< BlockVector<double>> yd; // desired temperature Dof vector
	BlockVector<FieldVector<double,n>> Mpsi; // mass matrix multiplied with Dof vector of form functions
	vector< BlockVector<double>> Mf; // mass matrix multiplied with Dof vector of function f
	BlockVector<double> initial; // Dof vector of initial values
	MatrixIndexSet adjacencyPattern; // pattern to set up mass/stiffness matrix
private: 
	const GV gv; // spatial grid
	FEM fem; // finite element method
	GFS gfs; // grid function space
	const vector<double> dt; // temporal grid
	
public: 

	PDEElements(const GV gv_,const vector<double> dt_,HeatProblem& problem,AdjProblem& desired) :  gv(gv_), fem(FEM(gv)), gfs(GFS(gv,fem)), dt(dt_){
		getAdjacencyPattern();
		assembleStiff();
		assembleMass();
		assembleRhs(problem);
		assembleObj(desired);
	}
 
	// pattern to set up matrices
	void getAdjacencyPattern();
  
	// assembling methods
  	void assembleStiff();
  	void assembleRhs(HeatProblem& heat);
	void assembleObj(AdjProblem& des);
        void assembleMass(); 

	GFS gridFunctionSpace(){return gfs;}
	double deltaT(int i) { return dt[i];}
	int steps(){return dt.size();}
};


/* adjacency pattern of grid vertices */
template<class GFS, int k> 
void PDEElements<GFS,k>::getAdjacencyPattern(){

	// number of vertices of the leaf grid
	const int N=gv.size(dim); 

	adjacencyPattern.resize(N,N); 
	// global indices of entities in the leaf grid
	const LeafIndexSet& set = gv.indexSet(); 
	
	// loop over all grid element
	for(auto&&e : elements(gv)){
		// get reference element and number of vertices
		auto geo=e.geometry(); 
		auto ref=referenceElement(geo); 
		int vertexsize = ref.size(dim); 
		
		// store all pairs of vertices in adjacencyPattern
		for (int i=0; i<vertexsize;i++){
			auto indexi=set.subIndex(e,i,dim); // global index of i-th vertex of element e
			for (int j=0; j<vertexsize;j++){
				auto indexj=set.subIndex(e,j,dim); // global index of j-th vertex of element e
				adjacencyPattern.add(indexi,indexj); 
			}
		}
	
	}
}


/* assembler stiffness matrix */
template<class GFS, int k> 
void PDEElements<GFS, k>::assembleStiff(){

	// number of vertices of the leaf grid
	const int N=gv.size(dim); 
	
	// global indices of entities in the leaf grid
	const LeafIndexSet& set = gv.indexSet();
	
	// initialize stiffMatrix 
	stiffMatrix.setSize(N,N); 
	stiffMatrix.setBuildMode(Matrix::random); 
	getAdjacencyPattern();
	adjacencyPattern.exportIdx(stiffMatrix); 
	stiffMatrix=0.0; 

	// local function space and basis
	LFS lfs(gfs);
 	LFSCache lfs_cache(lfs);
 	BasisCache cache;
	
	// loop over all grid elements
	for(auto&&e : elements(gv)){
		// get reference element and number of vertices
		auto geo = e.geometry(); 
		auto ref = referenceElement(geo); 
		int vertexsize = ref.size(dim);
		
		// local function space on current element
		lfs.bind(e); 
       		lfs_cache.update();
		
		// choose appropiate quadrature rule
        	const int order =  2*lfs.finiteElement().localBasis().order();
       		auto rule = PDELab::quadratureRule(geo,order);

		// loop over quadrature points
		for (const auto& ip : rule)
		{
			
			// compute the jacobian inverse transposed to transform the gradients
			const auto jacobian = geo.jacobianInverseTransposed(ip.position());
			
			// get the weight at the current quadrature point
			ctype weight=ip.weight(); 
			
			// compute Jacobian determinant for the transformation formula
			ctype detjac = geo.integrationElement(ip.position());

			// compute jacobian at current quadrature point
			vector< JacobianType > phihat = cache.evaluateJacobian(ip.position(),
                	lfs.finiteElement().localBasis());
			 
			// compute transformed gradients 
			for (int i=0; i<vertexsize; i++){
			
				FieldVector<ctype,dim> grad1; 
				jacobian.mv(phihat[i][0],grad1);
				
				for(int j=0; j<vertexsize;j++){
				FieldVector<ctype,dim> grad2; 
				jacobian.mv(phihat[j][0],grad2);
				
				// gain global inidices of vertices i and j and update associated matrix entry
				stiffMatrix[set.subIndex(e,i,dim)][set.subIndex(e,j,dim)]
					+=(grad1*grad2)*weight*detjac;
				}
			}
		}
	}
	
	// Dirichlet boundary conditions:
  	// replace lines in stiffMatrix related to Dirichlet vertices by trivial lines
	for(auto&&e : elements(gv)){
		auto geo = e.geometry();
      		auto ref = referenceElement(geo);
	
		const IntersectionIterator isend=gv.iend(e); 
		for (IntersectionIterator is = gv.ibegin(e) ; is != isend ; ++is){
  
      			// check whether current intersection is on the boundary
      			if ( is->boundary() ){
			
        	   		// traverse all vertices the intersection consists of
        	   		for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++){
          				// and replace the associated line of stiffMatrix with a trivial one
          				int indexi = set.subIndex(e,ref.subEntity(is->indexInInside(),1,i,dim),dim);

          				stiffMatrix[indexi] = 0.0;           
          				stiffMatrix[indexi][indexi] = 1.0;
        	   		}
      			}	
		}
	}
}


/* assembler mass matrix */
template<class GFS, int k> 
void PDEElements<GFS,k>::assembleMass(){

	// number of vertices of the leaf grid
	const int N=gv.size(dim); 
	const LeafIndexSet& set = gv.indexSet();
	
	// initialize massMatrix 
	massMatrix.setSize(N,N); 
	massMatrix.setBuildMode(Matrix::random); 
	getAdjacencyPattern();
	adjacencyPattern.exportIdx(massMatrix); 
	massMatrix=0.0; 
	
	// local function space and basis
	LFS lfs(gfs);
 	LFSCache lfs_cache(lfs);
 	BasisCache cache;
	
	// loop over all grid elements
	for(auto&&e : elements(gv)){
		// get reference element and number of vertices
		auto geo = e.geometry(); 
		auto ref = referenceElement(geo); 
		int vertexsize = ref.size(dim);
		
		// local function space on current element
		lfs.bind(e); 
       		lfs_cache.update();
		
		// choose appropiate quadrature rule
        	const int order = 2*lfs.finiteElement().localBasis().order();
       		auto rule = PDELab::quadratureRule(geo,order);

		// loop over quadrature points
		for (const auto& ip : rule){
			// compute the jacobian inverse transposed to transform the gradients
			const auto jacobian = geo.jacobianInverseTransposed(ip.position());
			
			// get the weight at the current quadrature point
			ctype weight=ip.weight(); 
			
			// compute Jacobian determinant for the transformation formula
			ctype detjac = geo.integrationElement(ip.position());

			// compute jacobian at current quadrature point
			vector<RangeType > phihat = cache.evaluateFunction(ip.position(),
                	lfs.finiteElement().localBasis());
			 
			for (int i=0; i<vertexsize; i++){
				
				for(int j=0; j<vertexsize;j++){
				
				// gain global inidices of vertices i and j and update associated matrix entry
				massMatrix[set.subIndex(e,i,dim)][set.subIndex(e,j,dim)]				
				+=phihat[i]*phihat[j]*weight*detjac;
				}
			}
		}
	}

	// mass matrix with trivial lines for Dirichlet vertices
	mMatrix.setSize(N,N);
	mMatrix.setBuildMode(Matrix::random); 
	getAdjacencyPattern();
	adjacencyPattern.exportIdx(mMatrix); 
	mMatrix=massMatrix; 

	for(auto&&e : elements(gv)){
		auto geo = e.geometry();
      		auto ref = referenceElement(geo);
	
		const IntersectionIterator isend=gv.iend(e); 
		for (IntersectionIterator is = gv.ibegin(e) ; is != isend ; ++is){
  

      			// check whether current intersection is on the boundary
      			if ( is->boundary() ){
			
        	   		// traverse all vertices the intersection consists of
        	   		for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++){
          				// and replace the associated line of massMatrix with a trivial one
          				int indexi = set.subIndex(e,ref.subEntity(is->indexInInside(),1,i,dim),dim);
          				mMatrix[indexi] = 0.0;           
        	   		}
      			}	
		}
	}
}	


/* assembler right hand side */
template<class GFS, int k> 
void PDEElements<GFS,k>::assembleRhs(HeatProblem& heat){
	// number of vertices of the leaf grid
	const int N=gv.size(dim);  
	
	const LeafIndexSet& set = gv.indexSet();
	 
	// calculation of mass matrix multiplied with Dof vector of form functions
	Mpsi.resize(N);
	for(int j=0; j < n; j++){
		BlockVector<double> Psi(N);
		// lambda function for i-th form function function
 		auto psilambda = [&](const auto& e, const auto& x)
    		{
			return heat.Psi(e,x)[j];
		};
		auto g = PDELab::makeGridFunctionFromCallable(gv,psilambda);


		// fill the Dof vector
		Y v(gfs); 
		PDELab::interpolate(g,gfs,v);

		// save Dof Vector 
		for(int i=0; i < N; i++)
			Psi[i]=v.block(i);

		// calculate rhs 
		BlockVector<double> h(N);
		massMatrix.mv(Psi,h); 
		
		// save result 
		for(int i=0; i < N; i++) 
			Mpsi[i][j]=h[i]; 
	}
	// calculation of mass matrix multiplied with Dof vector of f
	Mf.clear();
	for(int i=0; i< this->steps()+1; i++){
		double time=accumulate(dt.begin(),dt.begin()+i,0.0);
		heat.setTime(time);
	        // lambda function for function f at time t
 		auto flambda = [&](const auto& e, const auto& x)
    	    	{
		        return heat.q(e,x);
		};
		// make grid function for f at time t
 	 	auto f = PDELab::makeGridFunctionFromCallable(gv,flambda);

		// lagrange interpolation of function f at time t
		Y v(gfs);
		PDELab::interpolate(f,gfs,v);

		// save Dof vector as BlockVector
		BasisCache cache; 
		BlockVector<double> F(N);
		for(int l=0; l < N; l++) 
			F[l]=v.block(l);
		// calculate Mf for current time
		BlockVector<double> h(N);
		massMatrix.mv(F,h); 
		// save value
		Mf.push_back(h);	
	}
	
	// Dirichlet boundary conditions:
  	// replace entry in right hand side related to Dirichlet vertices by 0
	for(auto&&e : elements(gv)){
		auto geo = e.geometry();
      		auto ref = referenceElement(geo);
		const IntersectionIterator isend=gv.iend(e); 
		for (IntersectionIterator is = gv.ibegin(e) ; is != isend ; ++is){
  
      			// check whether current intersection is on the boundary
      			if ( is->boundary() ){
			
        	   		// traverse all vertices the intersection consists of
        	   		for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++){
          				// and replace the associated entry of right hand side by 0
          				int indexi = set.subIndex(e,ref.subEntity(is->indexInInside(),1,i,dim),dim);
          				Mpsi[indexi] = 0.0;
					for(int j=0; j < this->steps()+1; j++) 
						Mf[j][indexi]=0.0;           
        	   		}
      			}	
		}
	}


	// lambda function for initial function
 	auto inilambda = [&](const auto& e, const auto& x)
    	{
		return heat.ini(e,x);
	};

	// make grid function for initial values
 	auto ini= PDELab::makeGridFunctionFromCallable(gv,inilambda);

	// fill the Dof vector of initial values
	Y v(gfs); 
	PDELab::interpolate(ini,gfs,v);
	initial.resize(N); 
	for(int i=0; i < N; i++)
		initial[i]=v.block(i);

}


/* assembler objective (contribution of yd) */
template<class GFS, int k> 
void PDEElements<GFS,k>::assembleObj(AdjProblem& des){

	// number of vertices of the leaf grid
	const int N=gv.size(dim);  
	
	const LeafIndexSet& set = gv.indexSet();
	
	Myd.clear();
	yd.clear();
	des.setdeltaT(dt);
	for(int i=0; i < this->steps()+1; i++){
	        // lambda function for desired temperature at time t
 		auto desiredlambda = [&](const auto& e, const auto& x)
    	    	{
		        return des.rhs(e,x,i);
		};
		// make grid function for desired temperature at time t
 	 	auto desired_temp = PDELab::makeGridFunctionFromCallable(gv,desiredlambda);

		// lagrange interpolation of desired temperature at time t
		Y v(gfs);
		PDELab::interpolate(desired_temp,gfs,v);

		// save Dof vector as BlockVector
		BasisCache cache; 
		BlockVector<double> d(N);
		for(int l=0; l < N; l++) 
			d[l]=v.block(l);
		yd.push_back(d);
		// calculate Myd for current time
		BlockVector<double> h(N);
		massMatrix.mv(d,h); 
		// save value
		Myd.push_back(h);
		
	}
	
}




/* Gurobi solver to solve discretization of 
 min J(y,u)=1/2|y-y_d|_{L^2(Q)}^2 
 *    s.t.  \partial_t y - \Delta y = Gu + f = u(t)*Psi(x) + f	in Q=\Omega x (0,T)
 *                       		   y = 0 		on \partial \Omega x (0,T)
 *                        		  y(0) = y0 		
 * 
 * with switching constraint |u|_{BV(0,T)} \leq smax
 */ 
template<class GridView,class Constraints, class Problem, int l>
class Gurobi{
public:

    	enum {dim=GridView::dimension}; // dimension of spatial grid
	enum {degree = l}; // degree of spatial ansatz and test functions
	typedef GridView GV; // spatial grid view 
	typedef Constraints CON; // constraints (e.g. conforming dirichlet constraints)

	typedef Problem CGProblem; 

	typedef typename GV::ctype DF; // coordinate type of the spatial grid
	typedef typename conditional<dim>=2,PDELab::PkLocalFiniteElementMap<GV,DF,double,degree>,PDELab::QkLocalFiniteElementMap<GV,DF,double,degree>>::type FEM; // finite element map (for simplices in case dim>=2)
  	typedef PDELab::ISTL::VectorBackend<> VBE; // ISTL backend to storage Dofs
  	typedef typename PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS; 	// grid function space for spatial test and ansatz functions (Galerkin method

	typedef typename PDELab::Backend::Vector<GFS,double> Y;


private:  	
	PDEElements<GFS,1> pde; // data of discretized problem
	HeatProblemInterface<GV,1>& heat; // heat problem with box constraints
	AdjointProblemInterface<GV>& desired; // problem class for yd∫

	/* needed to recalculate exact objective on finer discretization */
	vector<Y> Sf; // right hand side contribution Sf on finer discretization

	const GV gv; // finer spatial grid
	const FEM fem; // finer finite element method
	const GFS gfs; // finer grid function space
	const vector<double> dt;
public: 

	Gurobi(const GV& sgrid,const vector<double>& tgrid, HeatProblemInterface<GV,1>& heat_, AdjointProblemInterface<GV>& desired_, const GV gv_, const vector<double> dt_): 
		heat(heat_),
		desired(desired_),
		pde(PDEElements<GFS,1>(sgrid,tgrid,heat_,desired_)),
		gv(gv_),
		fem(FEM(gv)),
		gfs(GFS(gv,fem)), 
		dt(dt_)
	{
		if(!heat.isHomogeneous())
			heatdriver(gfs,heat, dt, Sf);
		else 
			Sf=vector<Y>(dt.size()+1,Y(gfs));
	}

	Gurobi(const GV& sgrid,const vector<double>& tgrid, HeatProblemInterface<GV,1>& heat_, AdjointProblemInterface<GV>& desired_, const GV gv_, const int Nt): 
		heat(heat_),
		desired(desired_),
		pde(PDEElements<GFS,1>(sgrid,tgrid,heat_,desired_)),
		gv(gv_),
		fem(FEM(gv)),
		gfs(GFS(gv,fem)), 
		dt(vector(Nt,accumulate(tgrid.begin(),tgrid.end(),0.0)/Nt))
	{	
		if(!heat.isHomogeneous())
			heatdriver(gfs,heat, dt, Sf);
		else 
			Sf=vector<Y>(dt.size()+1,Y(gfs));
	}
		

	double objective(const double*);
	void optimize(int mod, string logfile, string resfile, int smax=2, bool heurisic=false,string dir="");

};

/* calculation of objective value on finer grid 
 */
template<class GridView,class Constraints, class Problem, int l>
double Gurobi<GridView,Constraints,Problem,l>::objective(const double* v)
{
	
	BlockVector<FieldVector<double,1>> u(dt.size());
	/* time points of finer temporal grid */
	vector<double> t_new(dt.size());
	t_new[0]=dt[0]; 
	for(int i=1; i < dt.size(); i++)
		t_new[i]=t_new[i-1]+dt[i];

	/* time points of coarser temporal grid */
	vector<double> t(pde.steps());
	t[0]=pde.deltaT(0); 
	for(int i=1; i < pde.steps(); i++)
		t[i]=t[i-1]+pde.deltaT(i);

	for(int i=0; i < dt.size(); i++){
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

	// calculate state on finer grid  
	vector<Y> y; 
	CGProblem problem(u,dt); 
	heatdriver(gfs, problem, dt, y); 
	
	double bound=0; 
    	// distribution 1/2 int_\Omega \int_[0,T] (y+Sf-y_d)^2 dt dx
    	// loop over all elements of leaf grid
	for(int i=0; i < dt.size()+1; i++)
		y[i]+=Sf[i]; 
	L2Deviation(gfs,dt,y,desired,bound);

	return bound;
}


/* Call GUROBI Optimizer 
 */
template<class GridView,class Constraints, class Problem, int l>
void Gurobi<GridView,Constraints,Problem,l>::optimize(int mod, 
string logfile, // logfile name for gurobi
string resfile, // result file name 
int smax, // number of allowed switchings
bool heuristic, // calculation of heuristic solution (if desired)
string  dir) // folder name to store heuristic solutions
{ 
	/* mod overview 
	 * 0: integer problem with additional variables for tv norm
	 * 1: naive linear relaxation of tv norm
	 * 2: tailored convexification
	 * 3: no tv constraint
	 */

	GRBEnv* env=0; 
	GRBVar** y=0; 
	GRBVar* u =0;
	GRBVar* z=0;
	env=new GRBEnv(); 
	GRBModel model=GRBModel(*env); 

	const int N=pde.gridFunctionSpace().gridView().size(dim);  
	const int timesteps=pde.steps();

	cout << "=== Set up variables u,y " << endl;
	// initialize state variables 
	// y[i][j] represents state at time i*deltaT and vertex j of spatial grid
	y=new GRBVar*[timesteps+1];
	for(int i=0; i < timesteps+1; i++){
		y[i]=model.addVars(N,GRB_CONTINUOUS); 
		for(int j=0; j < N; j++){
			stringstream yname; 
			yname << "y_" << i << "_" << j; 
			y[i][j].set(GRB_StringAttr_VarName, yname.str()); 
			y[i][j].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
		}
	}

	// initialize control variable
	// u[i] at time (i+1)*deltaT 
	u=model.addVars(timesteps,GRB_CONTINUOUS); 
	double time=0;
	for(int i=0; i < timesteps; i++){
		stringstream uname; 
		uname << "u_" << i; 
		u[i].set(GRB_StringAttr_VarName, uname.str());
		time+=pde.deltaT(i);
		double lb=heat.ua(time)[0];
		double ub=heat.ub(time)[0];
		u[i].set(GRB_DoubleAttr_LB, lb);
		u[i].set(GRB_DoubleAttr_UB, ub);
		if(mod==0) 
			u[i].set(GRB_CharAttr_VType, GRB_BINARY);				
	}

	cout << "=== Set up objective " << endl;
	// set up objective function terms 
	GRBQuadExpr obj(0.0);
	for (int i=0; i < timesteps; i++){ // loop over time points
		// objective terms in y
		// loop over spatial grid
		for(int j=0; j< N; j++){
			for(int k=0; k < N; k++){
				if(pde.massMatrix.exists(j,k)){ // quadratic terms in y	
					obj+=pde.deltaT(i)*(1.0/3.0)*y[i][j]*pde.massMatrix[j][k]*y[i][k]; 
					obj+=pde.deltaT(i)*(1.0/3.0)*y[i][j]*pde.massMatrix[j][k]*y[i+1][k]; 
					obj+=pde.deltaT(i)*(1.0/3.0)*y[i+1][j]*pde.massMatrix[j][k]*y[i+1][k];
				}
			}

			// linear terms in y
			obj-=pde.deltaT(i)*(2.0/3.0)*y[i][j]*pde.Myd[i][j]; 
			obj-=pde.deltaT(i)*(2.0/3.0)*y[i+1][j]*pde.Myd[i+1][j];
			obj-=pde.deltaT(i)*(1.0/3.0)*y[i][j]*pde.Myd[i+1][j]; 
			obj-=pde.deltaT(i)*(1.0/3.0)*y[i+1][j]*pde.Myd[i][j]; 

			// constant terms 
			obj+=pde.deltaT(i)*(1.0/3.0)*pde.yd[i][j]*pde.Myd[i][j]; 
			obj+=pde.deltaT(i)*(1.0/3.0)*pde.yd[i][j]*pde.Myd[i+1][j];				
			obj+=pde.deltaT(i)*(1.0/3.0)*pde.yd[i+1][j]*pde.Myd[i+1][j];
		}
	}
	obj*=0.5; 
	model.setObjective(obj,GRB_MINIMIZE);

	cout << "=== Set up PDE constraints " << endl;
	// set up pde constraints
	for(int i=0; i < timesteps; i++){
		for (int j=0; j < N; j++){
			 GRBLinExpr lhs = 0;
			 for (int k=0; k < N; k++){
				if(pde.mMatrix.exists(j,k)){
					lhs+=pde.mMatrix[j][k]*(y[i+1][k]-y[i][k]);
					lhs+=pde.deltaT(i)*pde.stiffMatrix[j][k]*y[i+1][k];
				}
			}
			lhs-=pde.deltaT(i)*pde.Mpsi[j][0]*u[i];
			model.addConstr(lhs,GRB_EQUAL,pde.deltaT(i)*pde.Mf[i+1][j]);
		}
	}
	for(int j=0; j < N; j++)
		model.addConstr(y[0][j],GRB_EQUAL,pde.initial[j]);


	if( mod==0 || mod==1){
		cout << "=== Set up z variables + constraints " << endl;
		// additional variables for tv norm
		z=model.addVars(timesteps-1,GRB_CONTINUOUS);
		for(int i=0; i < timesteps-1; i++){
			ostringstream zname; 
			zname << "z_" << i; 
			z[i].set(GRB_StringAttr_VarName, zname.str());
			// z_i=abs(u_{i+1}-u_i)
			model.addConstr(-z[i]+u[i+1]-u[i],GRB_LESS_EQUAL,0.0);
			model.addConstr(-z[i]-u[i+1]+u[i],GRB_LESS_EQUAL,0.0);
		}
		// tv norm constraint u_0+\sum z_i \leq smax
		GRBLinExpr lhs = 0;
		for(int i=0; i < timesteps-1; i++)
			lhs+=z[i]; 
		lhs+=u[0];
		model.addConstr(lhs,GRB_LESS_EQUAL,smax);
	}
	if(mod==2)
		model.set(GRB_IntParam_NumericFocus,1);

	model.set(GRB_StringParam_LogFile, logfile);
	model.set(GRB_IntParam_Method,1);
	model.set(GRB_IntParam_Threads,1);
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);

  	if ((status == GRB_INF_OR_UNBD) || (status == GRB_INFEASIBLE) || (status == GRB_UNBOUNDED) )
    		cout << "The model cannot be solved " << "because it is infeasible or unbounded" << endl;
	else if(status == GRB_TIME_LIMIT){
    		cout << "Optimization was stopped because time limit was reached " << endl;
		fstream result;
		result.precision(10);
		result.open(resfile, ios::out | ios::app);
		double* v=new double[timesteps];
		result << "\n[Current_control]" << endl; 
		for(int i=0; i < timesteps; i++){
			v[i]=u[i].get(GRB_DoubleAttr_X);
			result << u[i].get(GRB_StringAttr_VarName) << ":" << v[i] << endl;
		}
		result << "[Current_objective]" << endl; 
		double obj=objective(v);
		result << obj << endl;	
		result.close();	
		delete[] v;
	}
	else if((status == GRB_OPTIMAL) && (mod==2)){ 
		int count =0; // counts number of successive iterations with less than 0.1% relative change
		int cut=0; // total number of added cuts
		double time=0; // summed up runtime
		double* v = new double[timesteps]; // old solution
		double obj; // objective in current iteration
		double oldobj; //objective in previous iteration 
		
		obj=model.get(GRB_DoubleAttr_ObjVal);
		fstream result;
		result.precision(10);
		result.open(resfile, ios::out | ios::app); 
 		result << "\n[Gurobi_Objective]" << endl;
		while(model.get(GRB_IntAttr_Status)==GRB_OPTIMAL && count < 3){
			if(heuristic){
				cout << "\n=== Calculation of heuristic solution" << endl;
				double* w=new double[timesteps]; 
				for(int i=0; i < timesteps; i++)
					w[i]=u[i].get(GRB_DoubleAttr_X);
				callheuristic(timesteps,w,smax); 
				double hval=objective(w);
				stringstream hfile; 
				hfile << dir << "hsol_"<<cut << ".txt";
				fstream file; 
				file.precision(10);
				file.open(hfile.str(),  ios::out | ios::trunc); 
				file << "[Control]" << endl;
				for(int i=0; i < timesteps; i++)
					file << "w_" << i << ":" << w[i] << endl;
				file << "[Objective]" << endl;
					file << hval << endl;
				file.close();
				delete[] w;
			}
			result << obj << endl;
			time+=model.get(GRB_DoubleAttr_Runtime);
			// prove feasibility
			int flag = callback(&model,u,timesteps,smax);
			if(flag ==0){  // control feasible
				result << "control feasible" << endl;
				break;
			} 
			else{	// control infeasible
				cut++;
				for(int i=0; i < timesteps; i++)
					v[i]=u[i].get(GRB_DoubleAttr_X);
				oldobj=obj;
				model.optimize();

				obj=model.get(GRB_DoubleAttr_ObjVal);;
				if((obj-oldobj)/oldobj<0.001)
					count++; 
				else
					count =0;
			}
		}

		result << "cuts:" << cut << endl;
		if(model.get(GRB_IntAttr_Status)==GRB_OPTIMAL){
			result << "[End_control]" << endl; 
			for(int i=0; i < timesteps; i++){
				v[i]=u[i].get(GRB_DoubleAttr_X);
				result << u[i].get(GRB_StringAttr_VarName) << ":" << v[i] << endl;
			}
			obj=objective(v);
			result << "[End_objective]" << endl;
				result << obj << endl;
			result << "[Runtime]" << endl;
			result << time+model.get(GRB_DoubleAttr_Runtime) << endl;
			if(count==3)
				result << "Ended due to low relative change" << endl;
		}
		else{
			cout << "Optimization was stopped with status " << status << endl;
			result <<"Optimization was stopped with status " << status << endl;
			result << "[End_control]" << endl; 
			for(int i=0; i < timesteps; i++)
				result << u[i].get(GRB_StringAttr_VarName) << ":" << v[i]<< endl;
			obj=objective(v);
			result << "[End_objective]" << endl;
				result << obj << endl;
			result << "[Runtime]" << endl;
			result << time << endl;
		}
		result.close();	
		delete[] v;
	}		
	else if(status==GRB_OPTIMAL && mod!=2){
		fstream result;
		result.precision(10);
		result.open(resfile, ios::out | ios::app); 
		double* v=new double[timesteps]; 
		result << "\n[End_control]" << endl; 
		for(int i=0; i < timesteps; i++){
				v[i]=u[i].get(GRB_DoubleAttr_X);
				result << u[i].get(GRB_StringAttr_VarName) << ":" << v[i] << endl;
		}
		result << "[End_objective]" << endl; 
		double obj=objective(v);
		result << obj << endl;
		
		result << "[Runtime]" << endl;
		result << model.get(GRB_DoubleAttr_Runtime) << endl;
		result.close();	
		delete[] v;
	}
	else
		cout << "Optimization was stopped with status " << status << endl;

	
	delete[] y; 
	delete[] u;
	delete[] z; 
	delete env;
}

#endif // GUROBI_HH
