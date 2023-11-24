#include <type_traits>

using namespace Dune; 
using namespace std; 

/* parameter structure for outer approximation algorithm 
*/
template<class GridView, class Constraints,class Problem, int k, int number>
struct Param
{
	enum {dim=GridView::dimension}; // dimension of spatial grid
	enum {degree = k}; // degree of spatial ansatz and test functions
	enum {n=number}; // number of switches
	typedef GridView GV; // spatial grid view 
	typedef Constraints CON; // constraints (e.g. conforming dirichlet constraints)

	typedef Problem CGProblem; 

	typedef typename GV::ctype DF; // coordinate type of the spatial grid
	
	typedef typename PDELab::PkLocalFiniteElementMap<GV,DF,double,degree> FEM; // finite element map (for simplices)
  	typedef PDELab::ISTL::VectorBackend<> VBE; // ISTL backend to storage Dofs
  	typedef typename PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS; // grid function space for spatial ansatz and test functions (Galerkin method assumed)

	typedef typename PDELab::Backend::Vector<GFS,double> vtype; // type of Dof vectors
	typedef typename PDELab::DiscreteGridFunction<GFS,vtype> DGF; // discrete grid function
	typedef typename DGF::Traits::RangeType RangeType; // range type discrete grid function

     	typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
	typedef PDELab::LFSIndexCache<LFS> LFSCache; // local function space cache
	typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
	typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache

	typedef BlockVector<FieldVector<double,n>> ControlVector; // type of control vector
	/* values of all n switches j at all time points t_i saved as
	 * ( (u1(t_0),...,un(t_0)), (u1(t_1),...,un(t_1)),..., (u1(t_end),...,un(t_end)) )
	*/

	typedef BCRSMatrix<FieldMatrix<double,1,1>>Matrix; // matrix type for switching constraints

	typedef AdjointProblem<GFS> AdjProblem; // problem interface for S*SGv

	GV& gv;  // spatial grid
	FEM fem; // finite element map
	GFS gfs; // grid function space for spatial test and ansatz function (Galerkin method assumed)

	HeatProblemInterface<GridView,1>& heat; //  problem class for Sf (with boundary data)
	AdjointProblemInterface<GridView>& YdProblem; // problem class for yd (with boundary data)
	vector<double> dt;
	double alpha; // Thikonov term
	double rho; // parameter to update active switching constraints 
	int iter_outer; // maximal number of outer approximation iterations for each subproblem
	int iter_activeset; // maximal number of active set iterations
	int iter_lin; // maximal number of iterations for linear system 
	double tol; // absolute tolerance
	double reduction; // error reduction factor
	double timeouter; // time limit for outer approximation of each subproblem

	bool reopt; // reoptimization flag 

	/* constructor equidistant time grid */
	Param(GV gv_, 
	HeatProblemInterface<GridView,1>& heat_,
 	AdjointProblemInterface<GridView>& desired, 
	double alpha_, 
	double rho_,
	double T,
	int timesteps, 
	int iter_active, int iter, double tol_, double factor,
	int iter_o=500, 
	double limit=84600, 
	bool ropt=1) : 	
		gv(gv_),  
		fem(FEM(gv)), 
		gfs(GFS(gv,fem)), 
		heat(heat_),	
		YdProblem(desired),  
		dt(vector<double>(timesteps,T/timesteps)), 
		alpha(alpha_),
		rho(rho_), 
		iter_activeset(iter_active), 
		iter_lin(iter),
		tol(tol_),
		reduction(factor),
		iter_outer(iter_o), 
		timeouter(limit), 
		reopt(ropt)
	{}

	/* constructor non equidistant time grid */
	Param(GV gv_, 
	HeatProblemInterface<GridView,1>& heat_,
 	AdjointProblemInterface<GridView>& desired, 
	double alpha_, 
	double rho_,
	vector<double> dt_, 
	int iter_active, int iter, double tol_, double factor,
	int iter_o=500, 
	double limit=84600, 
	bool ropt=1) : 	
		gv(gv_),  
		fem(FEM(gv)), 
		gfs(GFS(gv,fem)), 
		heat(heat_),	
		YdProblem(desired),  
		dt(dt_), 
		alpha(alpha_), 
		rho(rho_),
		iter_activeset(iter_active), 
		iter_lin(iter),
		tol(tol_),
		reduction(factor),
		iter_outer(iter_o), 
		timeouter(limit), 
		reopt(ropt)
	{}

	
	
};

/* structure for output data of outer approximation algorithm 
*/
struct OutputData{
	int status; 
	/* status overview: 
	 * -3: maximum iterations of active set method reached
	 * -2: time limit reached
	 * -1: maximum iterations of outer approximation reached
	 * 0: optimum found
	 */
	int iter_lin; 
	int iter_activeset; 
	int iter_cut;
	double time_preprocess;
	vector<double> lower_bound; 
	vector<double> time_outer_solve; 
	vector<double> time_lin_solve; 
	vector<double> time_lin_assemble; 
	vector<double> time_outer_separation; 
	OutputData() {
		iter_lin = 0; 
		iter_activeset = 0; 
		iter_cut = 0;
		time_preprocess=0;
		lower_bound=vector<double>(); 
		time_outer_solve=vector<double>();
		time_outer_separation=vector<double>();
		time_lin_solve=vector<double>();
		time_lin_assemble=vector<double>();
	}
};

