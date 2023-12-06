// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef HELPFUNCTIONS_HH
#define HELPFUNCTIONS_HH

using namespace Dune; 
using namespace std;


/* contribution G*v =(\int Psi_j(x)v(t,x) dx )_{1<=j<=n} 
*/
template<class GFS, class Problem, class Vector>
void Gstar(const GFS& gfs, // grid functions space 
const Problem& problem, // problem class containing form functions
const vector<typename PDELab::Backend::Vector<GFS,double>>& v, // Dof vector v(t) at every time point
Vector& y) // vector to store result
{
    typedef typename GFS::Traits::GridView GV; // spatial grid view
    typedef typename GFS::Traits::FiniteElementMap FEM; // finite element map
    typedef typename PDELab::Backend::Vector<GFS,double> Y; // type of Dof vector
    typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF; // discrete grid function

  
    typedef typename DGF::Traits::RangeType RangeType; // range type discrete grid function
    typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
    typedef PDELab::LFSIndexCache<LFS> LFSCache; // local function space  cache
    typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
    typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache
	
    // local function space and basis
    LFS lfs(gfs);
    LFSCache lfs_cache(lfs);
    BasisCache cache;

    // clear y
    y.resize(v.size()-1);
    y=0;

    // loop over all elements of leaf grid
    for (const auto& eg : elements(gfs.gridView())){

	// local function space on current element
	lfs.bind(eg); 
       	lfs_cache.update();
        auto geo = eg.geometry();
	auto ref = referenceElement(geo);

        // choose appropiate quadrature rule
        const int order = -1+2*gfs.finiteElementMap().maxLocalSize();
        auto rule = PDELab::quadratureRule(geo,order);

	// loop over quadrature points
	for (const auto& ip : rule)
	{
		decltype(ip.weight()) factor = ip.weight()* geo.integrationElement(ip.position());

		// evaluate basis functions
		auto phihat = cache.evaluateFunction(ip.position(),
               lfs.finiteElement().localBasis());

            	// evaluate Psi
        	auto psi=problem.Psi(eg,ip.position());

            	for (int i=0; i < v.size()-1; i++){
                	RangeType val(0.0); 
                	// evaluate v(t_i,x)
                	for (size_t j=0; j<phihat.size(); j++) {
				val+=v[i][lfs_cache.containerIndex(j)]*phihat[j]; 
			} 
			val*=factor;
                	y[i].axpy(val,psi);
            	}
        }
   }
}


/* contribution \Chi_I G*v =\Chi_I (\int Psi_j(x)v(t,x) dx ) 
*/
template<class GFS, class Problem, class Vector>
void Gstar(const GFS& gfs, // grid functions space 
const Problem& problem, // problem class containing form functions
const vector<typename PDELab::Backend::Vector<GFS,double>>& v, // Dof vector v(t) at every time point
const vector<pair<int,int>> I, // set of inactive box constraints
Vector& y ) // vector to store result
{
    typedef typename GFS::Traits::GridView GV; // spatial grid view
    typedef typename GFS::Traits::FiniteElementMap FEM; // finite element map
    typedef typename PDELab::Backend::Vector<GFS,double> Y; // type of Dof vector
    typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF;  // discrete grid function

    typedef typename DGF::Traits::RangeType RangeType; // range type discrete grid function
    typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
    typedef PDELab::LFSIndexCache<LFS> LFSCache; // local function space cache
    typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
    typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache
	
    // local function space and basis
    LFS lfs(gfs);
    LFSCache lfs_cache(lfs);
    BasisCache cache;

    // clear y
    y.resize(v.size()-1);
    y=0;

    // loop over all elements of leaf grid
    for (const auto& eg : elements(gfs.gridView())){

	// local function space on current element
	lfs.bind(eg); 
       	lfs_cache.update();
        auto geo = eg.geometry();

        // choose appropiate quadrature rule
        const int order = -1+2*gfs.finiteElementMap().maxLocalSize();
        auto rule = PDELab::quadratureRule(geo,order);

	// loop over quadrature points
	for (const auto& ip : rule)
	{
		decltype(ip.weight()) factor = ip.weight()* geo.integrationElement(ip.position());

		// evaluate basis functions
		auto phihat = cache.evaluateFunction(ip.position(),
                lfs.finiteElement().localBasis());

            	// evaluate Psi
        	auto psi=problem.Psi(eg,ip.position());

            	for (auto const& p: I){
                	RangeType val(0.0); 
                	// evaluate v(t_i,x)
                	for (size_t j=0; j<phihat.size(); j++) 
				val+=v[p.first][lfs_cache.containerIndex(j)]*phihat[j];
                	val*=psi[p.second]; 
                	val*=factor; 
                	y[p.first][p.second]+=val;
            	}
        }
    }
}


/**************************************************************************
 *
 * calculation of 1/2 ||y-y_d||^2_{L_2(Q)} 
 *
***************************************************************************/
template<class GFS>
void L2Deviation(const GFS& gfs, // spatial grid function space
const vector<double>& dt, // temporal grid 
const vector<typename PDELab::Backend::Vector<GFS,double>>& y, // Dof vector 
AdjointProblemInterface<typename GFS::Traits::GridView>& YdProblem, // problem class for yd
double& bound)
{
	typedef typename GFS::Traits::GridView GV; // spatial grid view
    	typedef typename GFS::Traits::FiniteElementMap FEM; // finite element map
    	typedef typename PDELab::Backend::Vector<GFS,double> Y; // type of Dof vector
   	 typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF;  // discrete grid function

    	typedef typename DGF::Traits::RangeType RangeType; // range type discrete grid function
    	typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
    	typedef PDELab::LFSIndexCache<LFS> LFSCache; // local function space cache
    	typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
    	typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache

	YdProblem.setdeltaT(dt);

    	LFS lfs(gfs);
    	LFSCache lfs_cache(lfs);
    	BasisCache cache;
	
    	bound=0; 
    	// 1/2 int_\Omega \int_[0,T] (y-y_d)^2 dt dx
    	// loop over all elements of leaf grid
    	for (const auto& eg : elements(gfs.gridView())){
		
        	lfs.bind(eg); 
        	lfs_cache.update();
       	 	auto geo = eg.geometry();

        	// choose appropiate quadrature rule
       		const int order =-1+2*gfs.finiteElementMap().maxLocalSize();
        	auto rule = PDELab::quadratureRule(geo,order);

        	// loop over quadrature points
    		for (const auto& ip : rule)
     		{
            		decltype(ip.weight()) factor = ip.weight()* geo.integrationElement(ip.position());

            		// evaluate basis functions
            		auto phihat = cache.evaluateFunction(ip.position(),
                             lfs.finiteElement().localBasis());

                    	// loop over all time points
			// linear ansatz functions in time: y_h(t)=y(t_i,x)+(t-t_i)/dt*(y(t_{i+1},x)-y(t_i,x))
           		for (int i=0; i < dt.size(); i++){
                        	// evaluate desired temperature
				double des[2]; 
				des[0]=YdProblem.rhs(eg,ip.position(),i);
				des[1]=YdProblem.rhs(eg,ip.position(),i+1);
               	 		RangeType val0(0.0); 
				RangeType val1(0.0); 
                		// evaluate y(t_i,x), y(t_{i+1},x)
                		for (size_t j=0; j<lfs.size(); j++){
				        val0+=y[i][lfs_cache.containerIndex(j)]*phihat[j];
					val1+=y[i+1][lfs_cache.containerIndex(j)]*phihat[j];
                		}
				
                		double term; 
				term=pow(val0,2)+pow(val1,2)+val0*val1-2*val0*des[0]-2*val1*des[1]-val0*des[1]-val1*des[0]+pow(des[0],2)+pow(des[1],2)+des[0]*des[1];
                		term*=factor; 
				term*=dt[i];
                        	bound+=term; 
            		}
        	}
	}  
	bound*=1.0/6.0;

}
#endif //ifned HELPFUNCTIONS_HH
