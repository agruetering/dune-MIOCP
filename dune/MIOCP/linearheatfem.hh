// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef LINEARHEATFEM_HH
#define LINEARHEATFEM_HH

// dune pdelab includes
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>

using namespace Dune; 

/* a local operator for solving the heat equation
 *
 *  \partial_t y - \Delta y = q   	in Q=\Omega x (0,T)
 *                        y = g   	on \Sigma_D=\Gamma_D x (0,T)
 *         \nabla y \cdot n = j   	on \Sigma_N=\Gamma_N x (0,T)
 *                     y(0) = y0
 *
 *
 * residual formulation 
 *
 * r_h(y_h,\Phi_i) = d/dt \int_\Omega y_h \Phi_i dx 
 * 			+ \int_\Omega \nabla y_h \nabla \Phi_i dx 
 * 			+ \int_{\Gamma_N} j \Phi_i dx
 *			- \int_\Omega q \Phi_i dx 
 *			= 0					
 */

 /* spatial local operator
 * 			\int_\Omega \nabla y_h \nabla \Phi_i dx 
 * 			+ \int_{\Gamma_N} j \Phi_i dx
 *			- \int_\Omega q \Phi_i dx 
 *
 * \tparam Problem   Problem class providing all coefficient functions in the PDE
 * \tparam FEM      Type of a finite element map
 */
template<class Problem, class FEM> 
class LinearHeatFEM : 
	public PDELab::FullVolumePattern,
    	public PDELab::LocalOperatorDefaultFlags,
    	public PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
	 typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  	 PDELab::LocalBasisCache<LocalBasis> cache;
protected: 
	Problem& param; // problem class for boundary, right hand side data
public: 
	public:
  	// pattern assembly flags
 	 enum { doPatternVolume = true };

  	// residual assembly flags
  	enum { doLambdaVolume = true }; // local contribution right hand side integral
  	enum { doLambdaBoundary = true }; // local contribution of neuman boundary conditions
  	enum { doAlphaVolume = true };  // local contribution volume term

	
	LinearHeatFEM(Problem & param_) : 
		param(param_)
	{}


	// local contribution right hand side integral (- \int_\Omega q \Phi_i dx)
	// depending on test functions (LSFV finite element space)
	 template<class EG, class LFSV, class R>
  	void lambda_volume (const EG& eg, const LFSV& lfsv,
                      R& r) const
  	{
   		 // select quadrature rule
    		auto geo = eg.geometry();
    		const int order = 1+
      				2*lfsv.finiteElement().localBasis().order();
    		auto rule = PDELab::quadratureRule(geo,order);

    		// loop over quadrature points
    		for (const auto& ip : rule)
     		{
       			 // evaluate basis functions
        		auto& phihat = cache.evaluateFunction(ip.position(),
                             lfsv.finiteElement().localBasis());

        		// integrate -q \Phi_i
        		decltype(ip.weight()) factor = ip.weight()* geo.integrationElement(ip.position());
        		auto value=param.q(eg.entity(),ip.position());
        		for (size_t i=0; i<lfsv.size(); i++)
         		 	r.accumulate(lfsv,i,-value*phihat[i]*factor);
      		}
  	}

	// local contribution of neuman boundary condition (+ \int_{\Gamma_N} j \Phi_i dx) 
	// depending on test functions (LSFV finite element space)
	template<class IG, class LFSV, class R>
  	void lambda_boundary (const IG& ig, const LFSV& lfsv,
                        	R& r) const
  	{
    		// evaluate boundary condition type
    		auto localgeo = ig.geometryInInside();
    		auto facecenterlocal =referenceElement(localgeo).position(0,0);
    		bool isdirichlet=param.b(ig.intersection(),facecenterlocal);

    		// skip rest if we are on Dirichlet boundary
    		if (isdirichlet) return;

    		// select quadrature rule
    		auto globalgeo = ig.geometry();
    		const int order = 1+2*lfsv.finiteElement().localBasis().order();
    		auto rule = PDELab::quadratureRule(globalgeo,order);

    		// loop over quadrature points 
    		for (const auto& ip : rule)
      		{
        		// quadrature point in local coordinates of element
        		auto local = localgeo.global(ip.position());

        		// evaluate basis functions (assume Galerkin method)
        		auto& phihat = cache.evaluateFunction(local,
                             lfsv.finiteElement().localBasis());

        		// integrate j \Phi_i
       			 decltype(ip.weight()) factor = ip.weight()*globalgeo.integrationElement(ip.position());
        		 auto j = param.j(ig.intersection(),ip.position());
        		 for (size_t i=0; i<lfsv.size(); i++)
          			r.accumulate(lfsv,i,j*phihat[i]*factor);
      		}
	}

	// local contribution volume term (+ \int_\Omega \nabla y_h \nabla \Phi_i dx) 
	// depending on test functions (LSFV finite element space)
	// depending on ansatz functions (LSFU finite element space)
	template<class EG, class LFSU, class X,class LFSV, class R>
  	void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  	{
    		// dimension of grid
    		const int dim = EG::Entity::dimension;
		// range type of ansatz functions
    		typedef decltype(PDELab::makeZeroBasisFieldValue(lfsu)) RF;

    		// select quadrature rule
   		 auto geo = eg.geometry();
    		const int order = 1+2*lfsu.finiteElement().localBasis().order();
    		auto rule = PDELab::quadratureRule(geo,order);

    		// loop over quadrature points
    		for (const auto& ip : rule)
      		{

        		// evaluate gradient of shape functions
        		auto& gradphihat = cache.evaluateJacobian(ip.position(),
                             	lfsu.finiteElement().localBasis());

        		// transform gradients from reference element to grid element
        		const auto S = geo.jacobianInverseTransposed(ip.position());
        		auto gradphi = makeJacobianContainer(lfsu);
        		for (size_t i=0; i<lfsu.size(); i++)
          			S.mv(gradphihat[i][0],gradphi[i][0]);

        		// compute gradient of y_h
       			 FieldVector<RF,dim> grady(0.0);
        		 for (size_t i=0; i<lfsu.size(); i++)
          			grady.axpy(x(lfsu,i),gradphi[i][0]);

        		// integrate ( \nabla y_h)* \nabla \Phi_i 
        		auto factor = ip.weight()*geo.integrationElement(ip.position());
        		for (size_t i=0; i<lfsu.size(); i++)
          			r.accumulate(lfsu,i,grady*gradphi[i][0]*factor);
      		}
  	}

	// jacobian contribution of volume term
	// depending on test functions (LSFV finite element space)
	// depending on ansatz functions (LSFU finite element space)
	template<class EG, class LFSU,class X, class LFSV, class M>
  	void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        	M& mat) const
  	{	
		// value type of ansatz functions
    		typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

   		 // select quadrature rule
    		auto geo = eg.geometry();
    		const int order = 1+2*lfsu.finiteElement().localBasis().order();
    		auto rule = PDELab::quadratureRule(geo,order);

    		// loop over quadrature points
    		for (const auto& ip : rule)
      		{
       			// evaluate gradient of shape functions
        		auto& gradphihat = cache.evaluateJacobian(ip.position(),
                             	lfsu.finiteElement().localBasis());

        		// transform gradients from reference element to grid element
        		const auto S = geo.jacobianInverseTransposed(ip.position());
        		auto gradphi = makeJacobianContainer(lfsu);
        		for (size_t i=0; i<lfsu.size(); i++)
          			S.mv(gradphihat[i][0],gradphi[i][0]);

        		// integrate \nabla \Phi_j* \nabla \Phi_i
        		RF factor = ip.weight() * geo.integrationElement(ip.position());
        		for (size_t j=0; j<lfsu.size(); j++)
          			for (size_t i=0; i<lfsv.size(); i++)	
            				mat.accumulate(lfsv,i,lfsu,j,gradphi[j][0]*gradphi[i][0]*factor);
		}
	}	
}; 



/* temporal local operator (L_2 integral)
 *		d/dt \int_\Omega y_h \Phi_i dx
 *
 * \tparam FEM      Type of a finite element map
 */
template<class FEM>
class L2 : 
	public PDELab::FullVolumePattern,
    	public PDELab::LocalOperatorDefaultFlags,
    	public PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  	typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  	PDELab::LocalBasisCache<LocalBasis> cache;

public:
  	// pattern assembly flags
  	enum { doPatternVolume = true };

  	// residual assembly flags
  	enum { doAlphaVolume = true };

  	// volume integral (\int_T y_h \Phi_i dx) 
	// depending on test functions (LSFV finite element space)
	// depending on ansatz functions (LSFU finite element space)
  	template<class EG, class LFSU, class X,class LFSV, class R>
  	void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     	const LFSV& lfsv, R& r) const
  	{
   	 	// value type of ansatz functions
    		typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

    		// select quadrature rule
    		auto geo = eg.geometry();
    		const int order = 1+2*lfsu.finiteElement().localBasis().order();
    		auto rule = PDELab::quadratureRule(geo,order);

    		// loop over quadrature points
    		for (const auto& ip : rule)
      		{
       			 // evaluate basis functions
        		auto& phihat = cache.evaluateFunction(ip.position(),
                                 lfsu.finiteElement().localBasis());

        		// evaluate y_h
        		RF y=0.0;
        		for (size_t i=0; i<lfsu.size(); i++) 
				y += x(lfsu,i)*phihat[i];
			

        		// integrate y_h \Phi_i
        		RF factor = ip.weight() * geo.integrationElement(ip.position());
        		for (size_t i=0; i<lfsv.size(); i++)
          			r.accumulate(lfsv,i,y*phihat[i]*factor);
      		}
  	}

	// jacobian contribution of volume term
	// depending on test functions (LSFV finite element space)
	// depending on ansatz functions (LSFU finite element space)
  	template<class EG, class LFSU, class X, class LFSV, class M>
  	void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        	M& mat) const
  	{
    		// value type of ansatz functions
    		typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

   		 // select quadrature rule
    		auto geo = eg.geometry();
    		const int order = 1+2*lfsu.finiteElement().localBasis().order();
    		auto rule = PDELab::quadratureRule(geo,order);

    		// loop over quadrature points
    		for (const auto& ip : rule)
      		{
       			 // evaluate basis functions
        		auto& phihat =	cache.evaluateFunction(ip.position(),
                                 lfsu.finiteElement().localBasis());

        		// integrate \Phi_j*\Phi_i
        		RF factor = ip.weight() * geo.integrationElement(ip.position());
        		for (size_t j=0; j<lfsu.size(); j++)
			{
          			for (size_t i=0; i<lfsv.size(); i++)
				{
            				mat.accumulate(lfsv,i,lfsu,j,phihat[j]*phihat[i]*factor);
      				}
  			}
		}
	}
	
};


#endif // ifndef LINEARHEATFEM_HH
