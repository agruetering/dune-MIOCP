#ifndef ADJOINTHEATFEM_HH
#define ADJOINTHEATFEM_HH

#include "linearheatfem.hh"

using namespace Dune; 

/* a local operator for solving the adjoint equation
 *
 *  -\partial_t p - \Delta p = q   in Q=\Omega x (0,T)
 *                         p = g   on \Sigma_D=\Gamma_D x (0,T)
 *          \nabla y \cdot n = j   on \Sigma_N=\Gamma_N x (0,T)
 *                      p(T) = pT
 *
 * with conforming finite elements on all types of grids in any dimension
 *
 *
 * residual formulation 
 *
 * r_h(y_h,\Phi_i) = -d/dt \int_\Omega p_h \Phi_i dx 
 * 			+ \int_\Omega \nabla p_h \nabla \Phi_i dx 
 * 			+ \int_{\Sigma_N} j \Phi_i dx
 *			- \int_\Omega q \Phi_i dx 
 *			= 0					i=1,...,s
 */

 /* spatial local operator (identical with the one for the heat equation)
 */ 


/* temporal local operator (L_2 integral)
 *		-d/dt \int_\Omega p_h \Phi_i dx
 * \tparam FEM      Type of a finite element map
 */
template<class FEM>
class NL2 : public L2<FEM>
{
  	typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  	PDELab::LocalBasisCache<LocalBasis> cache;

public:
  	// pattern assembly flags
  	enum { doPatternVolume = true };

  	// residual assembly flags
  	enum { doAlphaVolume = true };

  	// volume integral (-\int_\Omega p_h \Phi_i dx) 
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

        		// evaluate p_h
        		RF p=0.0;
        		for (size_t i=0; i<lfsu.size(); i++) 
				p += x(lfsu,i)*phihat[i];

        		// integrate -p_h \Phi_i
        		RF factor = ip.weight() * geo.integrationElement(ip.position());
        		for (size_t i=0; i<lfsv.size(); i++)
          			r.accumulate(lfsv,i,-p*phihat[i]*factor);
      		}
  	}

	// jacobian contribution of volume term
	// depending on test functions (LSFV finite element space)
	// depending on ansatz functions (LSFU finite element space)
  	template<class EG, class LFSU, class X, class LFSV,class M>
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

        		// integrate -\Phi_j*\Phi_i
        		RF factor = ip.weight() * geo.integrationElement(ip.position());
        		for (size_t j=0; j<lfsu.size(); j++)
			{
          			for (size_t i=0; i<lfsv.size(); i++)
				{
            				mat.accumulate(lfsv,i,lfsu,j,-phihat[j]*phihat[i]*factor);
      				}
  			}
		}
	}
	
};


#endif // ifndef ADJOINTHEATFEM_HH
