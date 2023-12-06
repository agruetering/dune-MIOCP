// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#include <dune/MIOCP/callback.hh>
#include <dune/MIOCP/Dmax.hh>

#include<math.h>

int callback (GRBModel* model, GRBVar* vars, int timesteps, int smax) {
	/* return value:
	 * 0: control feasible
	 * 1: control infeasible; cut is added
	 */

	// call separation algorithm of switching polytope
	Dmax poly(smax);
	double u[timesteps];
	for(int i=0; i < timesteps; i++)
		u[i]=min(1.0,max(vars[i].get(GRB_DoubleAttr_X),0.0));

	double cut[timesteps];
	double b=floor(smax/2.0); 

	int status=poly.separate(timesteps,u,cut,b);
	
	if(status==1){
		GRBLinExpr lhs = 0;
		for(int i=0; i < timesteps; i++) {
			lhs+=cut[i]*vars[i];
		}
		model->addConstr(lhs,GRB_LESS_EQUAL,b);
		return 1;
	}
	else
		return 0;
}

void callheuristic (int timesteps,double* w, int smax)
{
	// call optimization algorithm of switching polytope
	Dmax poly(smax);
	double* c=new double[timesteps];
	for(int i=0; i < timesteps; i++)
		c[i]=0.5-w[i];
	poly.optimize(timesteps,c,w); 
	delete[] c;
}
	

