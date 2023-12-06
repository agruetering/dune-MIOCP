// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later

#ifndef CALLBACK_HH
#define CALLBACK_HH
#include <gurobi_c++.h>

using namespace std; 

/* helper functions for the gurobi optimizer 
*/

/* function to prove feasibility of control for the tailored convexification */
int callback (GRBModel* model, GRBVar* vars, int timesteps, int smax); 
/* function to calculate heuristic solutions */
void callheuristic (int timesteps, double* w, int smax); 

#endif
