// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
/* feasible switching pattern interface */ 
class SwitchPoly
{
public:
	virtual int separate(int,const double*, double*, double&) const=0; 
	/* compute a violated constraint, if any 
	 * parameters: 
	 * - current switching pattern
	 * - coefficients of violated constraint, if any
	 * - right-hand side of violated constraint, if any
	 * return value: 
	 * - 0: no violated constraint found 
	 * - 1: violated inequality found
	 */

	virtual int optimize(int,const double*, double*) const=0; 
	/* linear optimization over switching polytope
	 * parameters: 
	 * - objective coefficients 
	 * - computed solution
	 * return value:
 	 * - status of the optimization
  	 *   -1: error
  	 *    0: problem solved
	 */
};

