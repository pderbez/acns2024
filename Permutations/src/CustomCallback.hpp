#ifndef DEF_BOOMAES
#define DEF_BOOMAES

#include <vector>
#include <map>
#include "/opt/gurobi/linux64/include/gurobi_c++.h"
//#include "/Library/gurobi952/macos_universal2/include/gurobi_c++.h"
#include "Matrix.hpp"


class mycallback: public GRBCallback
{
public:
	Matrix mat;
	unsigned numvars;
	GRBVar * varsD;

	mycallback(int xnumvars, GRBVar* xvars, Matrix const & eqs) {
		numvars = xnumvars;
		varsD = xvars;
		mat = eqs;
	}
protected:
	void callback();
};

#endif
