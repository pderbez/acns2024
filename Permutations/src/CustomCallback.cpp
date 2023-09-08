#include "CustomCallback.hpp"

using namespace std;

void mycallback::callback()
{
	static unsigned mylazies = 0;
	try {
		if (where == GRB_CB_MIPSOL) {

			bool flagAddLazy = false;
			double * X = getSolution(varsD, numvars);
			unsigned l = mat.checkZ(X);
			if (l < mat.nblines) {
				GRBLinExpr e = 1-varsD[abs(mat.getFront(l))];
				for (unsigned c = 0; c < mat.nbcols; ++c) {
					if (mat(l, c) != 0) e += varsD[abs(mat.getColumns(c))];
				}
				addLazy(e >= 1);
				flagAddLazy = true;
			}
			delete[] X;
			if (flagAddLazy) cout << "\r" << ++mylazies;
		}

	} catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "Error during callback" << endl;
	}

}
