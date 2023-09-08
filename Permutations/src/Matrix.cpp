#include <set>
#include <map>

#include "Matrix.hpp"

using namespace std;

void printVar(ostream &flux, int x) {
  if (x >= 0) {
    int i = x%16;
    int v = (x/16)%3;
    int r = (x/16)/3;
    if (v == 0) flux << "X";
    else if (v == 1) flux << "Y";
    else flux << "K";
    flux << r << "[" << i << "]";
  }
  else {
    flux << "S(";
    printVar(flux, -x);
    flux << ")";
  }
}

Matrix::Matrix(vector<vector<pair<GFElement, int>>> const & sys) {
  set<int> variables;
  for (auto const & v : sys) {
    for (auto const & p : v) variables.emplace(p.second);
  }
  unsigned n = sys.size();
  unsigned m = 0;
  map<int, unsigned> mapvars;
  for (auto x : variables) {
    mapvars[x] = m++;
  }
  vector<int> v (variables.begin(), variables.end());
  vector<vector<GFElement>> mat (n, vector<GFElement> (m, 0));
  for (unsigned i = 0; i < n; ++i) {
    for (auto const & p : sys[i]) mat[i][mapvars[p.second]] = p.first;
  }

  // for (unsigned i = 0; i < n; ++i) {
  //   for (unsigned j = 0; j < m; ++j) {
  //     if (mat[i][j] == 0) continue;
  //     cout << " + " << mat[i][j] << ".";
  //     printVar(cout, v[j]);
  //   }
  //   cout << endl;
  // }
  // cout << " -------------- " << endl;

  for (unsigned p = 0; p < n; ++p) {
    unsigned c = p;
    while (c < m && mat[p][c] == 0) ++c;
    if (c == m) swap(mat[p--], mat[--n]);
    else {
      if (c != p) {
        for (unsigned i = 0; i < n; ++i) swap(mat[i][p], mat[i][c]);
        swap(v[p], v[c]);
      }
      auto const coefInv = mat[p][p].getInverse();
      for (c = p; c < m; ++c) mat[p][c] *= coefInv;
      for (unsigned l = 0; l < n; ++l) {
        if (mat[l][p] == 0 || l == p) continue;
        auto const coef = mat[l][p];
        for (c = p; c < m; ++c) mat[l][c] -= coef*mat[p][c];
      }
    }
  }

  // for (unsigned i = 0; i < n; ++i) {
  //   for (unsigned j = 0; j < m; ++j) {
  //     if (mat[i][j] == 0) continue;
  //     cout << " + " << mat[i][j] << ".";
  //     printVar(cout, v[j]);
  //   }
  //   cout << endl;
  // }
  // cout << " -------------- " << endl;

  front = vector<int> (n);
  for (unsigned i = 0; i < n; ++i) front[i] = v[i];
  columns = vector<int> (m-n);
  for (unsigned i = n; i < m; ++i) columns[i-n] = v[i];
  space = vector<GFElement> (n*(m-n));
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = n; j < m; ++j) space[i*(m-n) + j-n] = mat[i][j];
  }
  nbcols = m-n;
  nblines = n;
  lines = vector<GFElement *> (n);
  for (unsigned i = 0; i < n; ++i) lines[i] = space.data() + i*nbcols;
}

void Matrix::swapLineColumn(unsigned l, unsigned c) {
  auto coef = (*this)(l,c).getInverse();
  for (unsigned i = 0; i < nbcols; ++i) (*this)(l,i) *= coef;
  (*this)(l,c) = 0;
  for (unsigned i = 0; i < nblines; ++i) {
    if ((*this)(i,c) == 0) continue;
    auto coef2 = (*this)(i,c);
    for (unsigned j = 0; j < nbcols; ++j) (*this)(i,j) += coef2*(*this)(l,j);
    (*this)(i,c) = coef*coef2;
  }
  (*this)(l,c) = coef;
  swap(front[l], columns[c]);
}

unsigned Matrix::checkZ(double * X) {
  for (unsigned l = 0; l < nblines; ++l) {
    auto const u = abs(front[l]);
    if (X[u] < 0.5) {
      unsigned c = 0;
      while (c < nbcols && ((*this)(l,c) == 0 || X[abs(columns[c])] < 0.5)) ++c;
      if (c < nbcols) swapLineColumn(l,c);
    }
  }

  for (unsigned l = 0; l < nblines; ++l) {
    if (X[abs(front[l])] < 0.5) continue;
    unsigned c = 0;
    while (c < nbcols && ((*this)(l,c) == 0 || X[abs(columns[c])] < 0.5)) ++c;
    if (c == nbcols) return l;
  }

  return nblines;
}

bool Matrix::setAsPivot(int x, unsigned start) {
  for (unsigned l = start; l < nblines; ++l) {
    if (front[l] == x) {
      swap(front[l], front[start]);
      swap(lines[l], lines[start]);
      return true;
    }
  }
  for (unsigned c = 0; c < nbcols; ++c) {
    if (columns[c] == x) {
      for (unsigned l = start; l < nblines; ++l) {
        if ((*this)(l,c) != 0) {
          swapLineColumn(l,c);
          swap(front[l], front[start]);
          swap(lines[l], lines[start]);
          return true;
        }
      }
      return false;
    }
  }
  return false;
}

bool Matrix::isLinear(int x, unsigned start) {
  if (setAsPivot(x, start)) {
    return !setAsPivot(-x, start+1);
  }
  else return true;
}

void Matrix::printLine(unsigned i) {
  printVar(cout, front[i]);
  for (unsigned j = 0; j < nbcols; ++j) {
    if ((*this)(i,j) != 0) {
      cout << " + ";
      if ((*this)(i,j) != 1) cout << (*this)(i,j) << ".";
      //flux << mat.columns[j];
      printVar(cout, columns[j]);
    }
  }
  cout << endl;
}

ostream& operator<<( ostream &flux, Matrix const& mat)
{
	for (unsigned i = 0; i < mat.nblines; ++i) {
    //flux << mat.front[i];
    printVar(flux, mat.front[i]);
    for (unsigned j = 0; j < mat.nbcols; ++j) {
      if (mat(i,j) != 0) {
        flux << " + ";
        if (mat(i,j) != 1) flux << mat(i,j) << ".";
        //flux << mat.columns[j];
        printVar(flux, mat.columns[j]);
      }
    }
    flux << endl;
  }
	return flux;
}
