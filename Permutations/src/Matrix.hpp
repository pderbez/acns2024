#ifndef DEF_MATRIX
#define DEF_MATRIX

#include <vector>
#include <map>

#include "GFElement.hpp"


class Matrix
{
public:
  Matrix() = default;
  Matrix(std::vector<std::vector<std::pair<GFElement, int>>> const &);
  Matrix(Matrix const & m) : nbcols (m.nbcols), nblines (m.nblines), front(m.front), columns(m.columns), space(m.space) {
    lines.reserve(nblines);
    for (unsigned l = 0; l < nblines; ++l) lines.emplace_back(&space[0] + (m.lines[l]- &m.space[0]));
  };
  Matrix(Matrix &&) = default;

  Matrix & operator=(Matrix const & m) {
      nbcols = m.nbcols;
      nblines = m.nblines;
      front = m.front;
      columns = m.columns;
      space = m.space;
      lines = std::vector<GFElement*> ();
      lines.reserve(nblines);
      for (unsigned l = 0; l < nblines; ++l) lines.emplace_back(&space[0] + (m.lines[l]- &m.space[0]));
      return *this;
  };

  Matrix & operator=(Matrix &&) = default;

  unsigned checkZ(double * X);


  bool dim2(int x1, int x2) const;

  void echelonizeOn(std::vector<int> const &);

  Matrix extract(std::vector<int> const &);

  GFElement const & operator()(unsigned i, unsigned j) const { return lines[i][j];};
  GFElement & operator()(unsigned i, unsigned j) {return lines[i][j];};

  int getFront(unsigned i) const {return front[i];};
  int getColumns(unsigned i) const {return columns[i];};

  void printLine(unsigned i);

  friend std::ostream& operator<<( std::ostream &flux, Matrix const& mat);

  bool setAsPivot(int x, unsigned start);
  bool isLinear(int x, unsigned start);

  unsigned nbcols;
  unsigned nblines;

private:
  std::vector<int> front;
  std::vector<int> columns;

  std::vector<GFElement> space;
  std::vector<GFElement*> lines;

  void swapLineColumn(unsigned l, unsigned c);






};

#endif
