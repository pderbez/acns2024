#include "SysOfEqs.hpp"

using namespace std;


// system of equations of AES-128
Matrix AES128eqs(unsigned R, int *KPerm) {
  GFElement mc [4] = {2, 3, 1, 1};
  vector<unsigned> shiftRows ({0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12});
  vector<unsigned> shiftRowsInv (16);
  for (unsigned i = 0; i < 16; ++i) shiftRowsInv[shiftRows[i]] = i;

  vector<vector<pair<GFElement, int>>> sys;

  // X = 3r + 0
  // S(X) = -(3*r + 0)
  // Y (after MC) = 3*r + 1
  // K = 3*r + 2


  for (unsigned r = 1; r < R; ++r) {
    for (unsigned c = 0; c < 4; ++c) {
      for (unsigned i = 0; i < 4; ++i) {
        vector<pair<GFElement, int>> v;
        for (unsigned j = 0; j < 4; ++j) v.emplace_back(mc[(j + 4 - i)%4], -((3*r)*16 + shiftRowsInv[4*j + c]));
        v.emplace_back(1, (3*r+1)*16 + 4*i + c);
        sys.emplace_back(move(v));
      }
    }
  }

  for (unsigned r = 1; r < R; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      vector<pair<GFElement, int>> v;
      v.emplace_back(1, (3*r+1)*16 + i);
      v.emplace_back(1, (3*r+2)*16 + i);
      v.emplace_back(1, (3*r+3)*16 + i);
      sys.emplace_back(move(v));
    }
  }
    // Key schedule relations
  for (unsigned r = 2; r < R; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      if (KPerm[i] != 16) {
        vector<pair<GFElement, int>> v;
        v.emplace_back(1, 16*(3*r + 2) + KPerm[i]);
        v.emplace_back(1, 16*(3*(r-1) + 2) + i);
        sys.emplace_back(move(v));
      }
    }
  }
  return Matrix(sys);
}

// system of equations of AES-128
Matrix AES128eqs(unsigned R, int *KPerm, vector<vector<unsigned>> const & subkeys) {
  GFElement mc [4] = {2, 3, 1, 1};
  vector<unsigned> shiftRows ({0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12});
  vector<unsigned> shiftRowsInv (16);
  for (unsigned i = 0; i < 16; ++i) shiftRowsInv[shiftRows[i]] = i;

  vector<vector<pair<GFElement, int>>> sys;

  // X = 3r + 0
  // S(X) = -(3*r + 0)
  // Y (after MC) = 3*r + 1
  // K = 3*r + 2


  for (unsigned r = 1; r < R; ++r) {
    for (unsigned c = 0; c < 4; ++c) {
      for (unsigned i = 0; i < 4; ++i) {
        vector<pair<GFElement, int>> v;
        for (unsigned j = 0; j < 4; ++j) v.emplace_back(mc[(j + 4 - i)%4], -((3*r)*16 + shiftRowsInv[4*j + c]));
        v.emplace_back(1, (3*r+1)*16 + 4*i + c);
        sys.emplace_back(move(v));
      }
    }
  }

  for (unsigned r = 1; r < R; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      vector<pair<GFElement, int>> v;
      v.emplace_back(1, (3*r+1)*16 + i);
      v.emplace_back(1, (3*r+2)*16 + i);
      v.emplace_back(1, (3*r+3)*16 + i);
      sys.emplace_back(move(v));
    }
  }

  unsigned nb_cols_perm = subkeys[0].size()/4;

  for (unsigned c = 4; c < 4*R; ++c) {
    unsigned r_roundk = c/4;
    unsigned c_roundk = c%4;
    unsigned r_subk = (c-4)/nb_cols_perm;
    unsigned c_subk = (c-4)%nb_cols_perm;
    if (r_subk == 0) continue;
    for (unsigned l = 0; l < 4; ++l) {
      unsigned x = subkeys[r_subk][c_subk + nb_cols_perm*l];
      unsigned rx = 1;
      unsigned cx = x%nb_cols_perm;
      unsigned lx = x/nb_cols_perm;
      while (cx >= 4) {cx -= 4; rx += 1;}
      if (rx < R) {
        vector<pair<GFElement, int>> v;
        v.emplace_back(1, 16*(3*rx + 2) + 4*lx + cx);
        v.emplace_back(1, 16*(3*r_roundk + 2) + 4*l + c_roundk);
      }
    }
  }

  return Matrix(sys);
}
