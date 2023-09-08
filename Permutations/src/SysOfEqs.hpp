#ifndef DEF_SYSOFEQS
#define DEF_SYSOFEQS

#include "Matrix.hpp"
#include <vector>


Matrix AES128eqs(unsigned R, int *KPerm);

Matrix AES128eqs(unsigned R, int *KPerm, std::vector<std::vector<unsigned>> const & subkeys);

#endif
