#include <iostream>
#include <vector>
#include <set>
#include "CustomCallback.hpp"
#include "SysOfEqs.hpp"

//#include "/opt/gurobi/linux64/include/gurobi_c++.h"
//#include "/Library/gurobi952/macos_universal2/include/gurobi_c++.h"

using namespace std;

void printPerm(vector<int> & P) {
    for (unsigned i = 0; i < 16; ++i){
        if(P[i] == 16) cout << "* ";
        else cout << P[i] << " ";
    }
    cout << " " << endl;
}

unsigned modelAESMinSboxes(int R, vector<int> & KPerm, GRBEnv & env) {

    GRBModel model = GRBModel(env);

    auto const & size_perm = KPerm.size();
    auto const & nb_cols_perm = size_perm/4;


    // X = 3r + 0
    // S(X) = -(3*r + 0)
    // Y (after MC) = 3*r + 1
    // K = 3*r + 2

    // definition of variables
    vector<GRBVar> dX ((3*R + 1)*16); // before ARK

    for (unsigned i = 0; i < (3*R+1)*16; ++i) {
        dX[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }

    for (unsigned i = 0; i < 3*16; ++i) {
      model.addConstr(dX[i] == 0);
    }

    vector<vector<unsigned>> subkeys;

    {
      vector<unsigned> v (size_perm);
      for (unsigned x = 0; x < size_perm; ++x) v[x] = x;
      subkeys.emplace_back(v);
      while (size_perm*subkeys.size() < 16*R) {
        vector<unsigned> vv (size_perm);
        for (unsigned i = 0; i < size_perm; ++i) vv[KPerm[i]] = v[i];
        subkeys.emplace_back(vv);
        v = move(vv);
      }
    }

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
          model.addConstr(dX[16*(3*rx + 2) + 4*lx + cx] == dX[16*(3*r_roundk + 2) + 4*l + c_roundk]);
        }
      }
    }


    for (unsigned r = 1; r < R; ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 1)  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 0) + ((c + i)%4) + 4*i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);
      }
    }

    for (unsigned r = 1; r < R ; ++r) {
      //ARK (sX[r+1], kX[r+1] and zX[r+1])
      for (unsigned i = 0; i < 16; ++i) {
        model.addConstr(1-dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + 1-dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + 1-dX[16*(3*(r+1))  + i] >= 1);
      }
    }

    GRBLinExpr obj = 0;
    for (unsigned r = 1; r <= R; ++r) {
      for (unsigned i = 0; i < 16; ++i) obj += dX[16*(3*r) + i];
    }

    model.setObjective(obj, GRB_MINIMIZE);



    model.addConstr(obj >= 1);

    auto mat = AES128eqs(R, KPerm.data(), subkeys);
    //auto mat = AES128eqs(R, KPerm.data());

    mycallback cb ((3*R+1)*16, dX.data(), mat);
    model.setCallback(&cb);
    model.set(GRB_IntParam_OutputFlag , 0);
    model.set(GRB_IntParam_LazyConstraints , 1);
    // model.set(GRB_IntParam_PoolSolutions, 2000000);
    // model.set(GRB_DoubleParam_PoolGap, 0.001);
    // model.set(GRB_IntParam_PoolSearchMode, 2);

    //printPerm(KPerm);

    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);

    if (nSolutions > 0) {

      return model.getObjective().getValue();

    }

    return 0;
}

vector<vector<unsigned>> modelAES128(int R, vector<int> & KPerm, int nrSboxesWanted, GRBEnv & env) {

    GRBModel model = GRBModel(env);

    auto const & size_perm = KPerm.size();
    auto const & nb_cols_perm = size_perm/4;


    // X = 3r + 0
    // S(X) = -(3*r + 0)
    // Y (after MC) = 3*r + 1
    // K = 3*r + 2

    // definition of variables
    vector<GRBVar> dX ((3*R + 1)*16); // before ARK

    for (unsigned i = 0; i < (3*R+1)*16; ++i) {
        dX[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }

    for (unsigned i = 0; i < 3*16; ++i) {
      model.addConstr(dX[i] == 0);
    }


    // Contraints for permutation-based key schedule
    // number of active bytes per key must be the same
    // {
    //   GRBLinExpr e = 0;
    //   for (unsigned i = 0; i < 16; ++i) e += dX[16*(3*1 + 2) + i];
    //   for (unsigned r = 2; r < R; ++r) {
    //     GRBLinExpr e1 = 0;
    //     for (unsigned i = 0; i < 16; ++i) e1 += dX[16*(3*r + 2) + i];
    //     model.addConstr(e1 == e);
    //   }
    // }

    vector<vector<unsigned>> subkeys;

    {
      vector<unsigned> v (size_perm);
      for (unsigned x = 0; x < size_perm; ++x) v[x] = x;
      subkeys.emplace_back(v);
      while (size_perm*subkeys.size() < 16*R) {
        vector<unsigned> vv (size_perm);
        for (unsigned i = 0; i < size_perm; ++i) vv[KPerm[i]] = v[i];
        subkeys.emplace_back(vv);
        v = move(vv);
      }
    }

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
          model.addConstr(dX[16*(3*rx + 2) + 4*lx + cx] == dX[16*(3*r_roundk + 2) + 4*l + c_roundk]);
          //cout << r_roundk << ", " << c_roundk << ", " << l << " --> " << rx << ", " << cx << ", " << lx << endl;
        }
      }
    }

/*     for (auto const & v : subkeys) {
      for (unsigned i = 0; i < size_perm; ++i) {
        if (i % (size_perm/4) == 0) cout << endl;
        cout << v[i] << " ";
      }
      cout << endl;
    }
    getchar(); */


    for (unsigned r = 1; r < R; ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 1)  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 0) + ((c + i)%4) + 4*i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);
      }
    }

    for (unsigned r = 1; r < R ; ++r) {
      //ARK (sX[r+1], kX[r+1] and zX[r+1])
      for (unsigned i = 0; i < 16; ++i) {
        model.addConstr(1-dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + 1-dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + 1-dX[16*(3*(r+1))  + i] >= 1);
      }
    }

    GRBLinExpr obj = 0;
    for (unsigned r = 1; r <= R; ++r) {
      for (unsigned i = 0; i < 16; ++i) obj += dX[16*(3*r) + i];
    }

    GRBLinExpr obj2 = 0;
    //for (unsigned i = 0; i < 16; ++i) obj2 += dX[16*(3*1 + 2) + i];
    for (unsigned i = 0; i < size_perm; ++i) {
      unsigned r = 1;
      unsigned c = i%nb_cols_perm;
      unsigned l = i/nb_cols_perm;
      while (c >= 4) {c -= 4; r += 1;}
      if (r < R) obj2 += dX[16*(3*r + 2) + 4*l + c];
    }

    model.setObjective(obj2, GRB_MINIMIZE);

    //model.addConstr(obj2 == 2);




    /*GRBLinExpr obj3 = 0;

    for (unsigned r = 1; r < R; ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        GRBLinExpr e1 = 0;
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        GRBVar g = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        for (unsigned i = 0; i < 4; ++i) e1 += dX[16*(3*r + 0) + ((c + i)%4) + 4*i];
        model.addConstr(e1 <= 4*f + 4*g);
        GRBLinExpr e2 = 0;
        for (unsigned i = 0; i < 4; ++i) e2 += dX[16*(3*(r+1) + 0)  + c  + 4*i];
        model.addConstr(e1 + e2 >= 4*(1-f) - 4*(1-g));
        obj3 += f;
      }
    }

    model.setObjective(obj3, GRB_MINIMIZE);*/

    model.addConstr(obj >= 1);
    model.addConstr(obj <= nrSboxesWanted-1);

    auto mat = AES128eqs(R, KPerm.data(), subkeys);
    //auto mat = AES128eqs(R, KPerm.data());

    mycallback cb ((3*R+1)*16, dX.data(), mat);
    model.setCallback(&cb);
    model.set(GRB_IntParam_OutputFlag , 0);
    model.set(GRB_IntParam_LazyConstraints , 1);
    // model.set(GRB_IntParam_PoolSolutions, 2000000);
    // model.set(GRB_DoubleParam_PoolGap, 0.001);
    // model.set(GRB_IntParam_PoolSearchMode, 2);

    //printPerm(KPerm);

    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);


    vector<vector<unsigned>> v_res;
    if (nSolutions > 0) {
      for (auto e = 0; e < 1; e++) {
          model.set(GRB_IntParam_SolutionNumber, e);

          vector<unsigned> res;
          unsigned c_sub = 0;
          for (unsigned r = 1; r < R; ++r) {
            for (unsigned c = 0; c < 4; ++c) {
              for (unsigned l = 0; l < 4; ++l) {
                if (dX[16*(3*r + 2) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) {
                  res.emplace_back(nb_cols_perm*l + c_sub);
                  //cout << r << "," << c << "," << l << " -> " << nb_cols_perm*l + c_sub << endl;
                }
              }
              c_sub += 1;
              if (c_sub == nb_cols_perm) {
                c_sub = 0;
                v_res.emplace_back(res);
                res.clear();
              }
            }
          }
          /*
          while (c_sub%nb_cols_perm != 0) {
            for (unsigned l = 0; l < 4; ++l) res.emplace_back(nb_cols_perm*l + c_sub);
            ++c_sub;
          }*/
          if (!res.empty()) v_res.emplace_back(res);

          // for (unsigned i = 0; i < size_perm; ++i) cout << i << ": " << KPerm[i] << " - ";
          // cout << endl;
          //
          // for (unsigned r = 0; r < subkeys.size(); ++r) {
          //   for (unsigned l = 0; l < 4; ++l) {
          //     for (unsigned c = 0; c < nb_cols_perm; ++c) cout << subkeys[r][nb_cols_perm*l + c] << " ";
          //     cout << endl;
          //   }
          //   cout << endl << endl;
          // }
          // getchar();


          // cout << endl;
          //
          // for (unsigned r = 1; r < R; ++r) {
          //   for (unsigned l = 0; l < 4; ++l) {
          //     for (unsigned c = 0; c < 4; ++c) {
          //       if (dX[16*(3*r + 0) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
          //       else cout << " 0 ";
          //     }
          //     cout << "   ";
          //     for (unsigned c = 0; c < 4; ++c) {
          //       if (dX[16*(3*r + 1) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
          //       else cout << " 0 ";
          //     }
          //     cout << "  || ";
          //     for (unsigned c = 0; c < 4; ++c) {
          //       if (dX[16*(3*r + 2) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
          //       else cout << " 0 ";
          //     }
          //     cout << endl;
          //   }
          //   cout << endl << endl;
          // }
          // getchar();




      }
    }


    return v_res;

}

vector<vector<unsigned>> modelAESBoomerang(int R, vector<int> & KPerm, GRBEnv & env) {

    GRBModel model = GRBModel(env);

    auto const & size_perm = KPerm.size();
    auto const & nb_cols_perm = size_perm/4;


    // X = 3r + 0
    // S(X) = -(3*r + 0)
    // Y (after MC) = 3*r + 1
    // K = 3*r + 2

    // definition of variables
    vector<GRBVar> dX ((3*R + 1)*16); // before ARK
    vector<GRBVar> dX2 ((3*R + 1)*16); // before ARK

    vector<GRBVar> dF ((3*R + 1)*16); // before ARK
    vector<GRBVar> dF2 ((3*R + 1)*16); // before ARK

    for (unsigned i = 0; i < (3*R+1)*16; ++i) {
        dX[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        dX2[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        dF[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        dF2[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

        model.addConstr(dX[i] >= dF[i]);
        model.addConstr(dX2[i] >= dF2[i]);

        model.addConstr(dF[i] + dF2[i] <= 1);
    }

    for (unsigned i = 0; i < 3*16; ++i) {
      model.addConstr(dX[i] == 0);
      model.addConstr(dX2[i] == 0);
    }

    for (unsigned r = 1; r < R; ++r) {
      for (unsigned i = 0; i < 16; ++i) {
        model.addConstr(dF[16*(3*r + 2) + i] == 0);
        model.addConstr(dF2[16*(3*r + 2) + i] == 0);
      }
    }

    {
      GRBLinExpr e = 0;
      for (unsigned r = 1; r <= R; ++r) {
        for (unsigned i = 0; i < 16; ++i) e += dX[16*(3*r) + i];
      }
      model.addConstr(e >= 1);
    }

    {
      GRBLinExpr e = 0;
      for (unsigned r = 1; r <= R; ++r) {
        for (unsigned i = 0; i < 16; ++i) e += dX2[16*(3*r) + i];
      }
      model.addConstr(e >= 1);
    }

    for (unsigned i = 0; i < 16; ++i) model.addConstr(dF[16*(3*1) + i] == 0);
    for (unsigned i = 0; i < 16; ++i) model.addConstr(dF[16*(3*R) + i] == dX[16*(3*R) + i]);

    for (unsigned i = 0; i < 16; ++i) model.addConstr(dF2[16*(3*R) + i] == 0);
    for (unsigned i = 0; i < 16; ++i) model.addConstr(dF2[16*(3*1) + i] == dX2[16*(3*1) + i]);

    vector<vector<unsigned>> subkeys;

    {
      vector<unsigned> v (size_perm);
      for (unsigned x = 0; x < size_perm; ++x) v[x] = x;
      subkeys.emplace_back(v);
      while (size_perm*subkeys.size() < 16*R) {
        vector<unsigned> vv (size_perm);
        for (unsigned i = 0; i < size_perm; ++i) vv[KPerm[i]] = v[i];
        subkeys.emplace_back(vv);
        v = move(vv);
      }
    }

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
          model.addConstr(dX[16*(3*rx + 2) + 4*lx + cx] == dX[16*(3*r_roundk + 2) + 4*l + c_roundk]);
          model.addConstr(dX2[16*(3*rx + 2) + 4*lx + cx] == dX2[16*(3*r_roundk + 2) + 4*l + c_roundk]);
        }
      }
    }


    for (unsigned r = 1; r < R; ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 1)  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 0) + ((c + i)%4) + 4*i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);

        GRBLinExpr e1 = 0;
        for (unsigned i = 0; i < 4; ++i) e1 += dF[16*(3*r + 1)  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) model.addConstr(e1 >= 4*dF[16*(3*r + 0) + ((c + i)%4) + 4*i]);
      }

      for (unsigned c = 0; c < 4; ++c) {
        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += dX2[16*(3*r + 1)  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) e += dX2[16*(3*r + 0) + ((c + i)%4) + 4*i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);

        GRBLinExpr e1 = 0;
        for (unsigned i = 0; i < 4; ++i) e1 += dF2[16*(3*r + 0) + ((c + i)%4) + 4*i];
        for (unsigned i = 0; i < 4; ++i) model.addConstr(e1 >= 4*dF2[16*(3*r + 1)  + c  + 4*i]);
      }
    }

    for (unsigned r = 1; r < R ; ++r) {
      //ARK (sX[r+1], kX[r+1] and zX[r+1])
      for (unsigned i = 0; i < 16; ++i) {
        model.addConstr(1-dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + 1-dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + 1-dX[16*(3*(r+1))  + i] >= 1);

        model.addConstr(dF[16*(3*r + 1)  + i] <= dF[16*(3*(r+1))  + i]);

        model.addConstr(1-dX2[16*(3*r + 1)  + i] + dX2[16*(3*r + 2)  + i] + dX2[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX2[16*(3*r + 1)  + i] + 1-dX2[16*(3*r + 2)  + i] + dX2[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX2[16*(3*r + 1)  + i] + dX2[16*(3*r + 2)  + i] + 1-dX2[16*(3*(r+1))  + i] >= 1);

        model.addConstr(dF2[16*(3*r + 1)  + i] >= dF2[16*(3*(r+1))  + i]);

        model.addConstr(dF[16*(3*(r+1))  + i] + dF2[16*(3*r + 1)  + i] >= 1);

      }
    }

    GRBLinExpr obj = 0;
    for (unsigned r = 1; r < R; ++r) {
      for (unsigned i = 0; i < 16; ++i) {
        GRBVar pBCT = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(pBCT + 1-dF[16*(3*(r+1))  + i] + 1 - dF2[16*(3*r + 1)  + i] >= 1);
        GRBVar pDDT2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(pDDT2 + dF[16*(3*(r+1))  + i] + 1 - dX[16*(3*(r+1))  + i] + 1 - dX2[16*(3*r + 1)  + i] >= 1);
        model.addConstr(pDDT2 + dF2[16*(3*r + 1)  + i] + 1 - dX2[16*(3*r + 1)  + i] + 1 - dX[16*(3*(r+1))  + i] >= 1);
        GRBVar pDDT = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(pDDT + 1 - dX[16*(3*(r+1))  + i] + dF[16*(3*(r+1))  + i] + dX2[16*(3*(r+1))  + i] >= 1);
        model.addConstr(pDDT + 1 - dX2[16*(3*(r+1))  + i] + dF2[16*(3*r+1)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        obj += 5.4*pBCT + 12*pDDT2 + 6*pDDT;
      }
    }

    model.addConstr(obj <= 127);


    model.setObjective(obj, GRB_MINIMIZE);




    auto mat = AES128eqs(R, KPerm.data(), subkeys);

    mycallback cb ((3*R+1)*16, dX.data(), mat);
    model.setCallback(&cb);
    //model.set(GRB_IntParam_OutputFlag , 0);
    model.set(GRB_IntParam_LazyConstraints , 1);
    // model.set(GRB_IntParam_PoolSolutions, 2000000);
    // model.set(GRB_DoubleParam_PoolGap, 0.001);
    // model.set(GRB_IntParam_PoolSearchMode, 2);

    //printPerm(KPerm);

    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);


    vector<vector<unsigned>> v_res;
    if (nSolutions > 0) {
      for (auto e = 0; e < 1; e++) {
          model.set(GRB_IntParam_SolutionNumber, e);

          vector<unsigned> res;
          unsigned c_sub = 0;
          for (unsigned r = 1; r < R; ++r) {
            for (unsigned c = 0; c < 4; ++c) {
              for (unsigned l = 0; l < 4; ++l) {
                if (dX[16*(3*r + 2) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) {
                  res.emplace_back(nb_cols_perm*l + c_sub);
                  //cout << r << "," << c << "," << l << " -> " << nb_cols_perm*l + c_sub << endl;
                }
              }
              c_sub += 1;
              if (c_sub == nb_cols_perm) {
                c_sub = 0;
                v_res.emplace_back(res);
                res.clear();
              }
            }
          }
          /*
          while (c_sub%nb_cols_perm != 0) {
            for (unsigned l = 0; l < 4; ++l) res.emplace_back(nb_cols_perm*l + c_sub);
            ++c_sub;
          }*/
          if (!res.empty()) v_res.emplace_back(res);

          // for (unsigned i = 0; i < size_perm; ++i) cout << i << ": " << KPerm[i] << " - ";
          // cout << endl;
          //
          // for (unsigned r = 0; r < subkeys.size(); ++r) {
          //   for (unsigned l = 0; l < 4; ++l) {
          //     for (unsigned c = 0; c < nb_cols_perm; ++c) cout << subkeys[r][nb_cols_perm*l + c] << " ";
          //     cout << endl;
          //   }
          //   cout << endl << endl;
          // }
          // getchar();


          // cout << endl;
          //
          for (unsigned r = 1; r < R; ++r) {
            for (unsigned l = 0; l < 4; ++l) {
              for (unsigned c = 0; c < 4; ++c) {
                if (dX[16*(3*r + 0) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
                else cout << " 0 ";
              }
              cout << "   ";
              for (unsigned c = 0; c < 4; ++c) {
                if (dX[16*(3*r + 1) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
                else cout << " 0 ";
              }
              cout << "  || ";
              for (unsigned c = 0; c < 4; ++c) {
                if (dX[16*(3*r + 2) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
                else cout << " 0 ";
              }
              cout << endl;
            }
            cout << endl << endl;
          }
          getchar();

          for (unsigned r = 1; r < R; ++r) {
            for (unsigned l = 0; l < 4; ++l) {
              for (unsigned c = 0; c < 4; ++c) {
                if (dX2[16*(3*r + 0) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
                else cout << " 0 ";
              }
              cout << "   ";
              for (unsigned c = 0; c < 4; ++c) {
                if (dX2[16*(3*r + 1) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
                else cout << " 0 ";
              }
              cout << "  || ";
              for (unsigned c = 0; c < 4; ++c) {
                if (dX2[16*(3*r + 2) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) cout << " 1 ";
                else cout << " 0 ";
              }
              cout << endl;
            }
            cout << endl << endl;
          }
          getchar();




      }
    }


    return v_res;

}

void testBoomerang256(unsigned R) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  vector<int> KPerm ({27,16,9,25,11,13,14,18,22,21,19,23,28,31,29,3,2,15,8,24,17,1,26,0,7,20,10,4,6,30,12,5});
  modelAESBoomerang(R+1, KPerm, env);
}

void testBoomerang192(unsigned R) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  vector<int> KPerm ({2,17,19,9,13,12,23,0,4,21,18,16,10,20,22,1,11,3,7,5,15,6,14,8});
  modelAESBoomerang(R+1, KPerm, env);
}

void testPerm128(unsigned R) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  
  //vector<int> KPerm ({5, 6, 7, 4, 8, 9, 10, 11, 12, 13, 14, 15, 2, 3, 0, 1});
  //vector<int> KPerm ({15, 1, 9, 6, 0, 14, 3, 4, 8, 5, 2, 7, 12, 13, 10, 11});
  //vector<int> KPerm ({2, 10, 9, 5, 4, 1, 6, 3, 13, 8, 14, 11, 15, 12, 0, 7});
  //vector<int> KPerm ({6, 0, 4, 9, 13, 10, 8, 3, 7, 12, 15, 14, 11, 5, 1, 2});
  vector<int> KPerm ({3, 15, 11, 8, 2, 1, 10, 5, 4, 0, 9, 7, 6, 12, 13, 14});
  
  unsigned nsbox = modelAESMinSboxes(R, KPerm, env);
  cout << "nsbox: " << nsbox << endl;
  auto vres = modelAES128(R, KPerm, nsbox, env);
  for (auto const & v : vres) {
   for (auto const x : v) cout << x << " ";
   cout << "| ";
  }
  cout << endl;
}


void searchMILP(vector<pair<unsigned, unsigned>> const & myconstraints, int size_perm) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();


  vector<int> KPerm(size_perm, size_perm);

  GRBModel model = GRBModel(env);
  model.set(GRB_IntParam_OutputFlag , 0);

  unsigned max_R = 0;
  for (auto const & p : myconstraints) {
    if (max_R < p.first) max_R = p.first;
  }

  if (max_R > 2) max_R -= 1;

  vector<vector<vector<GRBVar>>> P (max_R, vector<vector<GRBVar>> (size_perm, vector<GRBVar> (size_perm)));

  for (unsigned r = 1; r < max_R; ++r) {
    for (unsigned i = 0; i < size_perm; ++i) {
      for (unsigned j = 0; j < size_perm; ++j) P[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
    for (unsigned i = 0; i < size_perm; ++i) {
      GRBLinExpr e = 0;
      for (unsigned j = 0; j < size_perm; ++j) e += P[r][i][j];
      model.addConstr(e == 1);
    }

    for (unsigned i = 0; i < size_perm; ++i) {
      GRBLinExpr e = 0;
      for (unsigned j = 0; j < size_perm; ++j) e += P[r][j][i];
      model.addConstr(e == 1);
    }
  }


  auto const n_cols_perm = size_perm/4;

  // remove symetries
  if (size_perm == 16) {
  
	  for (unsigned j0 = 0; j0 < size_perm; ++j0) {
	    for (unsigned i = 1; i < 4; ++i) {
	      for (unsigned j = 0; j < size_perm; ++j) {
		unsigned jj = n_cols_perm*(j/n_cols_perm) + (j%n_cols_perm - i + n_cols_perm)%n_cols_perm;
		if (jj < j0) model.addConstr(P[1][i][j] <= 1 - P[1][0][j0]);
	      }
	    }
	  }
  }

  for (unsigned r = 2; r < (4*max_R+3)/n_cols_perm; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      for (unsigned j = 0; j < 16; ++j) {
        for (unsigned k = 0; k < 16; ++k) model.addConstr(P[r][i][k] + 1 >= P[r-1][i][j] + P[1][j][k]);
      }
    }
  }




  bool found = false;
  while (!found) {
    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);
    if (nSolutions == 0) {cout << "pas de solutions" << endl; getchar();}
    model.set(GRB_IntParam_SolutionNumber, 0);
    for (unsigned i = 0; i < size_perm; ++i) {
      unsigned j = 0;
      while (P[1][i][j].get(GRB_DoubleAttr_Xn) < 0.5) ++j;
      KPerm[i] = j;
    }
    for (unsigned i = 0; i < size_perm; ++i) cout << KPerm[i] << " ";
    cout << endl;

    unsigned j = 0;

    unsigned R = myconstraints[0].first;
    unsigned nbsboxes = myconstraints[0].second;

    auto v_res = modelAES128(R, KPerm, nbsboxes, env);
    while (v_res.empty() && ++j < myconstraints.size()) {
    	R = myconstraints[j].first;
    	nbsboxes = myconstraints[j].second;
    	v_res = modelAES128(R, KPerm, nbsboxes, env);
    }
    if (v_res.empty()) {
      cout << "yeah!!" << endl;
      for (unsigned i = 0; i < size_perm; ++i) cout << KPerm[i] << " ";
      getchar();
      GRBLinExpr e = 0;
      for (unsigned i = 0; i < size_perm; ++i) e += P[1][i][KPerm[i]];
      model.addConstr(e <= size_perm-1);
    }
    else {
      cout << "v_res: " << v_res.size() << endl;

      // for (unsigned i = 0; i < size_perm; ++i) cout << KPerm[i] << " ";
      // cout << endl;
      // for (unsigned r = 0; r < v_res.size(); ++r) {
      //   for (auto x : v_res[r]) cout << x << " ";
      //   cout << endl;
      // }
      // getchar();

      unsigned bound = 0;
      GRBLinExpr e = 0;
      for (unsigned r = 1; r < v_res.size(); ++r) {
        for (auto i1 : v_res[r-1]) {
          for (auto i2 : v_res[r-1]) e += P[1][i1][KPerm[i2]];
          bound += 1;
          // bool test = false;
          // for (auto i2 : v_res[r]) if (KPerm[i1] == i2) test = true;
          // if (!test) cout << "weird" << endl;
        }
      }
      model.addConstr(e <= bound - 1);

    }
  }
}

bool modelAES128(int R, vector<int> & KPerm, int nrSboxesWanted, int active, GRBEnv & env) {

    GRBModel model = GRBModel(env);

    //model.set(GRB_IntParam_PoolSolutions, 2000000);
    //model.set(GRB_DoubleParam_PoolGap, 0.001);
    //model.set(GRB_IntParam_PoolSearchMode, 2);

    // X = 3r + 0
    // S(X) = -(3*r + 0)
    // Y (after MC) = 3*r + 1
    // K = 3*r + 2

    // definition of variables
    vector<GRBVar> dX ((3*R + 1)*16); // before ARK

    for (unsigned i = 0; i < (3*R+1)*16; ++i) {
        dX[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }

    for (unsigned i = 0; i < 3*16; ++i) {
      model.addConstr(dX[i] == 0);
    }


    // Contraints for permutation-based key schedule
    // number of active bytes per key must be the same
    {
      GRBLinExpr e = 0;
      for (unsigned i = 0; i < 16; ++i) e += dX[16*(3*1 + 2) + i];
      for (unsigned r = 2; r < R; ++r) {
        GRBLinExpr e1 = 0;
        for (unsigned i = 0; i < 16; ++i) e1 += dX[16*(3*r + 2) + i];
        model.addConstr(e1 == e);
      }
    }

    for (unsigned i = 0; i < 16; ++i) {
      unsigned x = i;
      unsigned r = 1;
      while (r < R && x != 16) {++r; x = KPerm[x];}
      if (r == R) { // OK
        x = i;
        for (r = 2; r < R; ++r) {
            model.addConstr(dX[16*(3*r + 2) + KPerm[x]] == dX[16*(3*(r-1) + 2) + x]);
            x = KPerm[x];
        }
      }
      else { // NOK
        model.addConstr(dX[16*(3*1 + 2) + i] == 0);
      }
    }

    for (unsigned r = 1; r < R; ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 1)  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 0) + ((c + i)%4) + 4*i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);
      }
    }

    for (unsigned r = 1; r < R ; ++r) {
      //ARK (sX[r+1], kX[r+1] and zX[r+1])
      for (unsigned i = 0; i < 16; ++i) {
        model.addConstr(1-dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + 1-dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + 1-dX[16*(3*(r+1))  + i] >= 1);
      }
    }

    GRBLinExpr obj = 0;
    for (unsigned r = 1; r <= R; ++r) {
      for (unsigned i = 0; i < 16; ++i) obj += dX[16*(3*r) + i];
    }

    // path should involve the new guess
    {
      GRBLinExpr e = 0;
      for (unsigned r = 1; r < R; ++r) e += dX[16*(3*r + 2) + active];
      model.addConstr(e >= 1);
    }


    model.addConstr(obj >= 1);
    model.addConstr(obj <= nrSboxesWanted-1);

    auto mat = AES128eqs(R, KPerm.data());

    mycallback cb ((3*R+1)*16, dX.data(), mat);
    model.setCallback(&cb);
    model.set(GRB_IntParam_OutputFlag , 0);
    model.set(GRB_IntParam_LazyConstraints , 1);

    //printPerm(KPerm);

    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);

    return (nSolutions != 0);
}

int PermComplete(unsigned *P){
    for(int i = 0; i<16; i++){
        if(P[i] == 16) return 0;
    }
    return 1;
}

unsigned firstAvailableElm(unsigned *P){
  /*  unsigned min = 16;
    for(int i = 0; i<16; i++){
        if((P[i] != 16) && (P[i] < min)) min = P[i];
    }
    return min;*/
}

void searchRec(unsigned R, int nrSboxesWanted, vector<int> & KPerm, vector<int> & images, unsigned n_images, unsigned n_forbidden, int x, unsigned length, unsigned min_length, GRBEnv & env) {

   if (KPerm[x] != 16) {
       if (length < min_length || (n_images > 0 && n_images < length)) return;
       // toujours évaluer le cycle
       if (modelAES128(R, KPerm, nrSboxesWanted, x, env)) return;

       // regarder si l'on peut commencer un autre cycle ou si on a tout fixé
       if (n_images == 0) {
         cout << "une permutation choisie: " << endl;
         for (const auto & y : KPerm) cout << y << " ";
         cout << endl;
         getchar();
       }
       else {
         if (n_forbidden == 0) {
           searchRec(R, nrSboxesWanted, KPerm, images, n_images, 0, images[0], 0, length, env);
         }
         else if (n_images >= length + 1) searchRec(R, nrSboxesWanted, KPerm, images, n_images, 0, images[0], 0, length+1, env);

         if (n_images >= 2*length + 1) {
           unsigned n_f = 1;
           if (n_images >= 2*length + 2) {
             while (n_f < n_forbidden && n_images - n_f >= length+1) {
               searchRec(R, nrSboxesWanted, KPerm, images, n_images, n_f, images[n_f], 0, length + 1, env);
               ++n_f;
             }
           }
           n_f = n_forbidden;
           while (n_images - n_f >= length) {
             searchRec(R, nrSboxesWanted, KPerm, images, n_images, n_f, images[n_f], 0, length, env);
             ++n_f;
           }
         }

       }
   }
   else {
       // regarder si l'on peut évaluer le début du cycle (si length est suffisamment grand)
       if ((length >= R-1) && modelAES128(R, KPerm, nrSboxesWanted, x, env)) return;
       // continuer le début du cycle
       for (unsigned i = n_forbidden; i < n_images; ++i) {
         KPerm[x] = images[i];
         images[i] = images[n_images - 1];
         searchRec(R, nrSboxesWanted, KPerm, images, n_images-1, n_forbidden, KPerm[x], length+1, min_length, env);
         images[i] = KPerm[x];
       }
       KPerm[x] = 16;
   }
}

void search(unsigned R, int nrSboxesWanted) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  vector<int> KPerm(16, 16);
  vector<int> images(16);
  for (int i = 0; i < 16; ++i) images[i] = i;
  searchRec(R, nrSboxesWanted, KPerm, images, 16, 0, 0, 0, 0, env);
  searchRec(R, nrSboxesWanted, KPerm, images, 16, 4, 4, 0, 0, env);
  searchRec(R, nrSboxesWanted, KPerm, images, 16, 8, 8, 0, 0, env);
  searchRec(R, nrSboxesWanted, KPerm, images, 16, 12, 12, 0, 0, env);
}

void searchCyclesRec(unsigned R, int nrSboxesWanted, vector<int> & KPerm, vector<int> & images, unsigned n_images, int x, unsigned length, unsigned length_obj, GRBEnv & env, bool & found) {
   if (KPerm[x] != 16) {
       if (length != length_obj) return;
       // toujours évaluer le cycle
       if (modelAES128(R, KPerm, nrSboxesWanted, x, env)) return ;

       int y = KPerm[x];
       cout << "\r" << x << " --> ";
       while (y != x) {
         cout << y << " --> ";
         y = KPerm[y];
       }
       cout << x << endl;

       found = true;

       if (length == 16) {cout << "permutation!!" << endl; getchar();}
   }
   else {
      if (length == length_obj) return;
       // regarder si l'on peut évaluer le début du cycle (si length est suffisamment grand)
       if ((length >= R-1) && modelAES128(R, KPerm, nrSboxesWanted, x, env)) return;
       // continuer le début du cycle
       for (unsigned i = 0; i < n_images; ++i) {
         KPerm[x] = images[i];
         images[i] = images[n_images - 1];
         searchCyclesRec(R, nrSboxesWanted, KPerm, images, n_images-1, KPerm[x], length+1, length_obj, env, found);
         images[i] = KPerm[x];
       }
       KPerm[x] = 16;
   }
}

void searchCycles(unsigned R, int nrSboxesWanted) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  vector<int> KPerm(16, 16);
  vector<int> images(16);
  for (int i = 0; i < 16; ++i) images[i] = 15-i;
  vector<bool> found (17, false);
  vector<bool> possible (17, false);
  possible[0] = true;
  for (unsigned l = 1; l <= 16; ++l) {
    bool found_cycle = false;
    if (l <= 8 || possible[16-l]) {
      searchCyclesRec(R, nrSboxesWanted, KPerm, images, 16, 0, 0, l, env, found_cycle);
      if (l <= 12) searchCyclesRec(R, nrSboxesWanted, KPerm, images, 12, 4, 0, l, env, found_cycle);
      if (l <= 8) searchCyclesRec(R, nrSboxesWanted, KPerm, images, 8, 8, 0, l, env, found_cycle);
      if (l <= 4) searchCyclesRec(R, nrSboxesWanted, KPerm, images, 4, 12, 0, l, env, found_cycle);
    }
    found[l] = found_cycle;
    if (found_cycle) {
      for (unsigned i = l; i <= 16; i += l) {
        for (unsigned j = 0; j <= 16; j += 1) {
          if (possible[j] && j + i <= 16) possible[j+i] = true;
        }
      }
    }
    cout << "cycles " << l << ": " << (found_cycle ? "possible" : "impossible") << endl;
  }

}



int main(int argc, char const *argv[]) {
    //searchMILP(vector<pair<unsigned, unsigned>>({make_pair(9, 22)}), 32);

    //testBoomerang256(11);
    
    testPerm128(stoi(argv[1]));

    vector<pair<unsigned, unsigned>> constraints;
    unsigned size_perm = stoi(argv[1]);
    for (unsigned i = 2; i < argc; i += 2) constraints.emplace_back(stoi(argv[i]), stoi(argv[i+1]));
    searchMILP(constraints, size_perm);




    return 0;
}
