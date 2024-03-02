#ifndef GurobiModelWrapper_HPP
#define GurobiModelWrapper_HPP


#include "gurobi_c++.h"
#include <gurobi_c.h>
#include <string>
#include <vector>

using namespace std;

class rand_LP{

public: 
  GRBVar* x;
  GRBConstr* con;

public:  
  GRBModel model;
  //constructor to initialize Gurobi model
  rand_LP(GRBEnv env) : model(env)
  {
  };  

  void add_var_constr(vector<vector<double>>& A){
    int n = A[0].size();
    int m = A.size();
    x = new GRBVar[n];
    for (int i=0; i<n; i++){
        x[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_"+to_string(i) );
    }
    con = new GRBConstr[m];
    for (int i=0; i<m; i++){
      GRBLinExpr expr;
      for (int j=0; j<n; j++){
        expr += A[i][j] * x[j];
      }
      con[i] = model.addConstr(expr <= 1 , "con_" + to_string(i));
    }
  }
  void add_obj(vector<double>& cost){
    int n = cost.size();
    GRBLinExpr obj;
    for (int i=0; i<n; i++){
        obj += cost[i]*x[i];
    }
    model.setObjective(obj, GRB_MAXIMIZE);

  }
  void cleanup_var_constr(){
    delete [] x;
    delete [] con;
  }

};

#endif
