#include "gurobi_c++.h"
#include <gurobi_c.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cstdlib> // for calling python program

#include "../include/GurobiModelWrapper.hpp"
#include "../include/functions.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  if (argc < 6){
    cout << "Usage: ./rlp subsection m n sampleSize distributionType" << endl;
    return 1;
  }
  string subsection = argv[1]; // corresponds to subsections of numerical experiments
  const int m = stoi(argv[2]); // number of constraints
  const int n = stoi(argv[3]); // number of variables
  const int max_iter = stoi(argv[4]);  // sample size 
  string type = argv[5]; // distribution of matrix A: normal or pmone or uniform or Bern_Gauss

  vector <double> c(n); // cost vactor, dynamic array
  vector<vector<double>> A(m, vector<double>(n)); // matrix of coefficients A, dynamic matrix
  
  double obj; // optimal objective value
  vector <double> objectives; // for storing sample optimal objective values
  
  c = generate_Rademacher_cost(n);

  // for the assumption on the cost vector ==================================
  int k;
  if (subsection == "costAssumption"){ 
    k = stoi(argv[6]); // the number of non zeros in the cost vector
    c = generate_compressible_cost(n, k);
  } 
  //========================================================================


  GRBEnv env = GRBEnv(true);
  //env.set("logFile", "rlp.log");
  env.start();
  env.set(GRB_IntParam_LogToConsole, 0);


  // for loop for sampling the random LP  ---------------------------------------------------
  for (int iter=0; iter < max_iter; iter++){

    if (iter % 5 == 0){
      cout << "iter : "<< iter << endl;
    }

    // if matrix A has normal distribution 
    if (type == "normal") {
      A = generate_normal_A(m,n);
    }
    // if matrix A has distribution {-1,1}, i.i.d. Rademacher distribution
    else if (type == "Redemacher"){
      A = generate_Rademacher_A(m,n); 
    }
    // if matrix A has real uniform distribution
    else if (type == "uniform"){
      A = generate_uniform_A(m,n);
    }
    // if matrix A is product of Bernoulli(1/2) and N(0,2)
    else if (type == "Bern_Gauss"){
      A = generate_Bern_Gauss_A(m,n);
    }
    else {
      cout << "Enter normal or Redemacher or uniform or Bern_Gauss for 5th argument" << endl;
    }

    rand_LP rlp = rand_LP(env);
    rlp.add_var_constr(A);
    rlp.add_obj(c);

    rlp.model.optimize();
    //rlp.model.write("rlp.lp");
    if (rlp.model.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
      obj = rlp.model.get(GRB_DoubleAttr_ObjVal);
      objectives.push_back( obj);

      rlp.cleanup_var_constr();
    } else {
      //cout << "\n no solution found for i = " <<iter << endl;
      ;
    }
  }
  // for loop ends  ------------------------------------------------------------------------

  double mean = sampleMean(objectives);
  double StdDev = sampleStdDev(objectives);
  double b = 1.0/(sqrt(2.0*log(m/n) ));
  double relative_gap = 100* abs(b - mean)/b;
  double sigma_sqrt_m = StdDev * sqrt(m);

   

  if (subsection == "objMagnitude"){  //for the magnitude of the opt objective value table 
      ofstream objFile("objMagnitude.txt", ios::app); 
      objFile<< fixed << setprecision(6);
      objFile<< m <<" & "<< n <<" & "<<b<<" & "<<mean<<" & "<< relative_gap <<" \\\\"<<endl;
      objFile.close();
  } else if (subsection == "costAssumption"){ // for assumption on the cost vector table
      double mu_hat_1000_100 = 0.50626; 
      double mu_gap;
      mu_gap = 100* abs(mu_hat_1000_100 - mean)/mu_hat_1000_100;
      ofstream outFile("costAssumption.txt", ios::app); 
      outFile << fixed << setprecision(6);
      outFile<<  k <<" & "<< mean << " & " << mu_gap << "\\\\"  <<endl;
      outFile.close();
  } else if (subsection == "objLimitingDistribution"){ //for limiting distribution of objectives 
      ofstream objFile("objLimitingDistribution.txt");
      objFile << fixed << setprecision(6);
      double last_element = objectives.back();
      objectives.pop_back();
      for (const auto& element : objectives) {
      objFile << element << ", ";
    }
      objFile << last_element << endl;
      objFile.close();
      system("python3 ../src/hypothesis_testing.py");
  } else if (subsection == "objStdDev"){ // for the Std Dev table 
      ofstream outFile("objStdDev.txt", ios::app); 
      outFile << fixed << setprecision(6);
      outFile<< m <<" & "<< n <<" & "<<b<< " & " << StdDev<<" & "<<sigma_sqrt_m<<" \\\\"  <<endl;
      outFile.close();
  }

  return 0;
}
