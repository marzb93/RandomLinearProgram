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

#include <Eigen/Dense>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "../include/GurobiModelWrapper.hpp"
#include "../include/functions.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  if (argc < 5){
    cout << "Usage: ./rlp m n sampleSize distributionType" << endl;
    return 1;
  }
  // for the assumption on the cost vector ==================================
  // if (argc < 6){
  //     cout << "Error: m, n, max_iter, type, K expected" << endl;
  //     return 1;
  // }
  // int n = stoi(argv[2]); // number of variables
  // int K = stoi(argv[6]); // Kth component of the cost vector are non zero
  // double c[n]; // cost vactor
  // for (int i =0; i<n; i++){
  // if (i<K){
  //     c[i] = static_cast<double>(1/sqrt(K));
  // }
  //      else{
  //   c[i] = 0;
  // }
  // }
  // ========================================================================
  const int m = stoi(argv[1]); // number of constraints
  const int n = stoi(argv[2]); // number of variables
  const int max_iter = stoi(argv[3]);  // sample size 
  string type = argv[4]; // distribution of matrix A: normal or pmone or uniform

  vector <double> c(n); // cost vactor, dynamic array
  vector<vector<double>> A(m, vector<double>(n)); // matrix of coefficients A, dynamic matrix
  // 
  double opt_value; // optimal objective value
  // Define accumulator set for mean and standard deviation of opt objective value
  using namespace boost::accumulators;
  accumulator_set<double, features<tag::mean, tag::variance>> obj_acc;

  vector <double> soln;
  // to compute angle between cost vector and optimal solution 
  double angle;
  vector <double> angles;

  // generate and normalize cost from i.i.d. Rademacher distribution
  for (int i =0; i<n; i++){
    c[i] = (static_cast<double>(rand()%2)-0.5)*2.0/sqrt(n);
  }
  // Eigen::Map<Eigen::VectorXd> c_eigen(c.data(),n); 

  GRBEnv env = GRBEnv(true);
  //env.set("logFile", "rlp.log");
  env.start();
  env.set(GRB_IntParam_LogToConsole, 0);



  // for loop for sampling the random LP  -----------------------------------------------------
  for (int iter=0; iter < max_iter; iter++){

    if (iter % 5 == 0){
      cout << "iter : "<< iter << endl;
    }

    // if matrix A has normal distribution 
    if (type == "normal") {
      A = generate_normal_A(m,n);
    }
    // if matrix A has distribution {-1,1}, i.i.d. Rademacher distribution
    else if (type == "pmone"){
      A = generate_Rademacher_A(m,n); 
    }
    // if matrix A has real uniform distribution
    else if (type == "uniform"){
      A = generate_uniform_A(m,n);
    }
    else if (type == "Bern_Gauss"){
      A = generate_Bern_Gauss_A(m,n);
    }
    else {
      cout << "Enter normal or pmone or uniform or Bern_Gauss for the 5th argument!" << endl;
    }

    rand_LP rlp = rand_LP(env);
    rlp.add_var_constr(A);
    rlp.add_obj(c);

    rlp.model.optimize();
    //rlp.model.write("rlp.lp");
    if (rlp.model.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
      opt_value = rlp.model.get(GRB_DoubleAttr_ObjVal);
      //cout << "obj value: " << opt_value << endl;
      obj_acc(opt_value);
      for (int i=0; i<n; i++){
        soln.push_back( rlp.x[i].get(GRB_DoubleAttr_X) );
      }
      // Eigen::Map<Eigen::VectorXd> soln_eigen(soln.data(), soln.size());
      // angle = std::acos(soln_eigen.dot(c_eigen) / (soln_eigen.norm() * c_eigen.norm())) * (180.0 / M_PI);
      // angles.push_back(angle);
      // // cout << angle << endl;
      // 
      // soln.clear();
      rlp.cleanup_var_constr();
    } else {
      //cout << "\n no solution found for i = " <<iter << endl;
      ;
    }
  }
  //for loop ends  -----------------------------------------------------------------------------
  //===============================================================================================  
  // Compute mean and standard deviation
  double mean = boost::accumulators::mean(obj_acc);
  cout << "mean : " << mean <<endl;
  double variance = boost::accumulators::variance(obj_acc);
  double std_dev = sqrt(variance);
  double l_b;
  double gap;
  double relative_gap;
  double sigma_sqrt_m;
  l_b = 1.0/(sqrt(2.0*log(m/n) ));
  // gap = abs(l_b - mean);
  relative_gap = 100* abs(l_b - mean)/l_b;
  // sigma_sqrt_m = StdDev * sqrt(m);

  //=========================================================================================
  // storing opt objectives for all samples in a file
  // ofstream objFile("obj.txt", ios::app); // 
  // objFile << fixed << setprecision(6);
  // objFile << "  distibution: " << type << "  --------------------------------" << endl;
  // objFile << m <<"  "<<n<< "  objectives : " << endl;
  // for (const auto& element : objectives) {
  // objFile << element << ", ";
  // }
  // objFile << endl;
  // objFile.close();
  // ============================================================================================
  //for the magnitude of the opt objective value table  
  ofstream obj_table("/home/marzieh/Desktop/research/random_lp/outputs/obj_table.txt", ios::app); 
  obj_table<< fixed << setprecision(6);
  obj_table<< "  distibution: " << type << "  ------------------------------" << endl;
  obj_table<< m <<" & "<< n <<" & "<<l_b<<" & "<<mean<<" & "<< relative_gap <<" \\\\"<<endl;
  obj_table.close();
  // ============================================================================================
  // for the Std Dev table 
  // ofstream outFile("out.txt", ios::app); // 
  // outFile << fixed << setprecision(6);
  // outFile << "  distibution: " << type << "  ------------------------------" << endl;
  // outFile<< m <<" & "<< n <<" & "<<l_b<<" & "<< mean << " & "<< relative_gap <<endl;
  //for std dev conjecture:
  // outFile<< m <<" & "<< n <<" & "<<l_b<< " & " << StdDev<<" & "<<sigma_sqrt_m<<" \\\\"  <<endl;
  // outFile.close();
  // ===========================================================================================
  // for assumption on the cost vector table
  // double mu_hat_1000_100 = 0.50626; 
  // double mu_gap;
  // mu_gap = 100* abs(mu_hat_1000_100 - mean)/mu_hat_1000_100;
  // ofstream outFile("cost_obj.txt", ios::app); 
  // outFile << fixed << setprecision(6);
  // outFile<<  K <<" & "<< mean << " & " << mu_gap << "\\\\"  <<endl;
  // outFile.close();
  // ============================================================================================
  //for the angle between the solution and the cost vector
  // double angle_mean = sampleMean(angles);
  // cout << "Angle between soln and cost: " << angle_mean << " degrees." << endl;
  // ofstream aFile("angle.txt", ios::app); // 
  // aFile << fixed << setprecision(4);
  // // aFile << "distibution: " << type << "  ------------------------------" << endl;
  // aFile<< m <<" & "<< n <<" & "<< angle_mean <<" \\\\"  <<endl;
  // aFile.close();
  //============================================================================================


  return 0;
}
