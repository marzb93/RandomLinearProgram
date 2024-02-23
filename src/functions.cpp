#include <cmath>
#include <vector>
#include <random>

#include "../include/functions.hpp"

using namespace std;

double calculateMagnitude(const vector<double>& vec) {
    double sumOfSquares = 0.0;
    for (double value : vec) {
        sumOfSquares += value * value;
    }
    return sqrt(sumOfSquares);
}
void normalizeVector(vector<double>& vec) {
    double magnitude = calculateMagnitude(vec);
    for (double& value : vec) {
        value /= magnitude;
    }
}
double sampleMean(const vector<double>& vec) {
    double sum = 0.0;
    for (double value : vec) {
        sum += value;
    }
    return sum / vec.size();
}

double sampleStdDev(const vector<double>& vec) {
    double mean = sampleMean(vec);
    double sum = 0.0;
    for (double value : vec) {
        sum += (value - mean) * (value - mean);
    }
    return sqrt(sum / (vec.size()-1) );
}
vector<double> generate_normal_cost(const int n){
    vector<double> cost(n);
    //const int seed = 123;
    std::random_device rd;
    std::mt19937 gen(rd());
    normal_distribution<double> ndist(0,1);
    for (int i =0; i<n; i++){
        cost[i] = ndist(gen);
    }
    normalizeVector(cost);
    return cost;
}
vector<double> generate_Rademacher_cost(const int n){
    vector<double> cost(n);
    //const int seed = 123;
    std::random_device rd;
    std::mt19937 gen(rd());
    uniform_int_distribution<int> u_int_dist(0,1);
    for (int i =0; i<n; i++){
        cost[i] = u_int_dist(gen);
    }
    normalizeVector(cost);
    return cost;
}
vector<double> generate_compressible_cost(const int n, const int k){
    vector<double> cost(n);
  for (int i =0; i<n; i++){
    if (i<k){
        cost[i] = static_cast<double>(1/sqrt(k));
    }
    else{
      cost[i] = 0;
    }
    }
    return cost;
}
vector<vector<double>> generate_normal_A(const int m, const int n){
    vector<vector<double>> A(m, vector<double>(n) ); // matrix of coefficients A
    //const int seed = 123;
    std::random_device rd;
    std::mt19937 gen(rd());
    normal_distribution<double> ndist(0,1);
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            A[i][j] = ndist(gen);
        }
    } 
    return A;
}
vector<vector<double>> generate_Rademacher_A(const int m, const int n){
  vector<vector<double>> A(m, vector<double>(n) ); // matrix of coefficients A
  std::random_device rd;
  std::mt19937 gen(rd()); 
  uniform_int_distribution<int> u_int_dist(0, 1);
  int randomNumber;
  for (int i=0; i<m; i++){
      for (int j=0; j<n; j++){
          randomNumber = u_int_dist(gen);
          if (randomNumber == 0){
              randomNumber = -1;
          }
          A[i][j] = randomNumber;
      }
  }
  return A;
}
vector<vector<double>> generate_uniform_A(const int m, const int n){
  vector<vector<double>> A(m, vector<double>(n) ); // matrix of coefficients A
  std::random_device rd;
  std::mt19937 gen(rd()); 
  normal_distribution<double> normal_dist(0.0, 1.0);
  const double min_val = -1.0;
  const double max_val = 1.0;
  for (int i=0; i<m; i++){
      for (int j=0; j<n; j++){
          double sample = normal_dist(gen);
          double scaled_value = (sample / 3.0); // Scale by sqrt(3) for variance 1
          scaled_value *= (max_val - min_val); // Scale to the desired range
          scaled_value += (max_val + min_val) / 2.0; // Shift mean to 0   
          A[i][j] = scaled_value;
      }
  }
  return A;
}
vector<vector<double>> generate_Bern_Gauss_A(const int m, const int n){
  vector<vector<double>> A(m, vector<double>(n) ); // matrix of coefficients A
  std::random_device rd;
  std::mt19937 gen(rd()); 
  // Bernoulli distribution with p = 0.5
  bernoulli_distribution bernoulliDist(0.5);
  // Gaussian distribution with variance 2
  normal_distribution<double> gaussianDist(0.0, std::sqrt(2.0));
  for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
          int bernoulliRandom = bernoulliDist(gen);
          double gaussianRandom = gaussianDist(gen);
          A[i][j] = bernoulliRandom * gaussianRandom;
      }
  }
  return A;
}
