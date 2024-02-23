#ifndef functions_HPP
#define functions_HPP

#include <vector>


using namespace std;

double calculateMagnitude(const vector<double>& vec);

void normalizeVector(vector<double>& vec);

double sampleMean(const vector<double>& vec);

double sampleStdDev(const vector<double>& vec);

vector<double> generate_normal_cost(const int n);

vector<double> generate_Rademacher_cost(const int n);

vector<double> generate_compressible_cost(const int n, const int k);

vector<vector<double>> generate_normal_A(const int m, const int n);

vector<vector<double>> generate_Rademacher_A(const int m, const int n);

vector<vector<double>> generate_uniform_A(const int m, const int n);

vector<vector<double>> generate_Bern_Gauss_A(const int m, const int n);

#endif
