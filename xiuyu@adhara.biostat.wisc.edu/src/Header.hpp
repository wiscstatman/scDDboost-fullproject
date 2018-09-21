//
//  Header.hpp
//  test for map<matrixXd>
//
//  Created by MaXiuyu on 6/9/17.
//  Copyright © 2017 MaXiuyu. All rights reserved.
//

#ifndef Header_hpp
#define Header_hpp

#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include <RcppEigen.h>
#include <map>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
typedef vector<int> V;

V unique(V);
struct Sum {
    MatrixXd count_sum;
    MatrixXd r_sum;
};
vector<V> partition(const int& n);
vector<int> which (const vector<int>& x, const int& index);
vector<int> which (const VectorXi x, const int& index);

Sum rc(MatrixXd x,MatrixXd r, vector<V> pos);



#endif /* Header_hpp */
