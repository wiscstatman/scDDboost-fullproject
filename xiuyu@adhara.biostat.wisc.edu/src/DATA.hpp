#ifndef DATA_hpp
#define DATA_hpp
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <map>
using namespace std;
using namespace Eigen;
typedef vector<int> V;
#include "Header.hpp"

struct r_q{
    MatrixXd rrr;
    VectorXd qqq;
};

struct hyperparam{
    double alp;
    VectorXd bta;
};


class DATA{
public:
    
    static vector<MatrixXd> r_d;
    
    static vector<MatrixXd> r_r;
    
    static VectorXd c_1;
    
    static VectorXd r_1;
    
    static size_t G;
    
    static MatrixXd r;
    
    static size_t K;
    
    static MatrixXd data;
    
    static VectorXd q;
    
    static MatrixXi pat;
    
    MatrixXd d_s;
    
    MatrixXd r_s;
    
    vector<MatrixXd> reorg(MatrixXd&,  const VectorXi&);
    
    vector<MatrixXd> reorgr(const MatrixXd&, const VectorXi&);
    
    DATA(MatrixXd&, const VectorXi&, const VectorXd&, const MatrixXi&);
    
    r_q cal_r(MatrixXd&, const VectorXi&, const VectorXd&);
    
    MatrixXd cal_hp(void);
    
    void cal_sum(void);
    
    ~DATA();
};

#endif /* DATA_hpp */
