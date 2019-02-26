//
//  derived.hpp
//  test for map<matrixXd>
//
//  Created by MaXiuyu on 6/9/17.
//  Copyright © 2017 MaXiuyu. All rights reserved.
//

#ifndef derived_hpp
#define derived_hpp

#include <stdio.h>
#include "DATA.hpp"
#include <boost/math/special_functions/digamma.hpp>

class DC:public DATA{
    
public:
    int gc;
    V gclus;
    MatrixXd gm;
    size_t PT;
    VectorXd p;
    std::unordered_map<size_t, Eigen::MatrixXd> converters;

    
public:
    
    DC(MatrixXd&, VectorXi&, VectorXd&, vector<int>&, MatrixXi&);
    
    //prior predictive function(PPF) on one group of subtypes
    VectorXd cb(double alpha, double beta, const MatrixXd& rs, const MatrixXd& cs);
        
    //derivative w.r.t. alpha and beta of PPF
    MatrixXd drv(double alpha, double beta, const MatrixXd& rs, const MatrixXd& cs);
    
    //derivate w.r.t. alpha and beta of log likelihood 
    double cal_drv(MatrixXd&);
    
    vector<MatrixXd> cal_gm(double []);
    
    MatrixXd cal_delta(MatrixXd&);
    
    ~DC();
};






#endif /* derived_hpp */
