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

class DC:public DATA{
public:
    static int gc;
    static V gclus;
    static MatrixXd gm;
    static size_t PT;
    static VectorXd p;
    
    DC(MatrixXd&, VectorXi&, VectorXd&,vector<int>&,double [],MatrixXi&);
    static inline double cb(double& alpha, double& beta, const VectorXd& rs, const VectorXd& cs){
        double res = 0;
        for(int i = 0; i < rs.size(); ++i){
            //        res+=R::lbeta(alpha+rs(i),beta+cs(i))-R::lbeta(alpha,beta);
            res += lgamma(alpha + rs(i)) + lgamma(beta + cs(i)) - lgamma(alpha + rs(i) + beta + cs(i)) - lgamma(alpha) - lgamma(beta) + lgamma(alpha + beta);
            
        }
        return res;
    }
    
    static inline double cb(double& alpha, double& beta, const double& rs,const double& cs){
      return lgamma(alpha + rs)+lgamma(beta + cs)-lgamma(alpha + rs + beta + cs) - lgamma(alpha) - lgamma(beta) + lgamma(alpha + beta);
        
    }
    static inline double LL(double hp[]){
        double res=0;
        double* thp = new double[gc+1];
        for(int i = 0; i < gc+1; ++i){
            thp[i] = exp(hp[i]);
            if(thp[i] > 1000)
                thp[i] = 1000;
        }
        size_t n = pat.rows();
        
        for(int i = 1;i < n;++i){
            int sub_k = pat.row(i).maxCoeff();
            MatrixXd tmp(G,sub_k);
            MatrixXd tmp_r(G,sub_k);
            tmp.fill(0);
            tmp_r.fill(0);
            for(int j = 0; j < sub_k; ++j){
                for(int t = 0; t < K; ++t){
                    if(pat(i,t) == j+1){
                        tmp.col(j) += r_d[t].rowwise().sum();
                        tmp_r.col(j) += r_r[t].rowwise().sum();
                    }
                }
            }
            
            for(int t1 = 0;t1<G;++t1){
                res += cb(thp[0], thp[gclus[t1]], tmp_r.row(t1), tmp.row(t1)) * gm(t1,i);
            }


            
        }

        for(int t1 = 0;t1<G;++t1){
            res += cb(thp[0], thp[gclus[t1]], r.row(t1), data.row(t1)) * gm(t1,0);
        }
        delete[] thp;
        return -1*res;
    }

    
    MatrixXd cal_gm(double []);
    
    MatrixXd cal_delta(MatrixXd&);
    
    ~DC();
};






#endif /* derived_hpp */
