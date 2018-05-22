//
//  derived.cpp
//  test for map<matrixXd>
//
//  Created by MaXiuyu on 6/9/17.
//  Copyright © 2017 MaXiuyu. All rights reserved.
//

#include "derived.hpp"

DC::DC(MatrixXd& data, VectorXi& conditions, VectorXd& sf, vector<int>& g_clus, double xmin[], MatrixXi& PP):DATA(data, conditions, sf, PP){
    this->gclus = g_clus;
    this->gc = *max_element(g_clus.begin(), g_clus.end());
    this->PT = PP.rows();
    VectorXd p1(PT);
    p1.fill(1/double(PT));
    this->p=p1;
}



MatrixXd DC::cal_gm(double xmin[]){
    MatrixXd tmpp(G,PT);
    for(int i = 1; i < PT; ++i){
        int sub_k = pat.row(i).maxCoeff();
        MatrixXd tmp(G,sub_k);
        MatrixXd tmp_r(G,sub_k);
        tmp.fill(0);
        tmp_r.fill(0);
        for(int j = 0; j < sub_k; ++j){
            for(int t = 0; t<K; ++t){
                if(pat(i,t) == j+1){
                    tmp.col(j) += d_s.col(t);
                    tmp_r.col(j) += r_s.col(t);
                }
            }
        }

        for(int j = 0;j < G; ++j){
            tmpp(j,i) = cb(xmin[0],xmin[gclus[j]],tmp_r.row(j),tmp.row(j));
        }
    }
    
    for(int j = 0;j < G; ++j){
        tmpp(j,0) = cb(xmin[0],xmin[gclus[j]],r.row(j).sum(),data.row(j).sum());
    }
    return tmpp;
}

MatrixXd DC::cal_delta(MatrixXd &A){
    MatrixXd tmpp(G,PT);
    for(int j=0;j<G;++j){
        double M=A.row(j).maxCoeff();
        double total=0;
        for(int i=0;i<PT;++i){
            total+=exp(A(j,i)-M)*p(i);
        }
        for(int i=0;i<PT;++i){
            tmpp(j,i)=exp(A(j,i)-M)/total*p(i);
        }
        
    }
    return tmpp;

}


DC::~DC(){
    gm.resize(0,0);
    vector<int>().swap(gclus);
    p.resize(0);
}




