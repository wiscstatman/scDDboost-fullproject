//
//  derived.cpp
//  test for map<matrixXd>
//
//  Created by MaXiuyu on 6/9/17.
//  Copyright © 2017 MaXiuyu. All rights reserved.
//

#include "derived.hpp"

DC::DC(MatrixXd& data, VectorXi& conditions, VectorXd& sf, vector<int>& g_clus, MatrixXi& PP):DATA(data, conditions, sf, PP){
    //initial gene cluster label (currently treat each gene independently)
    this->gclus = g_clus;
    
    //number of gene clusters (currently same as number of genes)
    this->gc = *max_element(g_clus.begin(), g_clus.end());
    
    //number of partitions under K groups
    this->PT = PP.rows();
    
    //initial proportions parameters (equal proportion)
    VectorXd p1(PT);
    p1.fill(1/double(PT));
    this->p=p1;
    
    //initial hash map of converter matrix (convert individual sum of each group to aggregated sum over multiple groups)
    for(size_t i = 0; i < PT; i++){
        int sub_k = PP.row(i).maxCoeff();
        MatrixXd converter(K, sub_k);
        converter.fill(0);
        for(int j = 0; j < K; j++){
            converter(j, PP(i,j) - 1) = 1;
        }
        (this->converters)[i] = converter;
    }
    
}



//optimizing via convert for iteration to matrix operation
MatrixXd DC::cal_gm(double xmin[]){
    MatrixXd tmpp(G,PT);
    
    
    for(size_t i = 0; i < PT; i++){
        int sub_k = pat.row(i).maxCoeff();
        MatrixXd tmp(G,sub_k);
        MatrixXd tmp_r(G,sub_k);
        
        tmp = d_s * converters[i];
        tmp_r = r_s * converters[i];
        
        tmpp.col(i) = cb(xmin[0], xmin[1], tmp_r, tmp);
        
    }
     
   
    return tmpp;
}



// MatrixXd DC::cal_gm(double xmin[]){
//     MatrixXd tmpp(G,PT);
//     for(int i = 1; i < PT; ++i){
//         int sub_k = pat.row(i).maxCoeff();
//         MatrixXd tmp(G,sub_k);
//         MatrixXd tmp_r(G,sub_k);
//         tmp.fill(0);
//         tmp_r.fill(0);
//         for(int j = 0; j < sub_k; ++j){
//             for(int t = 0; t<K; ++t){
//                 if(pat(i,t) == j+1){
//                     tmp.col(j) += d_s.col(t);
//                     tmp_r.col(j) += r_s.col(t);
//                 }
//             }
//         }

//         for(int j = 0;j < G; ++j){
//             tmpp(j,i) = cb(xmin[0],xmin[gclus[j]],tmp_r.row(j),tmp.row(j));
//         }
//     }
    
//     for(int j = 0;j < G; ++j){
//         tmpp(j,0) = cb(xmin[0],xmin[gclus[j]],r.row(j).sum(),data.row(j).sum());
//     }
//     return tmpp;
// }


MatrixXd DC::cal_delta(MatrixXd &A){
    MatrixXd tmpp(G,PT);
    
    VectorXd M = A.rowwise().maxCoeff();
    
    tmpp = A.colwise() - M;
    
    tmpp = tmpp.unaryExpr<double(*)(double)>(& exp);
    
    
    VectorXd total = tmpp * p;   
    
    total = (1 / total.array()).matrix();
    
    //outer product of total and p
    MatrixXd div = total * p.transpose();  
    
    tmpp = (tmpp.array() * div.array()).matrix();
         
    return tmpp;

}



// MatrixXd DC::cal_delta(MatrixXd &A){
//     MatrixXd tmpp(G,PT);
//     for(int j=0;j<G;++j){
//         double M=A.row(j).maxCoeff();
//         double total=0;
//         for(int i=0;i<PT;++i){
//             total+=exp(A(j,i)-M)*p(i);
//         }
//         for(int i=0;i<PT;++i){
//             tmpp(j,i)=exp(A(j,i)-M)/total*p(i);
//         }
        
//     }
//     return tmpp;

// }


DC::~DC(){
    gm.resize(0,0);
    vector<int>().swap(gclus);
    p.resize(0);
}




