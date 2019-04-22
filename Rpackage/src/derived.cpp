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


//prior predictive function(PPF) on one group of subtypes
VectorXd DC::cb(double alpha, VectorXd& beta, const MatrixXd& rs, const MatrixXd& cs)
{
    MatrixXd A = (rs.array() + alpha).matrix();
    
    MatrixXd B = cs.colwise() + beta;    
    
    MatrixXd res = ((A.unaryExpr<double(*)(double)>(&lgamma) + B.unaryExpr<double(*)(double)>(&lgamma) - (A + B).unaryExpr<double(*)(double)>(&lgamma)).array() - lgamma(alpha)).matrix();
    
    res =  res.colwise() - (beta.unaryExpr<double(*)(double)>(&lgamma) + (alpha + beta.array()).matrix().unaryExpr<double(*)(double)>(&lgamma));
    
    VectorXd RES = res.rowwise().sum();
                      
    return RES;
}



//derivative w.r.t. alpha and beta of PPF
MatrixXd DC::drv(double alpha, VectorXd& beta, const MatrixXd& rs, const MatrixXd& cs)
{
     MatrixXd A = (rs.array() + alpha).matrix();
    
     MatrixXd B = cs.colwise() + beta;    
    
     MatrixXd C = (A + B).unaryExpr<double(*)(double)>(&(boost::math::digamma));
    
     VectorXd D = (alpha + beta.array()).matrix().unaryExpr<double(*)(double)>(&(boost::math::digamma));
    
     MatrixXd resAlpha = ((A.unaryExpr<double(*)(double)>(&(boost::math::digamma)) - C).array() - boost::math::digamma(alpha)).matrix();
                          
     resAlpha = resAlpha.colwise() + D;

     MatrixXd resBeta = B.unaryExpr<double(*)(double)>(&(boost::math::digamma)) - C;
    
     resBeta = resBeta.colwise() - (beta.unaryExpr<double(*)(double)>(&(boost::math::digamma)) + D);
    
     MatrixXd res(G,2);
     
     res.col(0) = resAlpha.rowwise().sum();
     
     res.col(1) = resBeta.rowwise().sum();
    
     return res;
        
}
 

//derivative component w.r.t. alpha and beta of log likelihood
//log likelihood component
//optimizing via convert for iteration to matrix operation
vector<MatrixXd> DC::cal_gm(double alpha, VectorXd& beta){
    vector<MatrixXd> res(3);
    
    res[0].resize(G,PT);
    res[1].resize(G,PT);
    res[2].resize(G,PT);
    
    
//     MatrixXd tmpp(G,PT);
//     MatrixXd dAlpha(G,PT);
//     MatrixXd dBeta(G,PT);
    
    for(size_t i = 0; i < PT; i++){
        int sub_k = pat.row(i).maxCoeff();
        MatrixXd tmp(G,sub_k);
        MatrixXd tmp_r(G,sub_k);
        
        
        
        tmp = d_s * converters[i];
        tmp_r = r_s * converters[i];
        
        res[0].col(i) = cb(alpha, beta, tmp_r, tmp);
        
        MatrixXd DRV = drv(alpha, beta,  tmp_r, tmp);
        
        res[1].col(i) = DRV.col(0);
        res[2].col(i) = DRV.col(1);
    }
    
    
     
   
    return res;
}





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


void DC::go_drv(const MatrixXd& A, const MatrixXd& B, double& alpha, VectorXd& beta, double stepsize1, double stepsize2)
{
    
    double tmp1 = alpha + stepsize1 * (A * p).sum();
    VectorXd tmp2 = beta + stepsize2 * (B * p);
    for(int i = 0; i < G; i++)
    {
        if(tmp2[i] > 0)
        {
            beta[i] = tmp2[i];
        }
    } 
    if(tmp1 > 0)
    {
        alpha = tmp1;
    }
}
    



DC::~DC(){
    gm.resize(0,0);
    vector<int>().swap(gclus);
    p.resize(0);
}




