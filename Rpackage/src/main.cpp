//
//  main.cpp
//  EBSeq
//
//  Created by MaXiuyu on 5/16/17.
//  Copyright © 2017 MaXiuyu. All rights reserved.
//
#include <Rcpp.h>
#include <RcppEigen.h>
#include <map>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <sys/types.h>
#include <sys/sysctl.h>
#include "DATA.hpp"
#include "derived.hpp"
#include "asa047.hpp"
#include "kit.h"

using namespace std;
using namespace Rcpp;
using namespace Eigen;
typedef vector<int> V;
















// MatrixXd DATA::data;

// vector<MatrixXd> DATA::r_d;

// vector<MatrixXd> DATA::r_r;

// size_t DATA::G;

// size_t DATA::K;

// MatrixXd DATA::r;

// VectorXd DATA::q;

// VectorXd DC::p;;

// int DC::gc;

// size_t DC::PT;

// vector<int> DC::gclus;

// MatrixXi DATA::pat;

// MatrixXd DC::gm;


Rcpp::IntegerVector MCP(Rcpp::IntegerVector X, double MASS, Rcpp::NumericVector PARAM);
RcppExport SEXP MCP(SEXP X, SEXP MASS, SEXP PARAM) {
    BEGIN_RCPP
    IntegerVector XX(X);
    double mass=Rcpp::as<double>(MASS);
    NumericVector param(PARAM);
    int len = XX.size();
    double shape = param[0];
    double scale = param[1];
    double coeff1 = lgamma(shape);
    double coeff2 = shape * log(1/scale);
    double med = coeff2-coeff1;
    double iscale = 1/scale;
    
    vector<int> dat;
    
    dat.assign(XX.begin(),XX.end());
    vector<size_t> ascend_indexes(len);
    iota(ascend_indexes.begin(), ascend_indexes.end(), 0);
    sort(ascend_indexes.begin(), ascend_indexes.end(), [&](size_t it1,size_t it2){return dat[it1] < dat[it2];});
    V ascend(dat);
    sort(ascend.begin(),ascend.end());
    
    
    
    
    //IMP to store incomplete modal partitions
    vector<component> IMP(len);
    IMP[0].logdensity = lgamma(ascend[0] + shape) + med -(ascend[0] + shape) * log(1 + iscale);
    IMP[0].logprior = log(mass);
    IMP[0].tailsum = ascend[0];
    IMP[0].tailsize = 1;
    IMP[0].partition.push_back(0);
    
    
    for(int i = 1; i < len; ++i){
        //        temporary container to store subset partition situations
        
        vector<component> tmp(i + 1);
        for(int j = 0;j < i; ++j){
            if(j != 0){
                tmp[j].tailsum = tmp[j - 1].tailsum + ascend[i - j];
            }
            else{
                tmp[0].tailsum = ascend[i - j];
            }
            tmp[j].tailsize = j + 1;
            tmp[j].logdensity = IMP[i -j -1].logdensity + lgamma(tmp[j].tailsum+shape) + med -(tmp[j].tailsum+shape) * log(tmp[j].tailsize+iscale);
            tmp[j].logprior=IMP[i-j-1].logprior+lgamma(j+1)+log(mass);
            
            
            
            tmp[j].partition.assign(IMP[i - j - 1].partition.begin(),IMP[i -j - 1].partition.end() );
            for(int k = 0;k < j + 1; ++k){
                tmp[j].partition.push_back(IMP[i - j - 1].partition[i - j - 1] + 1);
            }
            
            
            
        }
        tmp[i].tailsum = tmp[i - 1].tailsum + ascend[0];
        tmp[i].tailsize = i + 1;
        tmp[i].logdensity = lgamma(tmp[i].tailsum+shape)+med-(tmp[i].tailsum+shape) * log(tmp[i].tailsize+iscale)
        ;
        tmp[i].logprior=lgamma(i + 1) + log(mass);
        tmp[i].partition.assign(i + 1,0);
        
        
        IMP[i] = *max_element(tmp.begin(), tmp.end(),myfn);
        
    }
    
    
    
    
    
    for(int i = 0; i < len; ++i){
        XX[ascend_indexes[i]] = IMP[len-1].partition[i];
    }
    
    return XX;
    
    END_RCPP
}




int isref(Rcpp::IntegerVector X, Rcpp::IntegerVector Y);
RcppExport SEXP isref(SEXP X, SEXP Y){
    BEGIN_RCPP
    IntegerVector XX(X);
    IntegerVector YY(Y);
    const int n=YY.size();
    V A(n),B(n);
    
    copy(XX.begin(),XX.end(),A.begin());
    copy(YY.begin(),YY.end(),B.begin());
    
    
    
    int k_1 = *max_element(A.begin(),A.end());
    
    vector<int> C(n);
    for(int i = 0; i < n;++i)
    C[i] = A[i] - B[i];
    vector<int> dif(k_1, -k_1);
    for(int i = 0; i < n;++i){
        if(dif[A[i] - 1] == -k_1)
        dif[A[i]-1] = C[i];
        if(dif[A[i]-1] != C[i])
        return wrap(0);
    }
    return wrap(1);
    
    END_RCPP
    
    
}



Rcpp::List new_D(Rcpp::IntegerVector X, Rcpp::IntegerMatrix part, Rcpp::IntegerVector Y);
RcppExport SEXP new_D(SEXP X, SEXP part, SEXP Y){
    BEGIN_RCPP
    IntegerMatrix PP(part);
    MatrixXi Part(PP.rows(),PP.cols());
    copy(PP.begin(),PP.end(),Part.data());
    
    IntegerVector XX(X);
    vector<int> cd(XX.length());
    copy(XX.begin(),XX.end(),cd.begin());
    
    IntegerVector YY(Y);
    vector<int> de_i(YY.length());
    copy(YY.begin(),YY.end(),de_i.begin());
    
    MatrixXi nD(de_i.size(),cd.size());
    MatrixXd res(cd.size(),cd.size());
    nD.fill(0);
    res.fill(0);
    int K=*max_element(cd.begin(),cd.end());
    vector<vector<int> > sub_i;
    for(int i=0;i<K;i++)
        sub_i.push_back(which(cd,i+1));
    for(int i=0;i<de_i.size();i++){
        int tmp=Part.row(de_i[i]-1).maxCoeff();
        for(int j=1;j<=tmp;j++){
            for(int t=0;t<K;t++){
                if(Part(de_i[i]-1,t)==j){
                    for(vector<int>::iterator it=sub_i[t].begin();it!=sub_i[t].end();it++)
                        nD(i,*it)=j;
                }
            }
        }
    }
    for(int i=0;i<cd.size();i++){
        for(int j=i+1;j<cd.size();j++){
            double count=0;
            for(int t=0;t<de_i.size();t++){
                if(nD(t,i)!=nD(t,j))
                    count++;
            }
            res(i,j)=count/de_i.size();
        }
    }
    return Rcpp::List::create(Named("nD")=nD,Named("new_dist")=res,Named("part")=Part);
    
    END_RCPP
}



Rcpp::IntegerMatrix g_ref(Rcpp::IntegerMatrix Posp);
RcppExport SEXP g_ref(SEXP Posp){
    BEGIN_RCPP
    IntegerMatrix PP(Posp);
    const int ng=PP.rows();
    const int nc=PP.cols();
    MatrixXi Part(ng,nc);
    MatrixXi res(ng,ng);
    copy(PP.begin(),PP.end(),Part.data());
    VectorXi K1(ng);
    K1 = Part.rowwise().maxCoeff();
    vector<int> C(nc);
    for(int i = 0; i < ng; i++){
        for(int j = 0; j < ng; j++){
            vector<int> dif(K1(i), K1(i));
            for(int t = 0; t < nc; t++)
                C[t] = Part(i,t) - Part(j,t);
            res(i,j) = 1;
            for(int t = 0; t < nc; t++){
                if(dif[Part(i,t) - 1] == K1(i))
                    dif[Part(i,t) - 1] = C[t];
                if(dif[Part(i,t) - 1] != C[t]){
                    res(i,j) = 0;
                    break;
                }
            }
        }
    }
    return wrap(res);
    END_RCPP
}

Rcpp::List pat(int K);
RcppExport SEXP pat(SEXP K){
    BEGIN_RCPP
    int k=as<int>(K);
    vector<vector<int> > res=partition(k);
    MatrixXi res_m(res.size(),k);
    for(int i=0;i<res.size();i++)
        for(int j=0;j<k;j++)
            res_m(i,j)=res[i][j];
    return List::create(Named("part")=res_m);
    END_RCPP
}


Rcpp::List EBS(Rcpp::NumericMatrix X, Rcpp::IntegerVector Y, Rcpp::IntegerVector Z, Rcpp::NumericVector W, int iter, Rcpp::NumericVector hyper, Rcpp::IntegerMatrix part);
RcppExport SEXP EBS(SEXP X, SEXP Y, SEXP Z, SEXP W, SEXP iter, SEXP hyper,SEXP part) {
    
    BEGIN_RCPP
    int itr=as<int>(iter);
    NumericMatrix XX(X);
    IntegerMatrix PP(part);
    //    cell clus
    IntegerVector YY(Y);
    //    gene clus
    IntegerVector ZZ(Z);
    NumericVector WW(W);
    NumericVector TT(hyper);
    
    const int ng=XX.rows();
    const int nc=XX.cols();
    
    MatrixXd data(ng,nc);
    
    copy(XX.begin(),XX.end(),data.data());
    
    MatrixXi Part(PP.rows(),PP.cols());
    
    copy(PP.begin(),PP.end(),Part.data());
    
    VectorXi conditions(nc);

    copy(YY.begin(),YY.end(),conditions.data());
    
    vector<int> gclus(ng);
    
    copy(ZZ.begin(),ZZ.end(),gclus.begin());
    
    int gcc=*max_element(gclus.begin(),gclus.end());
    
    VectorXd sf(nc);
    
    copy(WW.begin(),WW.end(),sf.data());
    
    VectorXd hyp(gcc+1);
    
    copy(TT.begin(),TT.end(),hyp.data());
    
    //initialize parameters
    double hp[gcc+1];
//        alpha
    hp[0] = hyp[0];
//        beta
    for(int i = 0; i < gcc;++i)
        hp[i + 1] = hyp[i + 1];
    
    double stepsize1 = 1e-6;
    
    double stepsize2 = 1e-5;
    
//        DATA init(data,conditions,sf);    
    DC init(data, conditions, sf, gclus, Part);
    
    cout<<"first init"<<endl;

    vector<MatrixXd> A = init.cal_gm(hp);
    
    init.gm = init.cal_delta(A[0]);
    
    double dAlpha = init.cal_drv(A[1]);
    
    double dBeta = init.cal_drv(A[2]);
    
    cout<<"cal_gm done"<<endl;
    
    VectorXd tm_p(init.PT);
    
    if(hp[0] + stepsize1 * dAlpha > 0)
        hp[0] = hp[0] + stepsize1 * dAlpha;
    
    if(hp[1] + stepsize2 * dBeta > 0)
        hp[1] = hp[1] + stepsize2 * dBeta;
    
    
    double tt = init.gm.sum();
    
    tm_p = init.gm.colwise().sum()/tt;
    
    init.p = tm_p;
    
    double diff = 1;
    
    int o = 0;
    
    while(diff > 0.01 && o < itr){
        stepsize1 /= 2;
        stepsize2 /= 2;
        
        A = init.cal_gm(hp);
        init.gm = init.cal_delta(A[0]);
        
        double dAlpha = init.cal_drv(A[1]);
    
        double dBeta = init.cal_drv(A[2]);
        
        if(hp[0] + stepsize1 * dAlpha > 0)
            hp[0] = hp[0] + stepsize1 * dAlpha;
    
        if(hp[1] + stepsize2 * dBeta > 0)
            hp[1] = hp[1] + stepsize2 * dBeta;
        
        tt = init.gm.sum();
        tm_p = init.gm.colwise().sum()/tt;
        diff = abs((tm_p - init.p).maxCoeff());
        init.p = tm_p;
        ++o;
    }
    
    cout<<"final done"<<endl;

    return Rcpp::List::create(Named("DEpattern")=init.gm,Named("r")=init.r,Named("q")=init.q,Named("p")=init.p,Named("iteration")=o,Named("Alpha") = hp[0], Named("Beta") = hp[1], Named("obj") = A[0] * init.p);
    
    END_RCPP
    
}

