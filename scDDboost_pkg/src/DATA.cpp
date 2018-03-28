//
//  DATA.cpp
//  EBSeq
//
//  Created by MaXiuyu on 5/22/17.
//  Copyright © 2017 MaXiuyu. All rights reserved.
//

#include "DATA.hpp"

vector<MatrixXd> DATA::reorg(MatrixXd& data, const VectorXi& conditions){
    vector<MatrixXd> tmp;
    int ngr=conditions.maxCoeff();
    for(int i=0;i<ngr;++i){
        vector<int> s=which(conditions,i+1);
        MatrixXd tc(G,s.size());
        for(int j=0;j<s.size();++j){
            tc.col(j)=data.col(s[j]);
        }
        tmp.push_back(tc);
    
    }
    return tmp;
}

vector<MatrixXd> DATA::reorgr(const MatrixXd& r, const VectorXi& conditions){
    vector<MatrixXd> tmp;
    int ngr=conditions.maxCoeff();
    for(int i=0;i<ngr;++i){
        vector<int> s=which(conditions,i+1);
        MatrixXd tr(G,s.size());
        for(int j=0;j<s.size();++j){
            tr.col(j)=r.col(s[j]);
        }
        tmp.push_back(tr);
    }
    return tmp;
}


r_q DATA::cal_r(MatrixXd& data,const VectorXi& conditions, const VectorXd& sf){
    VectorXd whole_mean(G);
    size_t nc = data.cols();
    //    whole_mean = data.rowwise().mean();
    
    MatrixXd data_dvd(G,nc);
    for(int i = 0; i < G; i++){
        for(int j = 0; j < nc; j++){
            data_dvd(i,j) = data(i,j) / sf(j);
        }
    }
    whole_mean = data_dvd.rowwise().mean();
    MatrixXd cond_var(G,K);
    int tmp_mean = 0;
    int tmp_var = 0;
    for(int i = 0; i < G; i++)
        for(int j = 0; j < K; j++){
            vector<int> s = which(conditions,j + 1);
            tmp_mean = 0;
            tmp_var = 0;
            for(int t = 0; t < s.size(); t++)
                tmp_mean += data_dvd(i,t);
            tmp_mean /= s.size();
            for(int t = 0; t < s.size(); t++)
                tmp_var += (data(i,t) - tmp_mean) * (data(i,t) - tmp_mean) / sf(t);
            cond_var(i,j) = tmp_var / s.size();
        }
    
    
    
    //    for(int i = 0; i < M; ++i){
    //        size_t sub_nc = reorg_data[i].cols();
    ////        if(sub_nc == 1)
    ////            sub_nc = 2;
    //        MatrixXd center = (reorg_data[i].colwise()-reorg_data[i].rowwise().mean());
    //        cond_var.col(i) = center.rowwise().squaredNorm() / sub_nc;
    //    }
    
    VectorXd var(G);
    var = cond_var.rowwise().mean();
    MatrixXd res(G,nc);
    VectorXd q(G);
    VectorXd I(G);
    I.fill(1.0);
    
    for(int i = 0; i < G; ++i){
        if(abs(var(i) - 0) < 0.0001)
            var(i) = 1;
        if(var(i) <= whole_mean(i))
            q(i) = 0.99;
        else
            q(i) = whole_mean(i) / var(i);
    }
    
    res = ((whole_mean.cwiseProduct(q)).array() / (I - q).array()).matrix() * sf.transpose();
    r_q res_r_q;
    res_r_q.rrr = res;
    res_r_q.qqq = q;
    return res_r_q;
}

MatrixXd DATA::cal_hp(void){
    MatrixXd hyper(G,2);
    MatrixXd cond_q(G,K);
    MatrixXd cond_mean(G,K);
    MatrixXd cond_var(G,K);
    for(int j=0;j<K;++j){
        cond_mean.col(j)=r_d[j].rowwise().mean();
        size_t sub_nc=r_d[j].cols();
        if(sub_nc==1)
            sub_nc=2;
        cond_var.col(j)=(r_d[j].colwise()-cond_mean.col(j)).rowwise().squaredNorm()/(sub_nc-1);
    }
    for(int i=0;i<G;++i){
        for(int j=0;j<K;++j){
            if(cond_mean(i,j)>=cond_var(i,j))
                cond_q(i,j)=0.9999;
            else
                cond_q(i,j)=cond_mean(i,j)/cond_var(i,j);
        }
        double tm;
        double tmean;
        double tvar;
        tmean=cond_q.row(i).mean();
        tvar=(cond_q.row(i).array()-tmean).matrix().squaredNorm()/(cond_q.row(i).cols()-1);
        if(tvar<0.0001){
            tm=100000;
        }
        else{
        tm=tmean*(1-tmean)/tvar;
        }
        hyper(i,0)=tmean*(tm-1);
        hyper(i,1)=(1-tmean)*(tm-1);
        
    }
    return hyper;
}





DATA::DATA(MatrixXd& data, const VectorXi& conditions, const VectorXd& sf, const MatrixXi& pp){
    MatrixXd r;
    this->data = data;
    this->G = data.rows();
    this->r_d = reorg(data,conditions);
    this->K = r_d.size();
    r_q res_rq;
    res_rq = cal_r(data,conditions,sf);
    this->r = res_rq.rrr;
    this->q = res_rq.qqq;
    this->r_r = reorgr(this->r,conditions);
    this->pat = pp;
    
    MatrixXd td(G,K);
    MatrixXd tr(G,K);
    for(int j=0;j<K;j++){
        td.col(j)=this->r_d[j].rowwise().sum();
        tr.col(j)=this->r_r[j].rowwise().sum();
    }
    this->d_s=td;
    this->r_s=tr;
}

DATA::~DATA(){
    q.resize(0);
    r.resize(0,0);
    vector<MatrixXd>().swap(r_d);
    vector<MatrixXd>().swap(r_r);
    pat.resize(0,0);
}


