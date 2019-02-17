#ifndef kit_h
#define kit_h
#include <Rcpp.h>
#include <RcppEigen.h>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <math.h>
using namespace std;
typedef vector<int>::iterator it;
typedef vector<int> V;



//template <typename T>
//vector<size_t> sort_indexes(const vector<T> &v) {
//    
//    // initialize original index locations
//    vector<size_t> idx(v.size());
//    iota(idx.begin(), idx.end(), 0);
//    
//
//    
//    // sort indexes based on comparing values in v
//    sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
//    
//    return idx;
//}


class component{
public:
    int tailsize;
    int tailsum;
    double logdensity;
    double logprior;
    V partition;
    component(){
        tailsize=0;
        logdensity=0;
        logprior=0;
        tailsum=0;
    }
    
    
    bool operator==(const component & a){
        return (logdensity+logprior)==(a.logdensity+a.logprior);
    }
    bool operator<(const component& a){
        return (logdensity+logprior)<(a.logdensity+a.logprior);
    }
    bool operator>(const component& a){
        return !(*this==a)&&(*this<a);
    }
    void operator=(const component& a){
        logdensity=a.logdensity;
        logprior=a.logprior;
        tailsum=a.tailsum;
        partition=a.partition;
    }
    
};


bool myfn(const component & a, const component & b) { return (a.logdensity+a.logprior)<(b.logdensity+b.logprior); }

//if c++ version not support iota function

//template <class ForwardIterator, class T>
// void iota (ForwardIterator first, ForwardIterator last, T val)
//{
//    while(first!=last){
//        *first=val;
//        ++first;
//        ++val;
//    }
//}



#endif