//
//  Header.cpp
//  test for map<matrixXd>
//
//  Created by MaXiuyu on 6/9/17.
//  Copyright © 2017 MaXiuyu. All rights reserved.
//

#include "Header.hpp"



V unique(V x){
    V ref;
    ref.push_back(x[0]);
    for(V::iterator iter=x.begin(); iter!=x.end();++iter){
        V::iterator curr;
        curr=find(ref.begin(),ref.end(),*iter);
        if(curr==ref.end())
            ref.push_back(*iter);
    }
    return ref;
}



//all possible partitions of a set
vector<V> partition(const int& n){
    vector<vector<int> > start(1);
    start[0].push_back(1);
    if(n==1){
        return start;
    }
    
    for(int i=1;i<n;++i){
        vector<V> new_p;
        size_t L=start.size();
        for(int j=0;j<L;++j){
            int M=*max_element(start[j].begin(),start[j].end());
            for(int k=0;k<M;++k){
                start[j].push_back(k+1);
                new_p.push_back(start[j]);
                start[j].pop_back();
            }
            start[j].push_back(M+1);
            new_p.push_back(start[j]);
            start[j].pop_back();
        }
        start=new_p;
    }
    return start;
}

//which function in r

vector<int> which (const vector<int>& x, const int& index){
    vector<int> v;
    for(vector<int>::size_type i=0,length=x.size();i<length; ++i){
        if(x[i]==index){
            v.push_back(int(i));
        }
    }
    return v;
}

vector<int> which (const VectorXi x, const int& index){
    vector<int> v;
    for(size_t i=0,length=x.size();i<length; ++i){
        if(x(i)==index){
            v.push_back(int(i));
        }
    }
    return v;
}

Sum rc(MatrixXd x,MatrixXd r, vector<V> pos){
    Sum res;
    size_t npt = pos.size();
    size_t G = x.rows();
    MatrixXd sum(G,npt);
    sum.fill(0);
    MatrixXd R(G,npt);
    R.fill(0);
    int num = 0;
    
    for(vector<V>::iterator it = pos.begin(); it != pos.end(); ++it){
        for(V::iterator it1 = (*it).begin(); it1 != (*it).end(); ++it1){
            sum.col(num) += x.col(*it1);
            R.col(num) += r.col(*it1);
        }
        ++num;
    }
    res.count_sum=sum.cast<double>();
    res.r_sum=R;
    return res;
}


