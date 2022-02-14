#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <math.h>
#include <limits>
using namespace Rcpp;
using namespace std;

typedef Eigen::SparseMatrix<double>::InnerIterator cit;

template <class num>
bool gt(num i, num j){ return (i > j); }

template <class num>
bool lt(num i, num j){ return (i < j); }

Eigen::VectorXd fraction_vec(Eigen::VectorXd x){
    double xsum = 0;
    for (int i = 0; i < x.size(); i++){
        xsum += x(i);
    }
    for (int i = 0; i < x.size(); i++){
        x(i) /= xsum;
    }
    return x;
}


// fraction of logs
// @param x numeric vector
Eigen::VectorXd fraction_log(Eigen::VectorXd x){
    int inf_ix = -1;
    for (int i = 0; i < x.size(); i++){
        if ( x(i) == std::numeric_limits<double>::infinity()){
            inf_ix = i;
            break;
        }
    }
    if (inf_ix >= 0){
        for (int i = 0; i < x.size(); i++){
            x(i) = 0;
        }
        x(inf_ix) = 1;
        return x;
    }

    double xc = x(0);
    for (int i = 1; i < x.size(); i++){
        if ( x(i) > xc )
            xc = x(i);
    }
    double xsum = 0;
    for (int i = 0; i < x.size(); i++){
        x(i) -= xc;
        x(i) = std::exp(x(i));
        xsum += x(i);
    }
    for (int i = 0; i < x.size(); i++){
        x(i) /= xsum;
    }
    return x;
}

Eigen::VectorXd sparse_colsum(const Eigen::SparseMatrix<double> &x){
    int N = x.cols();
    Eigen::VectorXd sizes(N);
    sizes.setZero();
    for (int i = 0; i < N; i++){
        for (cit it(x, i); it; ++it){
            sizes(i) += it.value();
        }
    }
    return sizes;
}

Eigen::SparseMatrix<double> col_prop(const Eigen::SparseMatrix<double> &x){
    int M = x.rows();
    int N = x.cols();
    Eigen::SparseMatrix<double> y(x);
    Eigen::VectorXd colsums = sparse_colsum(x);
    for (int i = 0; i < N; i++){
        for (cit it(y, i); it; ++it){
            it.valueRef() /= colsums(i);
        }
    }
    return y;
}


Eigen::VectorXi vec_cmplmnt(Eigen::VectorXi x, int n, int &xcn){
    // free and fixed observations
    int xn = x.size();
    if (xn > 0){
        std::sort(x.data(), x.data() + x.size(), lt<int>); // ascending
    } else {
        x = Eigen::VectorXi(1);
        x(0) = -1;
    }

    xcn = 0;
    Eigen::VectorXi xc(n);
    int j = 0;
    for (int i = 0; i < n; i++){
        if (i != x(j)){
            xc(xcn) = i;
            xcn++;
        }
        if (i == x(j)){
            j++;
        }
    }
    xc.conservativeResize(xcn);
    return xc;
}

