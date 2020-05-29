
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <progress.hpp>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]

#define EIGEN_DONT_PARALLELIZE

template <class num>
bool gt(num i, num j){ return (i > j); }

// Project y onto simplex, return x
Eigen::VectorXd proj_smplx(Eigen::VectorXd y){
    int D = y.size();
    Eigen::VectorXd u = y;
    
    std::sort(u.data(), u.data() + u.size(), gt<double>); // descending

    double csum = 0;
    Eigen::VectorXd p(D);
    int max_j = 0;
    for (double j = 0; j < D; j++){
        csum += u(j);
        p(j) = (1. - csum) / (j + 1.);
        double ui = u(j) + p(j);
        if (ui > 0) max_j = j;
    }
    double l = p(max_j);
    Eigen::VectorXd x = y;
    for (int i = 0; i < D; i++){
        x(i) = y(i) + l;
        x(i) = (x(i) < 0) ? 0 : x(i);
    }
    return x;
}



// Compute log density of Dirichlet-Multinomial PMF with sparse matrix
// x is a gene by droplet sparse matrix
// [[Rcpp::export]]
Eigen::MatrixXd fit_prop(Eigen::SparseMatrix<double> x, 
        Eigen::MatrixXd alpha, 
        Eigen::MatrixXd w, 
        double lrate = 0.1, 
        double eps = 0.001, 
        bool accelerate = false, 
        int threads = 1, 
        bool display_progress = true){
    int n_g = x.rows();
    int n_c = x.cols();
    if ( alpha.rows() != n_g ){
        stop("The number of rows in alpha must be the same as the number of rows in x");
    }

#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
    if (display_progress){
        REprintf("Number of threads=%i\n", threads);
    }
#endif

    int K = alpha.cols();
    if ( w.rows() != n_c ){
        stop("The number of rows in w must be the same as the number of columns in x");
    }
    if ( w.cols() != K ){
        stop("The number of columns in w must be the same as the number of columns in alpha");
    }

    Eigen::MatrixXd w_tm2(n_c, K);
    w_tm2.setZero();
    Eigen::MatrixXd w_new(n_c, K);
    Eigen::MatrixXd w_grad(n_c, K);

    Progress p(n_c, display_progress);


#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 10)
#endif
    for (int i = 0; i < n_c; ++i){
        if ( ! Progress::check_abort() )
            p.increment();
        double delta = eps + 1;
        double iter = 1;
        while (delta >= eps){
            Eigen::VectorXd v(K);
            if (accelerate){
                for (int k = 0; k < K; k++){
                    v(k) = w(i, k) + (w(i,k) - w_tm2(i,k))*(iter - 2)/(iter + 1);
                }
            }
            for (int k = 0; k < K; k++){
                w_grad(i,k) = 0;
                for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it){
                    double wa_sum = 0;
                    for (int k2 = 0; k2 < K; k2++){
                        if (accelerate)
                            wa_sum += v(k2) * alpha(it.index(), k2);
                        else
                            wa_sum += w(i,k2) * alpha(it.index(), k2);
                    }
                    w_grad(i,k) += (it.value() * alpha(it.index(), k)) / wa_sum;
                }
                // update w
                if (accelerate)
                    w_new(i, k) = v(k) + ((lrate / iter) * w_grad(i,k));
                else
                    w_new(i, k) = w(i,k) + ((lrate / iter) * w_grad(i,k));
            }
            Eigen::VectorXd oldpv = w_new.row(i);
            w_new.row(i) = proj_smplx(oldpv);

            delta = 0;
            for (int k = 0; k < K; k++){
                delta += abs(w_new(i,k) - w(i,k));
                w_tm2(i,k) = w(i,k);
                w(i,k) = w_new(i,k);
            }
            iter++;
        }
    }
    return w_new;
}

