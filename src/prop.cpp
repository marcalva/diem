
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

typedef Eigen::SparseMatrix<double>::InnerIterator cit;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator rit;

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
        Eigen::MatrixXd theta, 
        Eigen::MatrixXd w, 
        double lrate = 0.1, 
        double eps = 0.001, 
        double psc = 0,
        int threads = 1, 
        int max_iter = 1e4, 
        bool display_progress = true){
    int m = x.rows();
    int n = x.cols();
    if ( theta.rows() != m ){
        stop("The number of rows in theta must be the same as the number of rows in x");
    }

#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
    if (display_progress){
        REprintf("Number of threads=%i\n", threads);
    }
#endif

    int K = theta.cols();
    if ( w.rows() != n ){
        stop("The number of rows in w must be the same as the number of columns in x");
    }
    if ( w.cols() != K ){
        stop("The number of columns in w must be the same as the number of columns in theta");
    }

    Eigen::SparseMatrix<double> v(x);

    Progress p(n, display_progress);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 10)
#endif
    for (int i = 0; i < n; ++i){
        if ( ! Progress::check_abort() )
            p.increment();
        double delta = eps;
        double iter = 1;
        Eigen::VectorXd w_grad(K);
        Eigen::VectorXd w_iter(K);
        while (delta >= eps & iter <= max_iter){
            for (cit it(v,i); it; ++it){
                it.valueRef() = 0;
                for (int k = 0; k < K; k++){
                    it.valueRef() += w(i, k) * theta(it.index(), k);
                }
                it.valueRef() += psc;
            }
            for (int k = 0; k < K; k++){
                w_grad(k) = 0;
                for (cit it_x(x,i), it_v(v, i); it_x; ++it_x, ++it_v){
                    w_grad(k) += (it_x.value() * theta(it_x.index(), k)) / it_v.value();
                }
                // update w
                w_iter(k) = w(i,k) + ((lrate / iter) * w_grad(k));
            }
            // Rcout << "sum for iter " << iter << ": " << w_iter.sum() << "\n";
            w_iter = proj_smplx(w_iter);

            delta = 0;
            for (int k = 0; k < K; k++){
                delta += abs(w_iter(k) - w(i,k));
                w(i,k) = w_iter(k);
            }
            iter++;
        }
    }
    return w;
}

double delta_m(Eigen::MatrixXd x1, Eigen::MatrixXd x2){
    double delta = 0;
    for (int i = 0; i < x1.rows(); i++){
        for (int j = 0; j < x1.cols(); j++){
            delta += abs(x1(i,j) - x2(i,j));
        }
    }
    return delta;
}

Eigen::VectorXd sum_to_one(Eigen::VectorXd x, double psc = 1e-8){
    double xs = 0;
    for (int i = 0; i < x.size(); i++){
        x(i) += psc;
        xs += x(i);
    }
    for (int i = 0; i < x.size(); i++){
        x(i) /= xs;
    }
    return x;
}

/*****************************************
 *
 *
 * ***************************************/

/* Estimation of mixture and cluster parameters 
 *
 * @param x An M by N count matrix.
 * @param theta Initialized M x K parameter values.
 * @param w Initialized N x K weights matrix.
 *
 * */
//
// Compute log density of Dirichlet-Multinomial PMF with sparse matrix
// x is a gene by droplet sparse matrix
// [[Rcpp::export]]
List fitm(Eigen::SparseMatrix<double> x, 
        Eigen::MatrixXd theta, 
        Eigen::MatrixXd w, 
        Eigen::VectorXi f,
        double eps1 = 0.001,
        double eps2 = 0.001, 
        int max_iter1 = 200, 
        int max_iter2 = 100,
        double alpha = 0, 
        double beta = 0, 
        double lrate_w = 0.1, 
        double lrate_theta = 1e-8, 
        double psc = 1e-8, 
        int threads = 1, 
        bool display_progress = true){

    x.makeCompressed();
    Eigen::SparseMatrix<double,Eigen::RowMajor> xr(x);
    
    int m = x.rows();
    int n = x.cols();
    int n_free = f.size();

    Rcout << "f(5): " << f(4) << "\n";
    
    if ( theta.rows() != m ){
        stop("The number of rows in theta must be the same as the number of rows in x");
    }

#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
    if (display_progress){
        REprintf("Number of threads=%i\n", threads);
    }
#endif

    int K = theta.cols();
    if ( w.rows() != n ){
        stop("The number of rows in w must be the same as the number of columns in x");
    }
    if ( w.cols() != K ){
        stop("The number of columns in w must be the same as the number of columns in theta");
    }

    Eigen::SparseMatrix<double> v = x;
    Eigen::MatrixXd w_old(w);
    Eigen::MatrixXd theta_new(m, K);
    Eigen::MatrixXd theta_grad(m, K);
    
    double delta1 = eps1 + 1;
    int iter1 = 1;
    while (delta1 >= eps1 && iter1 <= max_iter1){
        if (display_progress)
            Rcout << "Iteration: " << iter1 << "\n";
        // Get v
        if (display_progress)
            Rcout << "Getting V \n";
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 10)
#endif
        for (int i = 0; i < n; i++){
            for (cit it(v,i); it; ++it){
                double wt_sum = 0;
                for (int k = 0; k < K; k++){
                    wt_sum += w(i,k) * theta(it.index(), k);
                }
                // Rcout << it.value() << " ";
                it.valueRef() = wt_sum + psc;
                // Rcout << it.value() << "\n";
            }
        }

        // Get w_i for free columns
        if (display_progress)
            Rcout << "Getting w \n";
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 10)
#endif
        for (int ix = 0; ix < n_free; ix++){
            int i = f(ix) - 1;
            double delta2 = eps2 + 1;
            double iter2 = 1;
            // Gradient ascent
            Eigen::VectorXd w_iter(K);
            Eigen::VectorXd w_grad(K);
            while (delta2 >= eps2 && iter2 <= max_iter2){
                // w_ik
                for (int k = 0; k < K; k++){
                    w_grad(k) = 0;
                    // For each gene (rows)
                    cit it_v(v,i);
                    for (cit it_x(x,i); it_x; ++it_x){
                        w_grad(k) += (it_x.value() * theta(it_x.index(), k)) / it_v.value();
                        ++it_v;
                    }
                    w_grad(K) += alpha / w(i,k);
                    w_iter(k) = w(i,k) + ((lrate_w / iter2) * w_grad(k));
                }
                w_iter = proj_smplx(w_iter);
                w_iter = sum_to_one(w_iter, psc);

                delta2 = 0;
                for (int k = 0; k < K; k++){
                    delta2 += abs(w_iter(k) - w(i,k));
                    w(i,k) = w_iter(k);
                }
                // Rcout << "iter2 " << iter2 << " delta2 " << delta2 << "\n";
                iter2++;
            }
        }


        // Get v
        if (display_progress)
            Rcout << "Getting V \n";
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 10)
#endif
        for (int i = 0; i < n; i++){
            for (cit it(v,i); it; ++it){
                double wt_sum = 0;
                for (int k = 0; k < K; k++){
                    wt_sum += w(i,k) * theta(it.index(), k);
                }
                it.valueRef() = wt_sum + psc;
            }
        }
        
        // RowMajor V matrix
        Eigen::SparseMatrix<double,Eigen::RowMajor> vr(v);

        Eigen::MatrixXd theta_old(theta);
        // Get theta_k
        if (display_progress)
            Rcout << "Getting theta \n";
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 10)
#endif
        for (int k = 0; k < K; k++){
            double delta2 = eps2 + 1;
            double iter2 = 1;
            // Gradient ascent
            Eigen::VectorXd theta_iter(m);
            Eigen::VectorXd theta_grad(m);
            // Rcout << "k: " << k << "\n";
            while (delta2 >= eps2 && iter2 <= max_iter2){
                // theta_gk
                for (int j = 0; j < m; j++){
                    theta_grad(j) = 0;
                    // For each sample (columns)
                    rit it_v(vr,j);
                    for (rit it_x1(xr,j); it_x1; ++it_x1){
                        theta_grad(j) += (it_x1.value() * w(it_x1.index(),k)) / it_v.value();
                        ++it_v;
                    }
                    theta_grad(j) += beta / theta(j,k) + psc;
                    // Rcout << "grad " << j << " " << theta_grad(j) << "\n";
                    theta_iter(j) = theta(j,k) + ((lrate_theta / iter2) * theta_grad(j));
                }
                //Rcout << "Sum: " <<theta_iter.sum() << "\n";
                theta_iter = proj_smplx(theta_iter);
                theta_iter = sum_to_one(theta_iter, psc);

                delta2 = 0;
                for (int j = 0; j < m; j++){
                    delta2 += abs(theta_iter(j) - theta(j,k));
                    theta(j,k) = theta_iter(j);
                }
                //Rcout << "iter2 " << iter2 << " delta2 " << delta2 << "\n";
                iter2++;
            }
        }

        delta1 = 0;
        delta1 += delta_m(w, w_old);
        Rcout << "delta1 1" << delta1 << "\n";
        w_old = w;
        delta1 += delta_m(theta, theta_old);
        delta1 = delta1 / (K * K);
        theta_old = theta;
        if (display_progress)
            Rcout << "delta: " << delta1 << "\n";

        // Calculate delta change in parameter values
        iter1++;
    }

    List L = List::create(w, theta);
    return L;

}

