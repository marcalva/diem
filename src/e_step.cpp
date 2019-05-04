
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double mvn_logllk_diagCPP(NumericVector x, 
		NumericVector mu, 
		NumericVector sgma){
	int np = sgma.size();
	if ( (x.size() != mu.size()) || (mu.size() != sgma.size()) ){
		Rcpp::stop("x, mu, and sgma must be the same length");
	}
	double logdet = 0;
	NumericVector sgma_inv(np);
	for (int i = 0; i < np; i++){
    	logdet += log(sgma(i));
    	sgma_inv(i) = 1.0/sgma(i);
	}
	
	NumericVector x_c = x - mu;
	
	double xe = 0;
	for (int i = 0; i < np; i++){
		xe += x_c(i) * x_c(i) * sgma_inv(i);
	}

    double con = np * log(2 * M_PI);

    double ret = -0.5 * (logdet + xe + con);
    return ret;
}

// [[Rcpp::export]]
NumericVector fraction_logCPP(NumericVector x){
	NumericVector x_c = x - max(x);
	x_c = exp(x_c);
	NumericVector frac = x_c / sum(x_c);
	return frac;
}

// [[Rcpp::export]]
double sum_logCPP(NumericVector x){
	double max_x = max(x);
	NumericVector x_c = x - max_x;
	double x_sum = log( sum( exp(x_c) ) ) + max_x;
    return x_sum;
}


// [[Rcpp::export]]
NumericMatrix e_stepCPP(NumericMatrix x, 
		int k, 
		List mu, 
		List sgma, 
		NumericVector tau,
		bool semisup = false,
		NumericVector labels = R_NilValue){
	int n = x.nrow();
	NumericMatrix Z(n, k);
	NumericVector f_i(k);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < k; j++){
			NumericVector x_i = x(i, _);
			NumericVector mu_i = mu[j];
			NumericVector sgma_i = sgma[j];
			f_i(j) = log(tau(j)) + mvn_logllk_diagCPP(x_i, mu_i, sgma_i);
		}
		// unlabeled data is 0. Shift to 0-based index
		if (semisup && (labels(i) != 0)){
			Z(i, labels(i)-1) = labels(i);
		}
		else{
			for (int j = 0; j < k; j++){
				NumericVector x_i = x(i, _);
				NumericVector mu_i = mu[j];
				NumericVector sgma_i = sgma[j];
				f_i(j) = log(tau(j)) + mvn_logllk_diagCPP(x_i, mu_i, sgma_i);
			}
			Z(i,_) = fraction_logCPP(f_i);
		}
	}
	return Z;
}

// [[Rcpp::export]]
double get_llkCPP(NumericMatrix x, 
		List mu, 
		List sgma, 
		NumericVector tau, 
		bool semisup = false, 
		NumericVector labels = R_NilValue){
	// double eps = 0.01;
	// if (sum(tau) < (1 - eps) | sum(tau) > (1+eps)) stop(paste0("Parameter tau must sum to 1, not ", as.character(tau)))
	bool same_len = (mu.size() == sgma.size()) && (sgma.size() == tau.size());
    if (!same_len) {
		Rcpp::stop("mu, sgma, and tau must have same length");
	}
    int n = x.nrow();
    int n_label = labels.size();
    int k = mu.size();
    double llk = 0;
	if (semisup){
		if (n_label != n){
			char buf[100];
			sprintf(buf, "Labels %i must same size as number of rows in x %i\n", n_label, n);
			Rcpp::stop(buf);
		}
	}
	for (int i = 0; i < n; i++){
		if ( semisup && labels(i) != 0 ){
			int ix = labels(i) - 1;
			NumericVector x_i = x(i,_);
			NumericVector mu_i = as<NumericVector>(mu[ix]);
			NumericVector sgma_i = as<NumericVector>(sgma[ix]);
			llk += log(tau(ix)) + mvn_logllk_diagCPP( x_i, mu_i, sgma_i);
		}
		else{
			NumericVector llk_i(k);
			for (int j = 0; j < k; j++){
				NumericVector x_i = x(i,_);
				NumericVector mu_ij = as<NumericVector>(mu[j]);
				NumericVector sgma_ij = as<NumericVector>(sgma[j]);
				llk_i(j) = ( log(tau(j)) + mvn_logllk_diagCPP( x_i, mu_ij, sgma_ij) ); 
			}
			llk += sum_logCPP(llk_i);
		}
	}
	return llk;
}
