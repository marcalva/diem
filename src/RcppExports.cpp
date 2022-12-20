// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// max_loo
Eigen::VectorXd max_loo(const Eigen::SparseMatrix<double>& x, const Eigen::VectorXd& sizes, const Eigen::VectorXd& weights, const Eigen::VectorXd& alpha, double eps, int max_iter, double alpha_prior, double psc, int threads, bool debug);
RcppExport SEXP _diem_max_loo(SEXP xSEXP, SEXP sizesSEXP, SEXP weightsSEXP, SEXP alphaSEXP, SEXP epsSEXP, SEXP max_iterSEXP, SEXP alpha_priorSEXP, SEXP pscSEXP, SEXP threadsSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prior(alpha_priorSEXP);
    Rcpp::traits::input_parameter< double >::type psc(pscSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(max_loo(x, sizes, weights, alpha, eps, max_iter, alpha_prior, psc, threads, debug));
    return rcpp_result_gen;
END_RCPP
}
// LlkDirMultSparsePar
Eigen::MatrixXd LlkDirMultSparsePar(const Eigen::SparseMatrix<double>& x, const Eigen::VectorXd& sizes, const Eigen::MatrixXd& alpha, Eigen::VectorXi ix, int ix_len, int threads, bool display_progress, bool debug);
RcppExport SEXP _diem_LlkDirMultSparsePar(SEXP xSEXP, SEXP sizesSEXP, SEXP alphaSEXP, SEXP ixSEXP, SEXP ix_lenSEXP, SEXP threadsSEXP, SEXP display_progressSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type ix(ixSEXP);
    Rcpp::traits::input_parameter< int >::type ix_len(ix_lenSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(LlkDirMultSparsePar(x, sizes, alpha, ix, ix_len, threads, display_progress, debug));
    return rcpp_result_gen;
END_RCPP
}
// dirmult_alpha_mom
Eigen::MatrixXd dirmult_alpha_mom(const Eigen::SparseMatrix<double>& x, const Eigen::MatrixXd& Z, double alpha_prior, double psc, int threads);
RcppExport SEXP _diem_dirmult_alpha_mom(SEXP xSEXP, SEXP ZSEXP, SEXP alpha_priorSEXP, SEXP pscSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prior(alpha_priorSEXP);
    Rcpp::traits::input_parameter< double >::type psc(pscSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(dirmult_alpha_mom(x, Z, alpha_prior, psc, threads));
    return rcpp_result_gen;
END_RCPP
}
// dirmult_alpha_loo
Eigen::MatrixXd dirmult_alpha_loo(const Eigen::SparseMatrix<double>& x, const Eigen::MatrixXd& alpha0, const Eigen::MatrixXd& Z, double eps, int max_iter, double alpha_prior, double psc, int threads);
RcppExport SEXP _diem_dirmult_alpha_loo(SEXP xSEXP, SEXP alpha0SEXP, SEXP ZSEXP, SEXP epsSEXP, SEXP max_iterSEXP, SEXP alpha_priorSEXP, SEXP pscSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prior(alpha_priorSEXP);
    Rcpp::traits::input_parameter< double >::type psc(pscSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(dirmult_alpha_loo(x, alpha0, Z, eps, max_iter, alpha_prior, psc, threads));
    return rcpp_result_gen;
END_RCPP
}
// dirmult_alpha_ml
Eigen::MatrixXd dirmult_alpha_ml(const Eigen::SparseMatrix<double>& x, const Eigen::MatrixXd& Z, double eps, int max_iter_loo, double alpha_prior, double psc, int threads);
RcppExport SEXP _diem_dirmult_alpha_ml(SEXP xSEXP, SEXP ZSEXP, SEXP epsSEXP, SEXP max_iter_looSEXP, SEXP alpha_priorSEXP, SEXP pscSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter_loo(max_iter_looSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prior(alpha_priorSEXP);
    Rcpp::traits::input_parameter< double >::type psc(pscSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(dirmult_alpha_ml(x, Z, eps, max_iter_loo, alpha_prior, psc, threads));
    return rcpp_result_gen;
END_RCPP
}
// dirmult_em
List dirmult_em(Eigen::SparseMatrix<double> X, Eigen::MatrixXd Alpha, Eigen::VectorXd Pi, Eigen::MatrixXd Z, Eigen::VectorXi fixed, double alpha_prior, double pi_prior, int threads, int max_iter, int max_iter_loo, double eps, double eps_loo, double psc, bool display_progress);
RcppExport SEXP _diem_dirmult_em(SEXP XSEXP, SEXP AlphaSEXP, SEXP PiSEXP, SEXP ZSEXP, SEXP fixedSEXP, SEXP alpha_priorSEXP, SEXP pi_priorSEXP, SEXP threadsSEXP, SEXP max_iterSEXP, SEXP max_iter_looSEXP, SEXP epsSEXP, SEXP eps_looSEXP, SEXP pscSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Alpha(AlphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type fixed(fixedSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prior(alpha_priorSEXP);
    Rcpp::traits::input_parameter< double >::type pi_prior(pi_priorSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter_loo(max_iter_looSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type eps_loo(eps_looSEXP);
    Rcpp::traits::input_parameter< double >::type psc(pscSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(dirmult_em(X, Alpha, Pi, Z, fixed, alpha_prior, pi_prior, threads, max_iter, max_iter_loo, eps, eps_loo, psc, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// LlkMultSparsePar
Eigen::MatrixXd LlkMultSparsePar(const Eigen::SparseMatrix<double>& x, const Eigen::VectorXd& sizes, const Eigen::MatrixXd& alpha, Eigen::VectorXi ix, int ix_len, int threads, bool display_progress, bool debug);
RcppExport SEXP _diem_LlkMultSparsePar(SEXP xSEXP, SEXP sizesSEXP, SEXP alphaSEXP, SEXP ixSEXP, SEXP ix_lenSEXP, SEXP threadsSEXP, SEXP display_progressSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type ix(ixSEXP);
    Rcpp::traits::input_parameter< int >::type ix_len(ix_lenSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(LlkMultSparsePar(x, sizes, alpha, ix, ix_len, threads, display_progress, debug));
    return rcpp_result_gen;
END_RCPP
}
// mult_alpha_ml
Eigen::MatrixXd mult_alpha_ml(const Eigen::SparseMatrix<double>& x, const Eigen::MatrixXd& Z, double prior, int threads, double psc);
RcppExport SEXP _diem_mult_alpha_ml(SEXP xSEXP, SEXP ZSEXP, SEXP priorSEXP, SEXP threadsSEXP, SEXP pscSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< double >::type psc(pscSEXP);
    rcpp_result_gen = Rcpp::wrap(mult_alpha_ml(x, Z, prior, threads, psc));
    return rcpp_result_gen;
END_RCPP
}
// mult_z_ml
double mult_z_ml(const Eigen::MatrixXd& llk, const Eigen::VectorXd& pi, const Eigen::VectorXi& ix, const int& ix_len, Eigen::MatrixXd& Z, int threads);
RcppExport SEXP _diem_mult_z_ml(SEXP llkSEXP, SEXP piSEXP, SEXP ixSEXP, SEXP ix_lenSEXP, SEXP ZSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type llk(llkSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type ix(ixSEXP);
    Rcpp::traits::input_parameter< const int& >::type ix_len(ix_lenSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mult_z_ml(llk, pi, ix, ix_len, Z, threads));
    return rcpp_result_gen;
END_RCPP
}
// mult_pi_ml
Eigen::VectorXd mult_pi_ml(const Eigen::MatrixXd& Z, double prior, int threads);
RcppExport SEXP _diem_mult_pi_ml(SEXP ZSEXP, SEXP priorSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mult_pi_ml(Z, prior, threads));
    return rcpp_result_gen;
END_RCPP
}
// mult_em
List mult_em(Eigen::SparseMatrix<double> X, Eigen::MatrixXd Alpha, Eigen::VectorXd Pi, Eigen::MatrixXd Z, Eigen::VectorXi fixed, double alpha_prior, double pi_prior, int threads, int max_iter, double eps, double psc, bool display_progress);
RcppExport SEXP _diem_mult_em(SEXP XSEXP, SEXP AlphaSEXP, SEXP PiSEXP, SEXP ZSEXP, SEXP fixedSEXP, SEXP alpha_priorSEXP, SEXP pi_priorSEXP, SEXP threadsSEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP pscSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Alpha(AlphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type fixed(fixedSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prior(alpha_priorSEXP);
    Rcpp::traits::input_parameter< double >::type pi_prior(pi_priorSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type psc(pscSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(mult_em(X, Alpha, Pi, Z, fixed, alpha_prior, pi_prior, threads, max_iter, eps, psc, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// fast_varCPP
Eigen::VectorXd fast_varCPP(const Eigen::SparseMatrix<double>& x, const Eigen::VectorXd& mu, int threads, bool display_progress);
RcppExport SEXP _diem_fast_varCPP(SEXP xSEXP, SEXP muSEXP, SEXP threadsSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_varCPP(x, mu, threads, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// fast_wvarCPP
Eigen::VectorXd fast_wvarCPP(const Eigen::SparseMatrix<double>& x, const Eigen::VectorXd& mu, const Eigen::VectorXd& weights, int threads, bool display_progress);
RcppExport SEXP _diem_fast_wvarCPP(SEXP xSEXP, SEXP muSEXP, SEXP weightsSEXP, SEXP threadsSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_wvarCPP(x, mu, weights, threads, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// fast_wmeanCPP
Eigen::VectorXd fast_wmeanCPP(const Eigen::SparseMatrix<double>& x, const Eigen::VectorXd& weights, int threads, bool display_progress);
RcppExport SEXP _diem_fast_wmeanCPP(SEXP xSEXP, SEXP weightsSEXP, SEXP threadsSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_wmeanCPP(x, weights, threads, display_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_diem_max_loo", (DL_FUNC) &_diem_max_loo, 10},
    {"_diem_LlkDirMultSparsePar", (DL_FUNC) &_diem_LlkDirMultSparsePar, 8},
    {"_diem_dirmult_alpha_mom", (DL_FUNC) &_diem_dirmult_alpha_mom, 5},
    {"_diem_dirmult_alpha_loo", (DL_FUNC) &_diem_dirmult_alpha_loo, 8},
    {"_diem_dirmult_alpha_ml", (DL_FUNC) &_diem_dirmult_alpha_ml, 7},
    {"_diem_dirmult_em", (DL_FUNC) &_diem_dirmult_em, 14},
    {"_diem_LlkMultSparsePar", (DL_FUNC) &_diem_LlkMultSparsePar, 8},
    {"_diem_mult_alpha_ml", (DL_FUNC) &_diem_mult_alpha_ml, 5},
    {"_diem_mult_z_ml", (DL_FUNC) &_diem_mult_z_ml, 6},
    {"_diem_mult_pi_ml", (DL_FUNC) &_diem_mult_pi_ml, 3},
    {"_diem_mult_em", (DL_FUNC) &_diem_mult_em, 12},
    {"_diem_fast_varCPP", (DL_FUNC) &_diem_fast_varCPP, 4},
    {"_diem_fast_wvarCPP", (DL_FUNC) &_diem_fast_wvarCPP, 5},
    {"_diem_fast_wmeanCPP", (DL_FUNC) &_diem_fast_wmeanCPP, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_diem(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
