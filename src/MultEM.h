#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


Eigen::MatrixXd LlkMultSparsePar(const Eigen::SparseMatrix<double> &x, 
        const Eigen::VectorXd &sizes,
        const Eigen::MatrixXd &alpha, 
        Eigen::VectorXi ix, 
        int ix_len, 
        int threads = 1, 
        bool display_progress = true,
        bool debug = false);

Eigen::MatrixXd mult_alpha_ml(const Eigen::SparseMatrix<double> &x, 
        const Eigen::MatrixXd &Z, 
        double prior = 0, 
        int threads = 1, 
        double psc = 1e-10);

double mult_z_ml(const Eigen::MatrixXd &llk, 
        const Eigen::VectorXd &pi, 
        const Eigen::VectorXi &ix, 
        const int &ix_len, 
        Eigen::MatrixXd &Z, 
        int threads = 1);

Eigen::VectorXd mult_pi_ml(const Eigen::MatrixXd &Z, 
        double prior = 0, 
        int threads = 1);

List mult_em(Eigen::SparseMatrix<double> X,
             Eigen::MatrixXd Alpha, 
             Eigen::VectorXd Pi,
             Eigen::MatrixXd Z,
             Eigen::VectorXi fixed,
             double Alpha_prior = 0,
             double Pi_prior = 0, 
             int threads = 1,
             int max_iter = 100, 
             double eps = 1e-2, 
             double psc = 1e-10, 
             bool display_progress = true);
