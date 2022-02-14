Eigen::VectorXd fraction_vec(Eigen::VectorXd x);
Eigen::VectorXd fraction_log(Eigen::VectorXd x);

template <class num>
bool gt(num i, num j);

template <class num>
bool lt(num i, num j);

Eigen::VectorXd sparse_colsum(const Eigen::SparseMatrix<double> &x);
Eigen::SparseMatrix<double> col_prop(const Eigen::SparseMatrix<double> &x);
Eigen::VectorXi vec_cmplmnt(Eigen::VectorXi x, int n, int &xcn);
