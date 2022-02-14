

Eigen::VectorXd fast_varCPP(const Eigen::SparseMatrix<double> &x, 
        const Eigen::VectorXd &mu, 
        int threads = 1, 
        bool display_progress = false);

Eigen::VectorXd fast_wvarCPP(const Eigen::SparseMatrix<double> &x, 
        const Eigen::VectorXd &mu, 
        const Eigen::VectorXd &weights,
        int threads = 1, 
        bool display_progress = false);

Eigen::VectorXd fast_wmeanCPP(const Eigen::SparseMatrix<double> &x, 
        const Eigen::VectorXd &weights,
        int threads = 1, 
        bool display_progress = false);

