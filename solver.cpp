// solver.cpp
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace Eigen;

struct SolveResult {
    VectorXd x; // lysi
    int m2; // # agnwstwn
    int n_nodes;
};


VectorXd compute_preconditioner_inverse(const MatrixXd& A) {
    int n = A.rows();
    VectorXd M_inv(n); // Μ ^ -1
    for (int i = 0; i < n; ++i) {
        double val = A(i, i);
        if (std::abs(val) < 1e-16) {
            M_inv(i) = 1.0; //aii=0, tote 1
        } 
        else {
            M_inv(i) = 1.0 / val; // 1/aii
        }
    }
    return M_inv;
}


VectorXd solve_cg(const MatrixXd& A, const VectorXd& b, double itol, int max_iter) {
    int n = A.rows(); // diastasi
    VectorXd x = VectorXd::Zero(n); // x0 = 0
    VectorXd r = b - A * x;         // arxiko ypoloipo 
    
    double b_norm = b.norm(); // gia elexgo sygklishs 
    if (b_norm == 0) b_norm = 1.0; // gia /0

    // preconditioner 
    VectorXd M_inv = compute_preconditioner_inverse(A); // M = diag A 
    VectorXd z = r.cwiseProduct(M_inv); // z = M-1 * r
    
    VectorXd p = z; // kateythinsi arxika
    double rho = r.dot(z); // rho =  r^T * z
    
    int iter = 0;
    while ( (r.norm() / b_norm > itol) && (iter < max_iter) ) { // oso sfalma > itol
        iter++;
        
        VectorXd q = A * p;
        double alpha = rho / p.dot(q); // rho / (p^T * A * p)
        
        x = x + alpha * p; // x new
        r = r - alpha * q; // new ypoloipo 
        
        
        // gia epomeno bhma
        VectorXd z_new = r.cwiseProduct(M_inv); // znew = M-1 * r
        double rho_new = r.dot(z_new); // rho new =  r^T * z
        
        double beta = rho_new / rho; // beta = rho new / rho old
        p = z_new + beta * p; // new direction
        
        rho = rho_new;
    }
    return x;
}


VectorXd solve_bicg(const MatrixXd& A, const VectorXd& b, double itol, int max_iter) {
    int n = A.rows();
    VectorXd x = VectorXd::Zero(n); // Initial guess
    
    // arxika
    VectorXd r = b - A * x;
    VectorXd r_tilde = r; // r_tilde_0 = r_0
    
    double b_norm = b.norm();
    if (b_norm == 0) b_norm = 1.0;

    VectorXd M_inv = compute_preconditioner_inverse(A);
    
    // Initial preconditioner solve
    VectorXd z = r.cwiseProduct(M_inv); // z = M^ -1 * r
    VectorXd z_tilde = r_tilde.cwiseProduct(M_inv); // M^T * z_tilde = r tilde
    
    double rho = r_tilde.dot(z); // rho = r_tilde^T * z
    if (std::abs(rho) < 1e-14) return x; // Failure 

    VectorXd p = z;
    VectorXd p_tilde = z_tilde; 
    
    int iter = 0;
    while ( (r.norm() / b_norm > itol) && (iter < max_iter) ) {
        iter++;
        
        VectorXd q = A * p;
        VectorXd q_tilde = A.transpose() * p_tilde;
        
        double omega = p_tilde.dot(q);
        if (std::abs(omega) < 1e-14) break; // Failure
        
        double alpha = rho / omega;
        
        x = x + alpha * p;
        r = r - alpha * q;
        r_tilde = r_tilde - alpha * q_tilde;
        
        // next step prep
        VectorXd z_new = r.cwiseProduct(M_inv);
        VectorXd z_tilde_new = r_tilde.cwiseProduct(M_inv);
        
        double rho_new = r_tilde.dot(z_new);
        if (std::abs(rho_new) < 1e-14) break;
        
        double beta = rho_new / rho;
        
        p = z_new + beta * p;
        p_tilde = z_tilde_new + beta * p_tilde;
        
        rho = rho_new;
    }
    return x;
}

// -------------------- Custom LU --------------------------------
VectorXd custom_solve_lu(const MatrixXd& A_in, const VectorXd& b_in) {
    int n = A_in.rows();
    MatrixXd A = A_in;
    VectorXd b = b_in;

    for (int k=0; k<n; ++k) {
        // pivot
        int pivot = k;
        double maxv = std::abs(A(k,k));
        for (int i=k+1; i<n; ++i) {
            if (std::abs(A(i,k)) > maxv) {
                maxv = std::abs(A(i,k));
                pivot = i;
            }
        }
        if (pivot != k) {
            A.row(k).swap(A.row(pivot));
            std::swap(b(k), b(pivot));
        }

        for (int i=k+1; i<n; ++i) {
            A(i,k) /= A(k,k);
            for (int j=k+1; j<n; ++j)
                A(i,j) -= A(i,k) * A(k,j);
        }
    }

    VectorXd y(n), x(n);

    for (int i=0; i<n; ++i) {
        double sum = b(i);
        for (int j=0; j<i; ++j) sum -= A(i,j)*y(j);
        y(i) = sum;
    }

    for (int i=n-1; i>=0; --i) {
        double sum = y(i);
        for (int j=i+1; j<n; ++j) sum -= A(i,j)*x(j);
        x(i) = sum / A(i,i);
    }
    return x;
}

// ------------------- Custom Cholesky ----------------------------
VectorXd custom_solve_cholesky(const MatrixXd& A_in, const VectorXd& b_in) {
    int n = A_in.rows();
    MatrixXd L = MatrixXd::Zero(n,n);

    for (int i=0; i<n; ++i) {
        for (int j=0; j<=i; ++j) {
            double sum = A_in(i,j);
            for (int k=0; k<j; ++k)
                sum -= L(i,k)*L(j,k);

            if (i == j)
                L(i,j) = std::sqrt(std::max(sum, 1e-15));
            else
                L(i,j) = sum / L(j,j);
        }
    }

    VectorXd y(n), x(n);

    for (int i=0; i<n; ++i) {
        double sum = b_in(i);
        for (int j=0; j<i; ++j) sum -= L(i,j)*y(j);
        y(i) = sum / L(i,i);
    }

    for (int i=n-1; i>=0; --i) {
        double sum = y(i);
        for (int j=i+1; j<n; ++j) sum -= L(j,i)*x(j);
        x(i) = sum / L(i,i);
    }
    return x;
}

// ----------------------- SOLVER ----------------------------------------
SolveResult solve_system(element* head, int num_nodes, const RunOptions& opts, MatrixXd A, VectorXd rhs) {
    
    int dim = A.rows();
    int m2 = 0;
    // ypologismos agnwstwn 
    for (element* e=head; e; e=e->next) {
        if (e->type == element::V || e->type == element::L) {
            m2++;
        }
    }

    VectorXd x;

    // ITERATIVE SOLVERS
    if (opts.use_iter) {
        int max_iter = dim + 100; 
        
        if (opts.use_spd) {
            // CG
            // cout << "MPHKA\n";
            x = solve_cg(A, rhs, opts.itol, max_iter);
        } else {
            // Bi-CG
            // cout << "MPHKA\n";
            x = solve_bicg(A, rhs, opts.itol, max_iter);
        }
    }
    // CUSTOM
    else if (opts.use_custom) { 
        if (opts.use_spd) {
            x = custom_solve_cholesky(A, rhs);
        }
        else {
            x = custom_solve_lu(A, rhs);
        }
    }
    // EIGEN
    else {
        if (opts.use_spd) { 
            Eigen::LLT<MatrixXd> LL(A);
            if (LL.info() == Eigen::Success)
                x = LL.solve(rhs);
            else    
                x = PartialPivLU<MatrixXd>(A).solve(rhs);
        }
        else { 
            x = PartialPivLU<MatrixXd>(A).solve(rhs);
        }
    }

    SolveResult sr;
    sr.x = x;
    sr.m2 = m2;
    sr.n_nodes = num_nodes;

    return sr;
}

// --------------------------DC ANALYSIS--------------------------------
void write_dc_op(const std::string& filename, element* head, const SolveResult& sol) {
    std::ofstream ofs(filename);

    for (int i=1; i<=sol.n_nodes; ++i)
        ofs << "NODE_" << i << " = " << sol.x(i-1) << "\n";

    int idx = sol.n_nodes;
    for (element* e=head; e; e=e->next) {
        if (e->type == element::V || e->type == element::L) {
            ofs << e->name << " = " << sol.x(idx) << "\n";
            idx++;
        }
    }
}
