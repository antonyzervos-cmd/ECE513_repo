// solver.cpp
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <suitesparse/cs.h>
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


// --------------------------
// ----------------------- SPARSE SOLVERS ----------------------
// --------------------------

// CSparse xrisimopoiei double* arrays. Opote metatrepo ta VectorXd-MatrixXd ths Eigen se double pinakes
void eigen_to_array(const VectorXd& vec, std::vector<double>& arr) {
    int n = vec.size();
    arr.resize(n);
    for(int i=0; i<n; i++) {
        arr[i] = vec(i);
    }
}   

VectorXd array_to_eigen(const std::vector<double>& arr, int n) {
    VectorXd v(n);
    for(int i=0; i<n; i++) v(i) = arr[i];
    return v;
}

// ---------------- SPARSE DIRECT SOLVER ---------------- 
VectorXd solve_sparse_direct(cs* A, const VectorXd& b_eig, bool use_spd) {
    int n = A->n;
    std::vector<double> x_val(n);
    std::vector<double> b_val(n);
    
    eigen_to_array(b_eig, b_val); 

    // Cholesky
    if (use_spd) {
        css* S = cs_schol(1, A); // to ordering
        csn* N = cs_chol(A, S); // ektelei chol  

        if (!N) { // an den einai thetika orismenos kane LU
            cs_sfree(S);
            return solve_sparse_direct(A, b_eig, false);
        }
        // Epilisi L kai U
        std::vector<double> y(n), x(n);
        
        cs_ipvec(S->pinv, b_val.data(), y.data(), n);   // y = P*b 
        cs_lsolve(N->L, y.data());                      // L*y = y 
        cs_ltsolve(N->L, y.data());                     // L'*x = y 
        cs_pvec(S->pinv, y.data(), x.data(), n);        // x = P'*y 

        cs_sfree(S);
        cs_nfree(N);
        return array_to_eigen(x, n);
    }
    // LU 
    else {
        // Ordering
        css* S = cs_sqr(2, A, 0);
        // LU
        csn* N = cs_lu(A, S, 1); 

        if (!N) { // fail LU
            cs_sfree(S);
            return VectorXd::Zero(n);
        }

        // Epilisi L kai U
        std::vector<double> y(n);
        cs_ipvec(N->pinv, b_val.data(), y.data(), n); // y = P*b 
        cs_lsolve(N->L, y.data());                    // L*y = y 
        cs_usolve(N->U, y.data());                    // U*x = y
        
        // Final Q 
        std::vector<double> final_x(n);
        cs_ipvec(S->q, y.data(), final_x.data(), n);  

        cs_sfree(S);
        cs_nfree(N);
        return array_to_eigen(final_x, n);
    }
}

// ---------------- SPARSE HELPERS ----------------
//  y = A*x me cs_gaxpy
VectorXd sparse_mat_vec(cs* A, const VectorXd& x) {
    int n = A->n;
    VectorXd y = VectorXd::Zero(n); // y se 0 arxikopoihsh
    // Arxikopoihseis kai metatropes se vectors
    std::vector<double> x_vec(n);
    std::vector<double> y_vec(n, 0.0); // y se 0 arxikopoihsh
    for(int i=0; i<n; ++i) {
        x_vec[i] = x(i);
    }

    cs_gaxpy(A, x_vec.data(), y_vec.data()); 

    for(int i=0; i<n; ++i) {
        y(i) = y_vec[i];
    }

    return y;
}

// Eyresi diagwniwn gia Jacobi preconditioner M = diag(A)
VectorXd get_sparse_diag_precond(cs* A) {
    int n = A->n;
    VectorXd M_inv = VectorXd::Ones(n); // Default 1 if missing
    
    for (int j = 0; j < n; j++) {
        for (int p = A->p[j]; p < A->p[j+1]; p++) {
             if (A->i[p] == j) { // found diagonal
                double val = A->x[p]; 
                if (std::abs(val) > 1e-16) {
                    M_inv(j) = 1.0 / val;
                }
            }
        }
    }
    return M_inv;
}

// ---------------- SPARSE ITERATIVE SOLVERS ---------------- 
VectorXd solve_sparse_cg(cs* A, const VectorXd& b, double itol, int max_iter) {
    int n = A->n;
    VectorXd x = VectorXd::Zero(n);
    VectorXd r = b - sparse_mat_vec(A, x); // r = b-A*x, GIA A*X
    
    double b_norm = b.norm();
    if (b_norm == 0) b_norm = 1.0;

    VectorXd M_inv = get_sparse_diag_precond(A); // Gia JACOBI PRECONDITIONER
    VectorXd z = r.cwiseProduct(M_inv); 
    
    VectorXd p = z; 
    double rho = r.dot(z); 
    
    int iter = 0;
    while ((r.norm() / b_norm > itol) && (iter < max_iter)) {
        iter++;
        
        VectorXd q = sparse_mat_vec(A, p); //A*p
        
        double alpha = rho / p.dot(q); 
        
        x = x + alpha * p; 
        r = r - alpha * q; 
        
        VectorXd z_new = r.cwiseProduct(M_inv); 
        double rho_new = r.dot(z_new); 
        
        double beta = rho_new / rho; 
        p = z_new + beta * p; 
        
        rho = rho_new;
    }
    return x;
}

VectorXd solve_sparse_bicg(cs* A, const VectorXd& b, double itol, int max_iter) {
    cs* At = cs_transpose(A, 1); // KANOYME KATASKEYH TOY A^t,gia na apofygoume ta for loops

    int n = A->n;
    VectorXd x = VectorXd::Zero(n);
    VectorXd r = b - sparse_mat_vec(A, x);
    VectorXd r_tilde = r; 
    
    double b_norm = b.norm();
    if (b_norm == 0) b_norm = 1.0;

    VectorXd M_inv = get_sparse_diag_precond(A);
    
    VectorXd z = r.cwiseProduct(M_inv);
    VectorXd z_tilde = r_tilde.cwiseProduct(M_inv);
    
    double rho = r_tilde.dot(z);
    
    VectorXd p = z;
    VectorXd p_tilde = z_tilde; 
    
    int iter = 0;
    while ((r.norm() / b_norm > itol) && (iter < max_iter)) {
        iter++;
        
        VectorXd q = sparse_mat_vec(A, p); // A * p
        VectorXd q_tilde = sparse_mat_vec(At, p_tilde); // A^T * p_tilde
        
        double omega = p_tilde.dot(q);
        if (std::abs(omega) < 1e-14) 
            break; 
        
        double alpha = rho / omega;
        
        x = x + alpha * p;
        r = r - alpha * q;
        r_tilde = r_tilde - alpha * q_tilde;
        
        VectorXd z_new = r.cwiseProduct(M_inv);
        VectorXd z_tilde_new = r_tilde.cwiseProduct(M_inv);
        
        double rho_new = r_tilde.dot(z_new);
        if (std::abs(rho_new) < 1e-14) 
            break;
        
        double beta = rho_new / rho;
        
        p = z_new + beta * p;
        p_tilde = z_tilde_new + beta * p_tilde;
        
        rho = rho_new;
    }
    cs_spfree(At); 
    return x;
}



// ----------------------- SOLVER ----------------------------------------
SolveResult solve_system(element* head, int num_nodes, const RunOptions& opts, std::variant<MatrixXd, cs*> A_var, VectorXd rhs) { 
    // me protsthiki varient gia na pernei kai tous 2 pinakes !!!

    int m2 = 0;
    for (element* e=head; e; e=e->next) {
        if (e->type == element::V || e->type == element::L) {
            m2++;
        }
    }

    VectorXd x;

    if (std::holds_alternative<cs*>(A_var)) { // elenxos poion pinaka exoyme
        cs* A = std::get<cs*>(A_var); // an petyxe o elenxos, pairnoume to A ws cs*

        int dim = A->n;
        int max_iter = dim; 

        if (opts.use_iter) {
            if (opts.use_spd) {
                x = solve_sparse_cg(A, rhs, opts.itol, max_iter);
            } else {
                x = solve_sparse_bicg(A, rhs, opts.itol, max_iter);
            }
        } 
        else {
            x = solve_sparse_direct(A, rhs, opts.use_spd);
        }
    }
    else {
        MatrixXd A = std::get<MatrixXd>(A_var);
        int dim = A.rows();
        
        // ITERATIVE SOLVERS
        if (opts.use_iter) {
            int max_iter = dim; 
            
            if (opts.use_spd) {
                // CG
                x = solve_cg(A, rhs, opts.itol, max_iter);
            } else {
                //Bi CG
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
