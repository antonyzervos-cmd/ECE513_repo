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
    if (opts.use_custom) { // custom solve chol or lu
        if (opts.use_spd) {
            x = custom_solve_cholesky(A, rhs);
        }
        else {
            x = custom_solve_lu(A, rhs);
        }
    }
    else {
        if (opts.use_spd) { // CHOLESKY EIGEN 
            Eigen::LLT<MatrixXd> LL(A);
            if (LL.info() == Eigen::Success) // Symmetrikos kai thetikos pinakas
                x = LL.solve(rhs);
            else    // eigen lu
                x = PartialPivLU<MatrixXd>(A).solve(rhs);
        }
        else { // eigen lu
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
