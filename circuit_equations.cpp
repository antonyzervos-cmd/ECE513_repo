#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <iomanip> 

using namespace std;
using namespace Eigen;

void build_MNA_DC(element* head, int num_nodes ) {
    int n = num_nodes;
    if (n < 0) {
        return;
    }

    // find how many V and L we have
    int m2 = 0;
    for (element* e = head; e != nullptr; e = e->next) {
        if (e->type == element::V || e->type == element::L) ++m2;
    }

    cout << "MNA-DC -- nodes (excluding ground): " << n << ", voltage sources (incl. L->0V): " << m2 << "\n";

    // create MNA blocks
    MatrixXd G = MatrixXd::Zero(n,n);    // array with resistors
    MatrixXd B = MatrixXd::Zero(n,m2);   // array with voltages
    VectorXd b = VectorXd::Zero(n);      // array with current sources
    VectorXd s = VectorXd::Zero(m2);     // array with voltage sources         


    int v_idx = 0;
    for (element* e = head; e != nullptr; e = e->next) {
        int n1 = e->nodes[0];
        int n2 = e->nodes[1];

        // epeidh o pinakas jekinaei apo 8esh (0,0) bazoyme ola ta nodes -1 gia na kanoun fit ston pinaka
        int i = (n1 == 0) ? -1 : (n1 - 1);
        int j = (n2 == 0) ? -1 : (n2 - 1);

        switch (e->type) {
            case element::R: {
                double R = static_cast<double>(e->value);
                double g = 1.0 / R;
                if (i >= 0) G(i, i) += g;
                if (j >= 0) G(j, j) += g;
                if (i >= 0 && j >= 0 && i != j) {
                    G(i, j) -= g;
                    G(j, i) -= g;
                }
                break;
            }

            case element::I: {
                double I = static_cast<double>(e->value);
                if (i >= 0) b(i) -= I; // current leaving node n1
                if (j >= 0) b(j) += I; // current entering node n2
                break;
            }

            case element::V: {
                if (i >= 0) B(i, v_idx) = 1.0;
                if (j >= 0) B(j, v_idx) = -1.0;
                s(v_idx) = static_cast<double>(e->value);
                ++v_idx;
                break;
            }

            case element::L: {
                // replace L with V with 0 value
                if (i >= 0) B(i, v_idx) = 1.0;
                if (j >= 0) B(j, v_idx) = -1.0;
                s(v_idx) = 0.0;
                ++v_idx;
                break;
            }

            case element::C: {
                // capacitor in DC -> open circuit
                break;
            }

            default:
                // D, Q, M ignored
                break;
        }
    }

    // create full MNA matrix
    int dim = n + m2;
    MatrixXd A = MatrixXd::Zero(max(0, dim), max(0, dim));
    A.topLeftCorner(n, n) = G;
    if (m2 > 0) A.topRightCorner(n, m2) = B;
    if (m2 > 0) A.bottomLeftCorner(m2, n) = B.transpose(); // anastrofos


    VectorXd rhs = VectorXd::Zero(max(0, dim));
    if (dim > 0) {
        if (n > 0 && m2 > 0) {
            // stack [b; s]
            rhs << b, s;
        } else if (n > 0 && m2 == 0) {
            rhs << b;
        }
    }

    // print the matrix
    cout << fixed << setprecision(3);
    cout << "\nAssembled MNA matrix A (" << dim << "x" << dim << "):\n" << A << "\n";
    cout << "\nRHS vector:\n" << rhs << "\n";

    if (dim == 0) {
        cout << "No unknowns to solve (dim == 0).\n";
        return;
    }

    // Solve system Ax = rhs
    VectorXd x;
    x = A.colPivHouseholderQr().solve(rhs);
    
    // Print node voltages
    cout << "\nNode voltages (node 0 = ground):\n";
    for (int k = 0; k < n; ++k) {
        cout << "v" << (k+1) << " = " << x(k) << " V\n";
    }

    // Print currents through voltage sources
    for (int k = 0; k < m2; ++k) {
        cout << "iV" << (k+1) << " = " << x(n + k) << " A\n";
    }
}
