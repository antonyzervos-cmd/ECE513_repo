#include <iostream>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <iomanip> 

using namespace std;
using namespace Eigen;

std::tuple<MatrixXd, VectorXd,  unordered_map<string,int>, unordered_map<string,int>> build_MNA_DC(element* head, int num_nodes) 
{
    int n = num_nodes;

    // count voltage sources and inductors (treated as V sources in DC)
    int m2 = 0;
    for (element* e = head; e != nullptr; e = e->next) {
        if (e->type == element::V || e->type == element::L)
            ++m2;
    }

    // Create MNA blocks
    MatrixXd G = MatrixXd::Zero(n, n);
    MatrixXd B = MatrixXd::Zero(n, m2);
    VectorXd b = VectorXd::Zero(n);
    VectorXd s = VectorXd::Zero(m2);

    // Maps
    unordered_map<string,int> vsrc_index_map;  // Voltage & inductors
    unordered_map<string,int> isrc_index_map;  // Current sources

    int v_idx = 0;

    for (element* e = head; e != nullptr; e = e->next)
    {
        int n1 = e->nodes[0];
        int n2 = e->nodes[1];

        int i = (n1 == 0 ? -1 : n1 - 1);
        int j = (n2 == 0 ? -1 : n2 - 1);

        switch (e->type)
        {
            // ----------------------------------------------------
            case element::R:
            {
                double R = static_cast<double>(e->value);
                double g = 1.0 / R;

                if (i >= 0) G(i,i) += g;
                if (j >= 0) G(j,j) += g;
                if (i >= 0 && j >= 0) {
                    G(i,j) -= g;
                    G(j,i) -= g;
                }
                break;
            }

            // ----------------------------------------------------
            case element::I:
            {
                double I = static_cast<double>(e->value);

                // stamp into b[]
                if (i >= 0) b(i) -= I;
                if (j >= 0) b(j) += I;

                // add to map
                isrc_index_map[e->name] = isrc_index_map.size();
                break;
            }

            // ----------------------------------------------------
            case element::V:
            {
                if (i >= 0) B(i, v_idx) = 1.0;
                if (j >= 0) B(j, v_idx) = -1.0;

                s(v_idx) = static_cast<double>(e->value);
                vsrc_index_map[e->name] = v_idx;

                ++v_idx;
                break;
            }

            // ----------------------------------------------------
            case element::L:
            {
                // DC: inductor -> short (behaves like ideal Vsource)
                if (i >= 0) B(i, v_idx) = 1.0;
                if (j >= 0) B(j, v_idx) = -1.0;

                s(v_idx) = 0.0; // DC: inductor has 0V
                vsrc_index_map[e->name] = v_idx;

                ++v_idx;
                break;
            }

            case element::C:
                // DC: capacitor , open , no stamping
                break;

            default:
                // ignore nonlinear DC elements (D, Q, M)
                break;
        }
    }

    // Build full matrix A   
    int dim = n + m2;
    MatrixXd A = MatrixXd::Zero(dim, dim);

    A.topLeftCorner(n, n) = G;

    if (m2 > 0) {
        A.topRightCorner(n, m2) = B;
        A.bottomLeftCorner(m2, n) = B.transpose();
    }


    // RHS = [ b ; s ]
    VectorXd rhs(dim);
    rhs << b, s;

    // Debug prints
    cout << fixed << setprecision(6);
    cout << "\n--- Assembled MNA matrix A (" << dim << "x" << dim << ") ---\n";
    cout << A << "\n";

    cout << "\n--- RHS vector ---\n";
    cout << rhs << "\n";


    // return A, rhs, and two maps
    return { A, rhs, vsrc_index_map, isrc_index_map };
}
