#include <iostream>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <iomanip> 
#include <suitesparse/cs.h>

using namespace std;
using namespace Eigen;


struct IStamp {
    int i;
    int j;
};


std::tuple<MatrixXd, VectorXd,  unordered_map<string,int>, unordered_map<string, IStamp>> build_MNA_DC(element* head, int num_nodes) {
    int n = num_nodes;

    // count voltage sources and inductors - treated as V sources in DC
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
    unordered_map<string, IStamp> isrc_index_map;  // Current sources

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

                if (i >= 0) b(i) -= I;
                if (j >= 0) b(j) += I;

               // cout << i << j;

                isrc_index_map[e->name] = { i, j };
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
                
                if (i >= 0) B(i, v_idx) = 1.0;
                if (j >= 0) B(j, v_idx) = -1.0;

                s(v_idx) = 0.0; // DC: inductor has 0V
                vsrc_index_map[e->name] = v_idx;

                ++v_idx;
                break;
            }

            case element::C:
                // DC capacitor 
                break;

            default:
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


    // cout << fixed << setprecision(6);
    // cout << "\n--- Assembled MNA matrix A (" << dim << "x" << dim << ") ---\n";
    // cout << A << "\n";

    // cout << "\n--- RHS vector ---\n";
    // cout << rhs << "\n";


    // return A, rhs, and two maps
    return { A, rhs, vsrc_index_map, isrc_index_map};
}


std::tuple<cs*, VectorXd, unordered_map<string,int>, unordered_map<string, IStamp>> build_MNA_sparse(element* head, int num_nodes) {
    int n = num_nodes;

    int m2 = 0;
    for (element* e = head; e != nullptr; e = e->next) {
        if (e->type == element::V || e->type == element::L)
            ++m2;
    }

    int dim = n + m2;

   
    cs* T = cs_spalloc(dim, dim, 1, 1, 1); 
    // typika orismata to 1 gia plithos nz pou tha ginoun realloc meta

    VectorXd b = VectorXd::Zero(n);
    VectorXd s = VectorXd::Zero(m2);

    unordered_map<string,int> vsrc_index_map;
    unordered_map<string, IStamp> isrc_index_map;

    int v_idx = 0; 

    for (element* e = head; e != nullptr; e = e->next)
    {
        int n1 = e->nodes[0];
        int n2 = e->nodes[1];

        int i = (n1 == 0 ? -1 : n1 - 1);
        int j = (n2 == 0 ? -1 : n2 - 1);

        switch (e->type)
        {
            case element::R:
            {
                double g = 1.0 / static_cast<double>(e->value);
                // use cs_entry to add to triplet
                if (i >= 0) cs_entry(T, i, i, g);
                if (j >= 0) cs_entry(T, j, j, g);
                if (i >= 0 && j >= 0) {
                    cs_entry(T, i, j, -g);
                    cs_entry(T, j, i, -g);
                }
                break;
            }

            case element::I:
            {
                double I = static_cast<double>(e->value);
                if (i >= 0) b(i) -= I;
                if (j >= 0) b(j) += I;
                isrc_index_map[e->name] = { i, j };
                break;
            }

            case element::V: // idio me L
            case element::L:
            {
              
                int row_k = n + v_idx; 

                double val = (e->type == element::V) ? static_cast<double>(e->value) : 0.0;
                s(v_idx) = val;

                if (i >= 0) {
                    cs_entry(T, i, row_k, 1.0);  // B
                    cs_entry(T, row_k, i, 1.0);  // B^T
                }
                if (j >= 0) {
                    cs_entry(T, j, row_k, -1.0); // B
                    cs_entry(T, row_k, j, -1.0); // B^T
                }

                vsrc_index_map[e->name] = v_idx;
                ++v_idx;
                break;
            }
            default: break;
        }
    }

    cs* A = cs_compress(T); // TRIPLET SE COMPRESSED COLUMN
    cs_spfree(T); // free triplet

    cs_dupl(A); // ATHROISMA APO PRIN

    // RHS
    VectorXd rhs(dim);
    rhs << b, s;

    return { A, rhs, vsrc_index_map, isrc_index_map };
}


MatrixXd build_MNA_C_Dense(element* head, int num_nodes, 
                           const unordered_map<string,int>& vsrc_map) {
    int n = num_nodes;
    int m2 = vsrc_map.size(); 
    int dim = n + m2;
    
    MatrixXd C_mat = MatrixXd::Zero(dim, dim);

    for (element* e = head; e != nullptr; e = e->next) {
        int n1 = e->nodes[0];
        int n2 = e->nodes[1];
        int i = (n1 == 0 ? -1 : n1 - 1);
        int j = (n2 == 0 ? -1 : n2 - 1);

        if (e->type == element::C) {
            double C_val = static_cast<double>(e->value);
            if (i >= 0) C_mat(i,i) += C_val;
            if (j >= 0) C_mat(j,j) += C_val;
            if (i >= 0 && j >= 0) {
                C_mat(i,j) -= C_val;
                C_mat(j,i) -= C_val;
            }
        }
        else if (e->type == element::L) {
            double L_val = static_cast<double>(e->value);
            // We need its index from the map
            if (vsrc_map.find(e->name) != vsrc_map.end()) {
                int k = vsrc_map.at(e->name); 
                int row = n + k;
                C_mat(row, row) -= L_val; 
            }
        }
    }
    return C_mat;
}

cs* build_MNA_C_Sparse(element* head, int num_nodes, 
                       const unordered_map<string,int>& vsrc_map) {
    int n = num_nodes;
    int m2 = vsrc_map.size();
    int dim = n + m2;

    cs* T = cs_spalloc(dim, dim, 1, 1, 1);

    for (element* e = head; e != nullptr; e = e->next) {
        int n1 = e->nodes[0];
        int n2 = e->nodes[1];
        int i = (n1 == 0 ? -1 : n1 - 1);
        int j = (n2 == 0 ? -1 : n2 - 1);

        if (e->type == element::C) {
            double C_val = static_cast<double>(e->value);
            if (i >= 0) cs_entry(T, i, i, C_val);
            if (j >= 0) cs_entry(T, j, j, C_val);
            if (i >= 0 && j >= 0) {
                cs_entry(T, i, j, -C_val);
                cs_entry(T, j, i, -C_val);
            }
        }
        else if (e->type == element::L) {
            double L_val = static_cast<double>(e->value);
            if (vsrc_map.find(e->name) != vsrc_map.end()) {
                int k = vsrc_map.at(e->name);
                int row = n + k;
                cs_entry(T, row, row, -L_val);
            }
        }
    }
    
    cs* C = cs_compress(T);
    cs_spfree(T);
    cs_dupl(C);
    return C;
}