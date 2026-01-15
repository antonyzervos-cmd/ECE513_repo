#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <variant>
#include <fstream>
#include <Eigen/Dense>
#include <suitesparse/cs.h>

using namespace Eigen;
using namespace std;

// ypologizei phgh timhs V/I SE ENA t
double eval_transient_src(element* e, double t) {
    const TransientSpec& spec = e->tran_spec;

    if (spec.type == NONE) 
        return e->value; // DC PERIPTWSH

    const auto& p = spec.params;
    
    if (spec.type == EXP) { 
        double i1=p[0], i2=p[1], td1=p[2], tc1=p[3], td2=p[4], tc2=p[5];
        if (t <= td1) 
            return i1;
        if (t <= td2) 
            return i1 + (i2-i1)*(1 - exp(-(t-td1)/tc1));
        return i1 + (i2-i1)*(exp(-(t-td2)/tc2) - exp(-(t-td1)/tc1)); 
    }
    else if (spec.type == SIN) { 
        double i1=p[0], ia=p[1], fr=p[2], td=p[3], df=p[4], ph=p[5];
        if (t <= td)
            return i1 + ia*sin(2*M_PI*ph/360.0);
        return i1 + ia * sin(2*M_PI*fr*(t-td) + 2*M_PI*ph/360.0) * exp(-(t-td)*df);
    }
    else if (spec.type == PULSE) { 
        double i1=p[0], i2=p[1], td=p[2], tr=p[3], tf=p[4], pw=p[5], per=p[6];
        if (t < td) 
            return i1;
        
        double t_cycle = fmod(t - td, per);
        if (t_cycle < tr) 
            return i1 + (i2-i1)*(t_cycle/tr);
        if (t_cycle < tr+pw) 
            return i2;
        if (t_cycle < tr+pw+tf) 
            return i2 - (i2-i1)*((t_cycle - tr - pw)/tf);
        return i1;
    }
    else if (spec.type == PWL) {
        const auto& pts = spec.pwl_points;
        if (pts.empty()) 
            return e->value;
        if (t <= pts.front().first) 
            return pts.front().second;
        if (t >= pts.back().first) 
            return pts.back().second;
        
        // Linear Interpolation
        for (size_t i = 0; i < pts.size()-1; ++i) {
            if (t >= pts[i].first && t <= pts[i+1].first) {
                double t1 = pts[i].first;
                double t2 = pts[i+1].first;
                double v1 = pts[i].second;
                double v2 = pts[i+1].second;
                return v1 + (v2-v1)*(t-t1)/(t2-t1);
            }
        }
    }
    return e->value;
}

// update RHS e(t) at time t
VectorXd update_rhs_transient(element* head, int n, int m2, double t,
                              const unordered_map<string,int>& vsrc_map,
                              const unordered_map<string, IStamp>& isrc_map) {
    
    VectorXd b = VectorXd::Zero(n);
    VectorXd s = VectorXd::Zero(m2);
    
    // vlepw ola ta stoixeia sth list
    for (element* e = head; e; e = e->next) {
        if (e->type == element::I || e->type == element::V) {
            double val = eval_transient_src(e, t);
            
            // enimeroso RHS
            if (e->type == element::I) {
                IStamp st = isrc_map.at(e->name);
                if (st.i >= 0) b(st.i) -= val;
                if (st.j >= 0) b(st.j) += val;
            }
            else { // V
                 int k = vsrc_map.at(e->name);
                 s(k) = val;
            }
        }
    }
    
    VectorXd rhs(n + m2);
    rhs << b, s;
    return rhs;
}

// Main Transient Loop
void run_transient_analysis(
    element* head, int num_nodes, const RunOptions& opts,
    std::variant<MatrixXd, cs*> G_var, // G 
    std::variant<MatrixXd, cs*> C_var, //  C - new
    VectorXd x_prev, // initial condition from DC OP
    const unordered_map<string,int>& vsrc_map,
    const unordered_map<string, IStamp>& isrc_map,
    const unordered_map<string,int>& node_map) {

    double h = opts.tran_cmd.step;
    double t_fin = opts.tran_cmd.fin_time;
    int m2 = vsrc_map.size();
    
    // ftiaxnei output
    vector<ofstream> files;
    for(auto& node : opts.tran_cmd.plot_nodes) {
        string name = node; 
        size_t p1 = name.find('('); size_t p2 = name.find(')');
        if(p1!=string::npos && p2!=string::npos) 
            name = name.substr(p1+1, p2-p1-1);
        files.emplace_back("tran_V(" + name + ").txt");
    }


    // TR: A = G + (2/h)*C
    // BE: A = G + (1/h)*C
    double coeff = (opts.tran_method == TR) ? (2.0/h) : (1.0/h);
    
    std::variant<MatrixXd, cs*> A_eff_var;
    
    // A_eff
    if (opts.use_sparse) {
        cs* G = std::get<cs*>(G_var);
        cs* C = std::get<cs*>(C_var);
        
        // A_eff = G + coeff * C
        cs* A_eff = cs_add(G, C, 1.0, coeff);
        A_eff_var = A_eff;
    } 
    else {
        MatrixXd G = std::get<MatrixXd>(G_var);
        MatrixXd C = std::get<MatrixXd>(C_var);
        MatrixXd A_eff = G + coeff * C;
        A_eff_var = A_eff;
    }

    // write sto file
    double t = 0;
    for(size_t i=0; i<opts.tran_cmd.plot_nodes.size(); ++i) {
        string node = opts.tran_cmd.plot_nodes[i];
         size_t p1 = node.find('('); size_t p2 = node.find(')');
        if(p1!=string::npos) {
            node = node.substr(p1+1, p2-p1-1);
        }
        transform(node.begin(), node.end(), node.begin(), ::tolower);
        int idx = node_map.at(node) - 1;
        files[i] << t << " " << x_prev(idx) << "\n";
    }

    // main loop time
    while (t < t_fin) {
        t += h;
        
        VectorXd e_t = update_rhs_transient(head, num_nodes, m2, t, vsrc_map, isrc_map);
        VectorXd rhs_final;

        if (opts.tran_method == BE) {
            // BEuler: rhs = e(t) + (1/h)*C*x_prev
            VectorXd Cx;
            if (opts.use_sparse) {
                Cx = sparse_mat_vec(std::get<cs*>(C_var), x_prev);
            }
            else {
                Cx = std::get<MatrixXd>(C_var) * x_prev;
            }
            
            rhs_final = e_t + coeff * Cx;
        } 
        else { 
            // TRap: rhs = e(t) + e(t-h) - (G - (2/h)C)*x_prev
           
            VectorXd e_old = update_rhs_transient(head, num_nodes, m2, t-h, vsrc_map, isrc_map);
            
            VectorXd Cx, Gx;
            if (opts.use_sparse) {
                Cx = sparse_mat_vec(std::get<cs*>(C_var), x_prev);
                Gx = sparse_mat_vec(std::get<cs*>(G_var), x_prev);
            } else {
                Cx = std::get<MatrixXd>(C_var) * x_prev;
                Gx = std::get<MatrixXd>(G_var) * x_prev;
            }
            
            rhs_final = e_t + e_old + coeff * Cx - Gx; 
        }

        // solve system A_eff * x = rhs_final
        SolveResult res = solve_system(head, num_nodes, opts, A_eff_var, rhs_final);
        x_prev = res.x;// i nea lysi ginetai palia gia next step

        // Plot
        for(size_t i=0; i<opts.tran_cmd.plot_nodes.size(); ++i) {
            string node = opts.tran_cmd.plot_nodes[i];
             size_t p1 = node.find('('); size_t p2 = node.find(')');
             if(p1!=string::npos) 
                node = node.substr(p1+1, p2-p1-1);
             transform(node.begin(), node.end(), node.begin(), ::tolower);
            int idx = node_map.at(node) - 1;
            files[i] << t << " " << x_prev(idx) << "\n";
        }
    }
    
    // free sparse EAN USE SPARSE
    if (opts.use_sparse) {
        cs_spfree(std::get<cs*>(A_eff_var));
    }
}