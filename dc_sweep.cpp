#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

void run_dc_sweep(
    element* head,
    int num_nodes,
    const RunOptions& opts,
    MatrixXd& A,
    VectorXd& rhs,
    unordered_map<string,int>& vsrc_index_map,
    unordered_map<string,int>& isrc_index_map)
{

    // check vsrc or isrc
    bool sweep_is_voltage = false;
    bool sweep_is_current = false;

    if (vsrc_index_map.find(opts.sweep_source) != vsrc_index_map.end())
        sweep_is_voltage = true;
    else if (isrc_index_map.find(opts.sweep_source) != isrc_index_map.end())
        sweep_is_current = true;
    else {
        cerr << "ERROR: sweep source not found \n" << opts.sweep_source;
        return;
    }

    // ta plot nodes
    vector<ofstream> outputs;
    for (const auto& p : opts.plot_nodes) {
        string fname = "dc_sweep_" + opts.sweep_source + "_V(" + p + ").txt";
        outputs.emplace_back(fname);
    }

    // LOOP SWEEP
    for (double v = opts.sweep_start; v <= opts.sweep_end; v += opts.sweep_step)
    {
        //SWEEP VOLTAGE SOURCE
        if (sweep_is_voltage) {
            int k = vsrc_index_map[opts.sweep_source];
            int rhs_pos = num_nodes + k;   // s[]

            rhs(rhs_pos) = v;  // change s[k]
        }
        // SWEEP CURRENT SOURCE
        else if (sweep_is_current) {
            // den exw 1 alla 2 indexes gia current source
            // vres element Isource wste na pareis tis times apo ta nodes parakatw
            element* src = nullptr;
            for (element* e = head; e; e = e->next) {
                string low = e->name;
                transform(low.begin(), low.end(), low.begin(), ::tolower);
                if (low == opts.sweep_source) {
                    src = e;
                    break;
                }
            }

            int n1 = src->nodes[0];
            int n2 = src->nodes[1];
            int i = (n1 == 0 ? -1 : n1 - 1);
            int j = (n2 == 0 ? -1 : n2 - 1);

            VectorXd b = rhs.segment(0, num_nodes); // new clean b

            // isource n1 to n2
            if (i >= 0) b(i) -= v;
            if (j >= 0) b(j) += v;

            //  kai update rhs
            rhs.segment(0, num_nodes) = b; 

        }

        // solution
        SolveResult sol = solve_system(head, num_nodes, opts, A, rhs);

        // sta arxeia grafw to output
        for (size_t i = 0; i < opts.plot_nodes.size(); ++i) {

            const string& tok = opts.plot_nodes[i];
            string digits;
            for (char c : tok)
                if (isdigit((unsigned char)c)) digits.push_back(c);

            int node = digits.empty() ? -1 : stoi(digits);
            double val = (node >= 1 && node <= sol.n_nodes)
                         ? sol.x(node - 1)
                         : 0;

            outputs[i] << v << " " << val << "\n";
        }
    }
}
