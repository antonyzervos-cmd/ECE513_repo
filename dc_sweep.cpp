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
    unordered_map<string, IStamp>& isrc_index_map,
    const unordered_map<string,int>& node_map) {

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


    // arxiko b_base gia dc_sweep
    VectorXd b_base = rhs.segment(0, num_nodes);

    // LOOP SWEEP
    for (double v = opts.sweep_start; v <= opts.sweep_end; v += opts.sweep_step)
    {
        //SWEEP VOLTAGE SOURCE
        if (sweep_is_voltage) {
            int k = vsrc_index_map[opts.sweep_source];
            int rhs_pos = num_nodes + k;   // s[]

            rhs(rhs_pos) = v;  // change s[k]

                // --- PRINT RHS ---
           // cout << "\n========== VOLTAGE SWEEP ==========\n";
           // cout << "Set " << opts.sweep_source << " = " << v << "\n";
           // cout << "Updated RHS:\n" << rhs << "\n";
           // cout << "===================================\n";
        }
        // SWEEP CURRENT SOURCE
        else if (sweep_is_current)
        {
            IStamp st = isrc_index_map[opts.sweep_source];

            rhs.segment(0, num_nodes) = b_base;

            if (st.i >= 0) rhs(st.i) -= v;
            if (st.j >= 0) rhs(st.j) += v;

                // --- PRINT RHS ---
             //   cout << "\n========== CURRENT SWEEP ==========\n";
            //    cout << "Set " << opts.sweep_source << " = " << v << "\n";
             //   cout << "Updated RHS:\n" << rhs << "\n";
            //    cout << "===================================\n";
        }

        // solution
        SolveResult sol = solve_system(head, num_nodes, opts, A, rhs);

        for (size_t i = 0; i < opts.plot_nodes.size(); ++i)
        {
            const string& tok = opts.plot_nodes[i];   // px "V(4)" ή "V(out)"

            size_t open = tok.find('(');
            size_t close = tok.find(')');

            string node_name = tok.substr(open + 1, close - open - 1);

            // lowercase 
            transform(node_name.begin(), node_name.end(),
                    node_name.begin(), ::tolower);


            int mna_id = node_map.at(node_name);   // map apo to HASH sto MNA ID

            double val = sol.x(mna_id - 1);
            outputs[i] << v << " " << val << "\n";
        }
    }
}
