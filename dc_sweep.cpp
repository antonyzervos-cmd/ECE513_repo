// dc_sweep.cpp
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

//struct element;
//struct RunOptions;
//struct SolveResult;

//SolveResult solve_system(element* head, int num_nodes, const RunOptions&);
//void write_dc_op(const std::string&, element*, const SolveResult&);

// find source sarwshs
static element* find_src(element* head, const std::string& name) { 
    std::string low = name;
    std::transform(low.begin(), low.end(), low.begin(), ::tolower); // to lower

    for (element* e=head; e; e=e->next) { // check list to find source name
        std::string en = e->name;
        std::transform(en.begin(), en.end(), en.begin(), ::tolower);
        if (en == low)
            return e;
    }
    return nullptr; // didnt find
}

void run_dc_sweep(element* head, int num_nodes, const RunOptions& opts) {
    
    element* src = find_src(head, opts.sweep_source);
    if (!src) {
        std::cerr << "source not found\n";
        return;
    }

    std::vector<std::ofstream> outputs; 
    for (const auto& p : opts.plot_nodes) { // poia nodes tha ginoun plot
        std::string fname = "dc_sweep_" + opts.sweep_source + "_V(" + p + ").txt";
        outputs.emplace_back(fname);
    }

    for (double v = opts.sweep_start; v <= opts.sweep_end + 1e-12; v += opts.sweep_step) {
        src->value = v; // timi source

        SolveResult sol = solve_system(head, num_nodes, opts); // lysi gia kathe timi

        for (size_t i = 0 ; i < opts.plot_nodes.size(); ++i) { //gia kathe node pou theloume na kanoume plot, pairnoume 
            const std::string& tok = opts.plot_nodes[i];    // MONO TIN TIMH 
            std::string digits;
            for (char c : tok){
                if (isdigit((unsigned char)c)) 
                    digits.push_back(c);
            }

            int node = (digits.empty() ? -1 : stoi(digits)); // ean den yparxei arithmos, -1

            double val;
            if (node >= 1 && node <= sol.n_nodes) { // check ekgyrou node
               val = sol.x(node-1);
            }
            else {
                val = 0;
            }
        
            outputs[i] << v << " " << val << "\n";
        }
    }
}
