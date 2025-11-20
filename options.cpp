// options.cpp
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

struct RunOptions {
    bool use_spd = false;
    bool use_custom = false;
    bool do_dc_sweep = false;

    std::string sweep_source; // onoma phghs sarwshs
    double sweep_start = 0;
    double sweep_end = 0;
    double sweep_step = 0;

    std::vector<std::string> plot_nodes;  // nodes to plot 
};

static inline std::string to_lower_str(const std::string& s) {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    return t;
}

RunOptions parse_options_from_file(const std::string& filename) {
    RunOptions opts;
    std::ifstream ifs(filename);


    std::string line;
    while (std::getline(ifs, line)) {
        std::string trimmed = line; // trimmed = afairesh kenwn kai tabs
        trimmed.erase(0, trimmed.find_first_not_of(" \t")); // prwtos char pou den einai space h tab

        if (trimmed.empty()) 
            continue; // keni
        if (trimmed[0] == '*') 
            continue; // sxolio

        std::string low = to_lower_str(trimmed);

        if (low.rfind(".options", 0) == 0) { // .options sthn thesi 0 tou string
            std::istringstream iss(low);
            std::string dot, tok;
            iss >> dot;
            while (iss >> tok) {
                if (tok == "spd") //chol
                    opts.use_spd = true;
                else if (tok == "custom") //custom LU or chol
                    opts.use_custom = true;
            }
        }
        else if (low.rfind(".dc", 0) == 0) {
            std::istringstream iss(trimmed);
            std::string dot, var_sweep_source;
            double start, end, step;
            if (iss >> dot >> var_sweep_source >> start >> end >> step) { // check for 5 tokens
                opts.do_dc_sweep = true;
                opts.sweep_source = to_lower_str(var_sweep_source);
                opts.sweep_start = start;
                opts.sweep_end = end;
                opts.sweep_step = step;
            }
        }
        else if (low.rfind(".plot", 0) == 0 || low.rfind(".print", 0) == 0) {
            std::istringstream iss(trimmed);
            std::string cmd, tok;
            iss >> cmd; // skip .plot, pairnw to epomeno 
            while (iss >> tok) { // ola ta V
                if (!tok.empty() && tok.back() == ',') {
                    tok.pop_back();
                }
                opts.plot_nodes.push_back(tok);
            }
        }
    }

    return opts;
}
