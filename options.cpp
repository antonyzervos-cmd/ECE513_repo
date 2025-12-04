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

    bool use_iter = false;     // iterative
    double itol = 1e-3;

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
    std::ifstream file(filename);
    RunOptions opts;

    opts.use_iter = false;
    opts.use_spd = false;
    opts.itol = 1e-3;
    opts.do_dc_sweep = false;

    std::string line;
    while (std::getline(file, line))
    {
        std::string low = line;
        to_lowercase(low);

        if (low.rfind(".options", 0) == 0) {
            std::istringstream iss(low);
            std::string token;
            iss >> token; // .options

            while (iss >> token) {
                if (token == "iter") 
                    opts.use_iter = true;
                else if (token == "spd") {
                    opts.use_spd = true;
                }
                else if (token.find("itol=") == 0)
                    opts.itol = std::stod(token.substr(5));
            }
        }
        else if (low.rfind(".print", 0) == 0 || low.rfind(".plot", 0) == 0)
        {
            std::istringstream iss(low);
            std::string token;
            iss >> token; // .print
            while (iss >> token)
                opts.plot_nodes.push_back(token);
        }
        else if (low.rfind(".dc", 0) == 0)
        {
            std::istringstream iss(low);
            std::string token;
            iss >> token;          // .dc
            iss >> opts.sweep_source;
            iss >> opts.sweep_start;
            iss >> opts.sweep_end;
            iss >> opts.sweep_step;

            opts.do_dc_sweep = true;
        }
    }

    return opts;
}
