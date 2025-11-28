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

std::vector<RunOptions> parse_all_options(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<RunOptions> all;

    RunOptions base; 
    base.do_dc_sweep = false;

    std::string line;
    while (std::getline(file, line))
    {
        std::string low = line;
        to_lowercase(low);

        // -----------------------
        // .print / .plot
        // -----------------------
        if (low.rfind(".print", 0) == 0 || low.rfind(".plot", 0) == 0)
        {
            std::istringstream iss(low);
            std::string cmd;
            iss >> cmd;         // .print
            std::string token;

            while (iss >> token)
                base.plot_nodes.push_back(token);

            continue;
        }

        // -----------------------
        // .dc
        // -----------------------
        if (low.rfind(".dc", 0) == 0)
        {
            RunOptions dc = base;  

            std::istringstream iss(low);
            std::string cmd;
            iss >> cmd;  // .dc

            iss >> dc.sweep_source;
            iss >> dc.sweep_start;
            iss >> dc.sweep_end;
            iss >> dc.sweep_step;

            dc.do_dc_sweep = true;

            all.push_back(dc);
            continue;
        }
    }

    if (all.empty())
        all.push_back(base);

    return all;
}