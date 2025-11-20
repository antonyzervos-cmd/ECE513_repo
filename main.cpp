// main.cpp
// g++ main.cpp -I /usr/include/eigen3 -lgsl -lgslcblas -O2 -std=c++17 -o main

#include <iostream>
#include "parser.cpp"
#include "circuit_equations.cpp"
#include "options.cpp"
#include "solver.cpp"
#include "dc_sweep.cpp"

//RunOptions parse_options_from_file(const std::string&);
//SolveResult solve_system(element*, int, const RunOptions&);
//void write_dc_op(const std::string&, element*, const SolveResult&);
//void run_dc_sweep(element*, int, const RunOptions&);
//void print_the_list(element*);

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "No args!\n";
        return 1;
    }

    std::string filename = argv[1];

    auto [head, num_nodes] = parse_netlist(filename);

    RunOptions opts = parse_options_from_file(filename);

    print_the_list(head);

    SolveResult sol = solve_system(head, num_nodes, opts);
    write_dc_op("dc_op.txt", head, sol);

    if (opts.do_dc_sweep) {
        run_dc_sweep(head, num_nodes, opts);
    }

    // free memory
    while (head) {
        element* t = head;
        head = head->next;
        delete t;
    }

    return 0;
}
