#include <iostream>
#include <variant>
#include <suitesparse/cs.h>
#include "parser.cpp"
#include "circuit_equations.cpp"
#include "options.cpp"
#include "solver.cpp"
#include "dc_sweep.cpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "No args! Usage: ./spice <netlist>\n";
        return 1;
    }

    std::string filename = argv[1];

    // Parse netlist
    auto [head, num_nodes, node_map] = parse_netlist(filename);
    // Parse .options, .dc, .plot
    RunOptions opts = parse_options_from_file(filename);
    // Print circuit elements 
    print_the_list(head);

  
    // Build System Sparse or Dense
    std::variant<MatrixXd, cs*> A_container;
    
    VectorXd rhs;
    std::unordered_map<std::string, int> vsrc_index_map;
    std::unordered_map<std::string, IStamp> isrc_index_map;

    if (opts.use_sparse) {

        auto [A_sp, b_sp, v_map, i_map] = build_MNA_sparse(head, num_nodes);
        
        A_container = A_sp;
        rhs = b_sp;
        vsrc_index_map = v_map;
        isrc_index_map = i_map;
    } 
    else {
        // Dense case
        auto [A_dn, b_dn, v_map, i_map] = build_MNA_DC(head, num_nodes);
        
        A_container = A_dn;
        rhs = b_dn;
        vsrc_index_map = v_map;
        isrc_index_map = i_map;
    }


    // Solve DC Operating Point
    SolveResult sol = solve_system(head, num_nodes, opts, A_container, rhs);
    write_dc_op("dc_op.txt", head, sol);


    // DC Sweep
    if (opts.do_dc_sweep) {
        VectorXd rhs_backup = rhs; 
        for (const auto& cmd : opts.dc_commands) {
            rhs = rhs_backup; // Reset RHS
            run_dc_sweep(head, num_nodes, cmd, opts, A_container, rhs, vsrc_index_map, isrc_index_map, node_map);
        }
    }


    // Cleanup
    if (opts.use_sparse && std::holds_alternative<cs*>(A_container)) {
        cout << "\n SPARSE \n";
        cs_spfree(std::get<cs*>(A_container));
    }

    // Free linked list elements
    while (head) {
        element* t = head;
        head = head->next;
        delete t;
    }

    return 0;
}