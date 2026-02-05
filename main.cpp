#include <iostream>
#include <variant>
#include <suitesparse/cs.h>
#include "parser.cpp"
#include "circuit_equations.cpp"
#include "options.cpp"
#include "solver.cpp"
#include "dc_sweep.cpp"
#include "transient.cpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "No args! ./spice \n";
        return 1;
    }

    std::string filename = argv[1];

    // Parse netlist
    auto [head, num_nodes, node_map] = parse_netlist(filename);
    // Parse .options, .dc, .plot
    RunOptions opts = parse_options_from_file(filename);
    // Print 
    print_the_list(head);

    std::variant<MatrixXd, cs*> G_container; // G sparse or dense 
    std::variant<MatrixXd, cs*> C_container; // C 
    VectorXd b_dc; // DC RHS

    
    VectorXd rhs;
    std::unordered_map<std::string, int> vsrc_index_map;
    std::unordered_map<std::string, IStamp> isrc_index_map;

    if (opts.use_sparse) {
        // SPARSE MODE ON
        
        auto [G_sp, rhs_local, v_m, i_m] = build_MNA_sparse(head, num_nodes); // Build G 
        
        cs* C_sp = build_MNA_C_Sparse(head, num_nodes, v_m); // Build C
        
        G_container = G_sp;
        C_container = C_sp;
        b_dc = rhs_local;
        vsrc_index_map = v_m; 
        isrc_index_map = i_m;
    } 
    else {
        // DENSE MODE ON
        auto [G_dn, rhs_local, v_m, i_m] = build_MNA_DC(head, num_nodes); // Build G 
        MatrixXd C_dn = build_MNA_C_Dense(head, num_nodes, v_m); // Build C
        
        G_container = G_dn;
        C_container = C_dn;
        b_dc = rhs_local;
        vsrc_index_map = v_m; 
        isrc_index_map = i_m;
    }

    // DC OP is calculated using t=0 transient values, arxiko RHS
    if (opts.do_transient) {
        b_dc.setZero(); 
        for (element* e = head; e; e = e->next) {
            if (e->type == element::I || e->type == element::V) {
                double val = eval_transient_src(e, 0.0);
                
                if (e->type == element::I) {
                    if (isrc_index_map.find(e->name) != isrc_index_map.end()) {
                        IStamp st = isrc_index_map.at(e->name);
                        if (st.i >= 0) b_dc(st.i) -= val;
                        if (st.j >= 0) b_dc(st.j) += val;
                    }
                }
                else if (e->type == element::V) {
                    if (vsrc_index_map.find(e->name) != vsrc_index_map.end()) {
                        int k = vsrc_index_map.at(e->name);
                        b_dc(num_nodes + k) = val;
                    }
                }
            }
        }
    }

    // Solve DC Operating Point
    SolveResult sol = solve_system(head, num_nodes, opts, G_container, b_dc);
    write_dc_op("dc_op.txt", head, sol);


    // Transient Analysis
    if (opts.do_transient) {
        run_transient_analysis(head, num_nodes, opts, G_container, C_container, sol.x, vsrc_index_map, isrc_index_map, node_map);
    }

    // DC Sweep
    if (opts.do_dc_sweep) {
        VectorXd rhs_backup = b_dc; 
        for (const auto& cmd : opts.dc_commands) {
            rhs = rhs_backup; // Reset RHS
            run_dc_sweep(head, num_nodes, cmd, opts, G_container, rhs, vsrc_index_map, isrc_index_map, node_map);
        }
    }


    // Cleanup
    if (opts.use_sparse) {
        if(std::holds_alternative<cs*>(G_container)) cs_spfree(std::get<cs*>(G_container));
        if(std::holds_alternative<cs*>(C_container)) cs_spfree(std::get<cs*>(C_container));
    }

    // Free linked list elements
    while (head) {
        element* t = head;
        head = head->next;
        delete t;
    }

    return 0;
}