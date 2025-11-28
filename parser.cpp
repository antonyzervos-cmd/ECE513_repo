// --------------------------------
// PART 1 - PARSING 
// --------------------------------
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <utility>
using namespace std;

struct element {
    std::string name, model_name; // element name     
    enum type {R, C, L, V, I, D, Q, M} type; // element type - enum kai pedio mazi
    float value, length, width, area; // element values
    std::vector<int> nodes; // nodes vector 
    element* next = nullptr; // pointer to next element
};


// -----------------------------------------------------
// HASH Table (string -> int) helping function, return node id
// -----------------------------------------------------
int get_the_node_id(std::unordered_map<std::string, int>& node_map, int& node_next_id, const std::string& node_string) {
    // "0" = gnd
    if (node_string == "0") {
        return 0;
    }
    else if (node_map.find(node_string) == node_map.end()) { // Find gia apofigi collisions
        // Ean den vrethike, prosthetei neo node
        node_map[node_string] = node_next_id; // new node ID
        node_next_id++; 
    }
    return (node_map[node_string]); // return node id
}

// -----------------------------------------------------
// PRINT HASH TABLE
// -----------------------------------------------------    
void print_hash_table(const std::unordered_map<std::string, int>& node_map) {
    std::cout << "\n HASH TABLE (Node name, Node ID):\n";

    // Print contents key,value
    std::unordered_map<std::string, int>::const_iterator iterator; // iterator se kathe zevgos
    for (iterator = node_map.begin(); iterator != node_map.end(); ++iterator) {
        std::string name = iterator->first;   // name-string node
        int id = iterator->second;            // ID node
        std::cout << " ( NAME = " << name << ", ID = " << id << ") \n";
    }
}

// -----------------------------------------------------
// Convert string to lowercase
// -----------------------------------------------------
void to_lowercase(std::string& str) {
    std::transform(str.begin(), str.end(), str.begin(), // deikths arxh , deikths telos kai deikths overwrite sto palio to neo
                   [](unsigned char character){ 
                    return std::tolower(character); 
                });
}

// -----------------------------------------------------
// Check if string is a number
// -----------------------------------------------------
bool is_number(const std::string& string) {
    // check ean ena string mporei na ginei number
    char* end = nullptr;
    std::strtod(string.c_str(), &end);
    return end != string.c_str() && *end == '\0';
}

// -----------------------------------------------------
// ADD TO LINKED LIST
// -----------------------------------------------------
void add_element_to_list(element*& head, const element& e) {
    element* copied_element = new element(e); //copy constructor tou dosmenou element sto copied_element
    copied_element->next = nullptr;

    if (head == nullptr) { // ean h lista einai adeia
        head = copied_element;
    } 
    else { // hdh yparxei node, opote prosthese sto telos
        element* temporary = head;
        // Psakse last stoixeio list
        while (temporary->next != nullptr) { // check ean exei ayto to node epomeno stoixeio syndedemeno
            temporary = temporary->next; // an exei proxwra se epomeno
        }
        temporary->next = copied_element; // sto teleutaio stoixeio prosthese to neo
    }
}


// -----------------------------------------------------
// Print List
// -----------------------------------------------------
void print_the_list(element* head) {
    std::cout << "\n LIST: \n";
    element* temp = head; // temp gia traversal
    while (temp != nullptr) {
        std::cout << "NAME:" << temp->name << " TYPE: ";
        switch (temp->type) {
            case element::R: std::cout << "R"; break;
            case element::C: std::cout << "C"; break;
            case element::L: std::cout << "L"; break;
            case element::V: std::cout << "V"; break;
            case element::I: std::cout << "I"; break;
            case element::D: std::cout << "D"; break;
            case element::Q: std::cout << "Q"; break;
            case element::M: std::cout << "M"; break;
        }

        std::cout << " NODES: ";
        for (int node : temp->nodes) {
            std::cout << node << " ";
        }

        if (temp->type != element::M && temp->type != element::Q && temp->type != element::D) { // to M kai to Q kai to D den exoun value
            std::cout << " VALUE: " << temp->value;
        }

        if (!temp->model_name.empty()) { // ean exei model name 
            std::cout << " MODEL: " << temp->model_name;
        }
        // IF MOS, print L and W
        if (temp->type == element::M) {
            std::cout << " L= " << temp->length << " W= " << temp->width;
        }

        temp = temp->next; // move to next element
        std::cout << "\n";
    }
}


std::tuple<element*, int, std::unordered_map<std::string,int>> parse_netlist(const std::string& filename) { 
    try {
        // Open the file
        std::ifstream file(filename); // stream apo file
        file.exceptions(std::ifstream::badbit); // exceptions enable 
        std::cout << "File opened !!!\n";

        // -------------------------------------
        // HASH Table (string TO int)
        // -------------------------------------
        std::unordered_map<std::string, int> node_map; 
        int node_next_id = 1;  // Τo 0 einai gia gnd, opote ksekinaw apo 1, ara ola thetika ints (to orizw sthn main, gia na menei i timh)

        // -------------------------------------
        // Read File
        // -------------------------------------
        std::string line_buf; // line buffer
        element* head = nullptr; // start linked list

        while (std::getline(file, line_buf)) {
            // Agnohse kenes, sxolia, entoles elenxou
            if (line_buf[0] == '*' || line_buf[0] == '.') {
                continue;
            }
            if (line_buf.empty()) {
                continue;
            }
            
            // Metatrepse se Lowercase, Non case sensitive
            to_lowercase(line_buf);

            // Pare ta TOKENS   
            std::istringstream instrstream(line_buf); // stream apo string, anti apo input
            std::string token; // variable gia kathe token
            std::vector<std::string> tokens_vector; // vector me ola ta tokens 

            while (instrstream >> token) { // Diavazei apo to stream mia mia tis lekseis, agnoontas kena 
                tokens_vector.push_back(token);
                // Aytomato EOF
            }
            if (tokens_vector.empty()) {
                // Agnohse kenes grammes
                continue;
            }

            // Dhmiourgia neou stoixeiou struct, to opoio tha paei sto list
            // ARXIKOPOIHSH
            element e;
            e.name = tokens_vector[0]; // px V1
            e.next = nullptr;
            e.value = 0.0;
            e.length = 0.0;
            e.width = 0.0;  
            e.area = 1.0; 
            e.model_name = "";

            // Find element type apo prwto stoixeio onomatos
            char type_token = e.name[0];
            int num_nodes = 0; // ton arithmo twn komvwn gia to current line buffer

            switch (type_token) {
                case 'r': 
                    e.type = element::R;
                    num_nodes = 2; 
                    break;
                case 'c': 
                    e.type = element::C;
                    num_nodes = 2; 
                    break;
                case 'l': 
                    e.type = element::L;
                    num_nodes = 2; 
                    break;
                case 'v': 
                    e.type = element::V;
                    num_nodes = 2; 
                    break;
                case 'i': 
                    e.type = element::I;
                    num_nodes = 2; 
                    break;
                case 'd': 
                    e.type = element::D; 
                    num_nodes = 2;
                    break;
                case 'q': 
                    e.type = element::Q; 
                    num_nodes = 3;
                    break;
                case 'm': 
                    e.type = element::M; 
                    num_nodes = 4;
                    break;
                default:
                    std::cerr << "Wrong element type" << "\n";
                    continue;
            }

            // Update E.NODES
            for (int i = 1; i <= num_nodes; ++i) { // ksekinaei apo thesi vector 1, giati to 0 einai to onoma 
                // calling the hash table function
                int current_node_id = get_the_node_id(node_map, node_next_id, tokens_vector[i]);
                e.nodes.push_back(current_node_id); // prostithetai sto e.nodes
            }

            // Next token is VALUE or MODEL NAME
            std::string last_token = tokens_vector[num_nodes + 1];

            if (is_number(last_token)) {
                e.value = std::stof(last_token); // string to float , ean px 1e-6
            } 
            else {
                e.model_name = last_token;
            }

            // If we have MOS TRANSISTOR, read additional parameters L,W
            if (e.type == element::M) {
                for (size_t j = 6; j <= 7; ++j) { // j = 6, giati eimaste sthn 6 thesi tou tokens_vector
                    std::string position = tokens_vector[j]; // px l=1e-6, to fortwnw se variable 
                    size_t equals_position = position.find('='); //  vriskw pou einai to =
                    // pairnw meros PRIN to =
                    std::string l_or_w = position.substr(0, equals_position); 
                    to_lowercase(l_or_w);
                    // pairnw meros META to =
                    std::string value = position.substr(equals_position + 1); 
                    // apothikeuw sto struct
                    if (l_or_w == "l") {
                        e.length = std::stof(value);
                    }
                    else if (l_or_w == "w") {
                        e.width = std::stof(value);
                    }
                }
            }
            // Add to list
            add_element_to_list(head, e);
        }
        file.close();

        // -------------------------------------
        // Print Hash Table
        // -------------------------------------
        print_hash_table(node_map);

        return {head, static_cast<int>(node_map.size()), node_map};
    } 
    catch (const std::ios_base::failure& e) {
        std::cerr << "Cannot open given file: " << e.what() << "\n";
        return {nullptr, 0, {}};
    }
}