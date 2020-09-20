// defined std::unique_ptr
#include <memory>
// defines Var and Lit
#include "minisat/core/SolverTypes.h"
// defines Solver
#include "minisat/core/Solver.h"


#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>
#include "prjece650.h"
#include <pthread.h>
#include <time.h>
//#include <errno.h>
//#include <error.h>

// store result of mini vertex cover
std::vector<int> vertex_cover_appro1 = {};
std::vector<int> vertex_cover_appro2 = {};
std::vector<int> vertex_cover_sat = {};

std::vector <Edge> edge_list_appro1 = {};
std::vector <Edge> edge_list_appro2 = {};
std::vector <Edge> edge_list_sat = {};
std::vector <Edge> edge_list = {};

std::vector <std::vector<int >> adj_list;

int n = 0;   //vertices number

// declare thread
//    pthread_t thread_io;
pthread_t thread_sat;
pthread_t thread_appro1;
pthread_t thread_appro2;

clockid_t clock_sat;
clockid_t clock_appro1;
clockid_t clock_appro2;

struct timespec ts_sat;
struct timespec ts_appro1;
struct timespec ts_appro2;

// input string
std::string command;
int vertices = 0;

std::string result_appro1;
std::string result_appro2;
std::string result_sat;

bool isValidVerticesCommandFormat(const std::string command) {
    std::regex re(R"(^\s*V\s*[0-9]+\s*$)");
    if (!std::regex_match(command, re)) {
        std::cerr << "Error: Vertices Command Error, should be 'E <number>' only.\n";
        std::cerr.flush();
        return false;
    }
    return true;
}

bool isValidEdgeCommandFormat(const std::string command) {
    std::regex re(R"(^\s*E\s*\{(<[0-9]+,[0-9]+>,*){0,}\}\s*$)");
    if (!std::regex_match(command, re)) {
        std::cerr << "Error: Edges Command Error, should be 'E "
                     "{<number,number><number,number>..}' only.\n";
        std::cerr.flush();
        return false;
    }
    return true;
}


int getVertices(const std::string v_command) {
    std::regex re(R"(^\s*V\s*([0-9]+)\s*$)");
    std::smatch m;
    std::regex_search(v_command, m, re);
    std::string result = m.str(1);
    return std::stoi(result);
}

bool compareEdge(Edge p1, Edge p2) {
    return (p1.x < p2.x) || (p1.x == p2.x && p1.y <p2.y);
//    return (i < j);
}

std::vector <Edge> getEdges(const std::string &e_command) {
    std::vector <Edge> EdgeList = {};
    std::regex re(R"(^\s*E\s*\{(.*)}\s*$)");
    std::smatch m;
    std::regex_search(e_command, m, re);
    std::string edge_pair_list = m.str(1);
    re = std::regex(R"(<(\d+,\d+)>)");
    std::sregex_iterator edge_begin(edge_pair_list.begin(), edge_pair_list.end(), re);
    std::sregex_iterator edge_end;

    for (std::sregex_iterator k = edge_begin; k != edge_end; k++) {
        std::smatch match = *k;
        std::string match_str = match.str();
        match_str = match_str.replace(0, 1, "")
                .replace(match_str.length() - 1, 1, "");
        int pos = match_str.find_first_of(',');
        int p1 = std::stoi(match_str.substr(0, pos));
        int p2 = std::stoi(match_str.substr(pos + 1));
        int x = std::min(p1, p2);
        int y = std::max(p1, p2);
        Edge e(x, y);
        EdgeList.push_back(e);
//        if (!EdgeList.empty()) {
//            // insert small coordinates into the front
//            if (EdgeList.front().x >= e.x && EdgeList.front().y >= e.y)
//                EdgeList.insert(EdgeList.begin(), e);
//            else
//                EdgeList.push_back(e);
//        } else
//            EdgeList.push_back(e);
    }
//    std::sort (EdgeList.begin(), EdgeList.end(), compareEdge);
//    for (unsigned int i = 0; i < EdgeList.size(); ++i) {
//        std::cout << EdgeList[i].x << "," << EdgeList[i].y << std::endl;
//    }
    return EdgeList;
}

bool find_duplicate_vetex(std::vector<int> list, int p) {
    if (list.empty())
        return false;

    for (unsigned int i = 0; i < list.size(); i++) {
        if (p == list[i]) {
            return true;
        }
    }
    return false;
}

void print_vertext(const std::string &s, std::vector<int> &vertex_cover) {
    unsigned int c = 0;
    std::cout << s << " ";
    if (!vertex_cover.empty()) {
        std::sort(vertex_cover.begin(), vertex_cover.end());
        for (; c < vertex_cover.size() - 1; ++c) {
            std::cout << vertex_cover[c] << ",";
        }
        std::cout << vertex_cover[c] << std::endl;
    } else {
        std::cout << "" << std::endl;
    }
}

std::vector <std::vector<int>> getAdjList(const int &vertics,
                                          const std::vector <Edge> &edge_list) {
    std::vector <std::vector<int >> AdjList(vertics);
    for (unsigned int i = 0; i < edge_list.size(); i++) {
        int p1 = edge_list[i].x;
        int p2 = edge_list[i].y;
        if (p1 > vertics - 1 || p2 > vertics - 1) {
            std::cout << "Error: Vertex exceeds upper boundary." << std::endl;
            adj_list.clear();
            break;
        }
        if (find_duplicate_vetex(AdjList[p2], p1) && find_duplicate_vetex(AdjList[p1], p2)) {
            std::cout << "Error: Edges has already been input before." << std::endl;
            AdjList.clear();
            break;
        } else {
            AdjList[p1].push_back(p2);
            AdjList[p2].push_back(p1);
        }
    }
    for (unsigned int j = 0; j < AdjList.size(); ++j) {
        sort(AdjList[j].begin(), AdjList[j].end());

    }
    return AdjList;

}

void *APPROX_VC_1(void *args) {
    vertex_cover_appro1.clear();
    while (!edge_list_appro1.empty()) {
        int max_index = -1;
        int max_len = -1;
        for (unsigned int i = 0; i < adj_list.size(); ++i) {
            int temp = adj_list[i].size();
            if (max_len < temp) {
                max_index = i;
                max_len = temp;
            }
        }
        // TODO, how to make sure it will eliminate
        std::vector<int>::iterator it;
        for (unsigned int e = 0; e < adj_list[max_index].size(); ++e) {
            // adjacent edge
            int a_e = adj_list[max_index][e];
            it = std::find(adj_list[a_e].begin(), adj_list[a_e].end(), max_index);
            adj_list[a_e].erase(it);
        }
        adj_list[max_index].clear();

        // delete incident edge
        bool deleted = false;
        for (unsigned int j = 0; j < edge_list_appro1.size(); ++j) {
            if (max_index == edge_list_appro1[j].x || max_index == edge_list_appro1[j].y) {
                edge_list_appro1.erase(edge_list_appro1.begin() + j);
                deleted = true;
                --j;
            }
        }
        if (deleted) {
            vertex_cover_appro1.push_back(max_index);
//            adj_list[max_index].clear();
        }

    }

    pthread_getcpuclockid(pthread_self(), &clock_appro1);
    result_appro1 = pclock("APPROX-VC-1 time(us) :", clock_appro1);
    return NULL;
}

void *CNF_SAT_VC(void *args) {
    vertex_cover_sat.clear();
    for (int k = 1; k < n + 1; ++k) {
        std::unique_ptr <Minisat::Solver> solver(new Minisat::Solver());
        Minisat::Lit n_matrix[n][k];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < k; ++j) {
                n_matrix[i][j] = Minisat::mkLit(solver->newVar());
            }
        }
        // 1. At least one vertex is the ith vertex in the vertex cover
        for (int i = 0; i < k; ++i) {
            Minisat::vec <Minisat::Lit> vertex_type_one;
            for (int j = 0; j < n; ++j) {
                vertex_type_one.push(n_matrix[j][i]);
            }
            solver->addClause(vertex_type_one);
        }

        // 2. No one vertex can appear twice in a vertex cover.
        for (int m = 0; m < n; ++m) {
            for (int p = 0; p < k - 1; ++p) {
                for (int q = p + 1; q < k; ++q) {
                    solver->addClause(~n_matrix[m][p], ~n_matrix[m][q]);
                }
            }
        }

        // 3. No more than one vertex appears in the mth position of the vertex cover.
        for (int m = 0; m < k; ++m) {
            for (int p = 0; p < n - 1; ++p) {
                for (int q = p + 1; q < n; ++q) {
                    solver->addClause(~n_matrix[p][m], ~n_matrix[q][m]);
//                solver->addClause(~n_matrix[p][m], ~n_matrix[p + 1][m]);
                }
            }
        }

        // 4. Every edge is incident to at least one vertex in the vertex cover.
        for (unsigned int i = 0; i < edge_list_sat.size(); i++) {
            int p1 = edge_list_sat[i].x;
            int p2 = edge_list_sat[i].y;
            Minisat::vec <Minisat::Lit> vertex_type_four;
            for (int c = 0; c < k; ++c) {
                vertex_type_four.push(n_matrix[p1][c]);
                vertex_type_four.push(n_matrix[p2][c]);
            }
            solver->addClause(vertex_type_four);
        }

        bool res = solver->solve();
//        std::cout << "New result is: " << res << "\n";
        if (res) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < k; ++j) {
//                    std::cout << Minisat::toInt(solver->modelValue(n_matrix[i][j])) << " ";
                    int cover = Minisat::toInt(solver->modelValue(n_matrix[i][j]));
                    if (cover == 0)
                        vertex_cover_sat.push_back(i);
                }
            }
            break;
        }
    }

    pthread_getcpuclockid(pthread_self(), &clock_sat);
    result_sat = pclock("CNF-SAT-VC time(us) :", clock_sat);
    return NULL;
}

void *APPROX_VC_2(void *args) {
    vertex_cover_appro2.clear();

    while (!edge_list_appro2.empty()) {
        Edge first_elem = edge_list_appro2.front();

        int p1 = first_elem.x;
        int p2 = first_elem.y;
        if (!(std::find(vertex_cover_appro2.begin(), vertex_cover_appro2.end(), p1) != vertex_cover_appro2.end()) &&
            !(std::find(vertex_cover_appro2.begin(), vertex_cover_appro2.end(), p2) != vertex_cover_appro2.end())) {
            vertex_cover_appro2.push_back(p1);
            vertex_cover_appro2.push_back(p2);
            edge_list_appro2.erase(edge_list_appro2.begin());
        } else {
            edge_list_appro2.erase(edge_list_appro2.begin());
        }
    }

    pthread_getcpuclockid(pthread_self(), &clock_appro2);
    result_appro2 = pclock("APPROX-VC-2 time(us) :", clock_appro2);
    return NULL;
}


std::string pclock(const std::string &msg, clockid_t cid) {
    struct timespec ts;
    std::stringstream s;
    s.precision(15);
    if (clock_gettime(cid, &ts) == -1)
        s << "clock_gettime failure: ";
    int microseconds = std::round((double) ts.tv_sec * 1000000 + (double) ts.tv_nsec / 1000);
    s << msg << std::to_string(microseconds);
    std::string re = s.str();
    return re;
}


int main(int argc, char **argv) {
    // get line until EOF
    while (std::getline(std::cin, command)) {
        std::regex re(R"(^\s*([VEs])\s.*$)");
        std::smatch m;
        if (std::regex_search(command, m, re)) {
            std::string c = m.str(1);
            if (c == "V") {
                // V handler
                if (isValidVerticesCommandFormat(command)) {
                    vertices = getVertices(command);
                    if (vertices < 1) {
                        std::cerr << "Error: Vertices can't be less than 1." << std::endl;
                        std::cerr.flush();
                    }
                }
            } else if (c == "E") {
                // E handler
                if (isValidEdgeCommandFormat(command)) {
                    if (vertices > 0) {
                        // assign values from from related Edge and Adjlist function.
                        edge_list = getEdges(command);
                        //todo: check if empty??
                        edge_list_appro1 = edge_list;
                        edge_list_appro2 = edge_list;
                        edge_list_sat = edge_list;

                        adj_list = getAdjList(vertices, edge_list);

                        n = vertices;

                        if (!edge_list.empty()) {
                            pthread_create(&thread_sat, NULL, CNF_SAT_VC, NULL);
                            pthread_create(&thread_appro1, NULL, APPROX_VC_1, NULL);
                            pthread_create(&thread_appro2, NULL, APPROX_VC_2, NULL);
                            struct timespec timeout_sat;
                            if (clock_gettime(CLOCK_REALTIME, &timeout_sat) == -1) {
                                /* Handle error */
                            }
                            timeout_sat.tv_sec += 30;
                            pthread_join(thread_appro1, NULL);
                            pthread_join(thread_appro2, NULL);
                            int re = pthread_timedjoin_np(thread_sat, NULL, &timeout_sat);

//                            output vertex cover
                            if (re != 0) {
                                pthread_cancel(thread_sat);
                                std::cout << "CNF-SAT-VC: timeout " << std::endl;

                            } else {
                                print_vertext("CNF-SAT-VC:", vertex_cover_sat);
//                                std::cout << result_sat << std::endl;
                            }
                            print_vertext("APPROX-VC-1:", vertex_cover_appro1);
                            print_vertext("APPROX-VC-2:", vertex_cover_appro2);

//                            output time
//                            if (re != 0) {
//                                pthread_cancel(thread_sat);
//                                std::cout << result_sat << std::endl;
////                                print_vertext("CNF-SAT-VC:", vertex_cover_sat);
//                            }else{
//                                std::cout << "CNF-SAT-VC: timeout " << std::endl;
//                            }
//                            std::cout << result_appro1 << std::endl;
//                            std::cout << result_appro2 << std::endl;

                            continue;
                        } else {
                            std::cout << "CNF-SAT-VC:" << std::endl;
                            std::cout << "APPROX-VC-1:" << std::endl;
                            std::cout << "APPROX-VC-2:" << std::endl;
                        }
                    } else {
                        std::cerr << "Error: Vertices is not defined yet." << std::endl;
                        std::cerr.flush();
                    }
                }
            } else {
                std::cerr << "Error: A5 Project Command not recognized." << std::endl;
            }
        } else {
            std::cerr << "Error: A5 Project Command not recognized." << std::endl;
        }
    }
    exit(EXIT_SUCCESS);
}
