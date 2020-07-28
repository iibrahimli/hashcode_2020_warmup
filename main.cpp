#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <random>

using namespace std;
using ulong = unsigned long;


ulong compute_score(const vector<ulong> &n_slices, const vector<ulong> &order, ulong max_slices){
    ulong slices = 0;

    for(ulong i=0; i<order.size(); ++i){
        slices += n_slices[order[i]];
    }

    if(slices > max_slices)
        return 0;
    else
        return slices;
}


ulong compute_score(const vector<ulong> &n_slices, const vector<bool> &order, ulong max_slices){
    ulong slices = 0;

    for(ulong i=0; i<order.size(); ++i){
        slices += order[i] * n_slices[i];
    }

    if(slices > max_slices)
        return 0;
    else
        return slices;
}


// greedy solution
vector<ulong> greedy_sol(vector<ulong> &n_slices, ulong max_slices){
    vector<ulong> sol;
    vector<bool> used(n_slices.size(), false);
    ulong score = 0;

    int n_loops = 0;
    while(score <= max_slices && n_loops < n_slices.size()){
        ulong temp = compute_score(n_slices, sol, max_slices);
        for(int i=n_slices.size()-1; i>=0; --i){
            if(score + n_slices[i] <= max_slices){
                sol.emplace_back(i);
                score += n_slices[i];
                used[i] = true;
            }
        }
        ++n_loops;
    }

    return sol;
}


// vector<vector<bool>> generate_neighbors(const vector<bool> &bin_sol, const vector<ulong> &n_slices, ulong max_slices){
//     vector<vector<bool>> neighbors;
//     neighbors.reserve(bin_sol.size());
//     vector<bool> candidate;
//     ulong score = compute_score(n_slices, bin_sol, max_slices);

//     for(int i=0; i<bin_sol.size(); ++i){
//         if(!bin_sol[i] && score + n_slices[i] > max_slices)
//             continue;
//         neighbors.push_back(bin_sol);
//         neighbors[neighbors.size()-1][i].flip();
//     }
    
//     return neighbors;
// }


void generate_neighbors(const vector<bool> &bin_sol, vector<vector<bool>> &neighbors, const vector<ulong> &n_slices, ulong max_slices){
    // ulong score = compute_score(n_slices, bin_sol, max_slices);

    // #pragma omp for
    for(int i=0; i<neighbors.size(); ++i){
        for(int j=0; j<bin_sol.size(); ++j){
            neighbors[i][j] = bin_sol[j];
        }
        neighbors[i][i].flip();
    }
    for(int j=0; j<1000; ++j)
        neighbors[rand() % neighbors.size()][rand() % neighbors[0].size()].flip();

}


vector<ulong> sol_bool2ulong(const vector<bool> &bin_sol){
    vector<ulong> sol;
    for(int i=0; i<bin_sol.size(); ++i){
        if(bin_sol[i]){
            sol.emplace_back(i);
        }
    }
    return sol;
}


vector<bool> sol_ulong2bool(const vector<ulong> &sol, ulong n_pizzas){
    vector<bool> bin_sol(n_pizzas, false);
    for(int i=0; i<sol.size(); ++i){
        bin_sol[sol[i]] = true;
    }
    return bin_sol;
}


// simulated annealing
vector<ulong> sa_sol(vector<ulong> &n_slices, ulong max_slices, string ofname){
    vector<bool> bin_sol(n_slices.size(), false);
    ulong score = 0;
    double INITIAL_TEMP;
    double FINAL_TEMP   = 10e-5;
    double TEMP_DECAY   = 0.95;
    int    N_CYCLES     = 250;
    int    MAX_ITER     = 10000;

    vector<vector<bool>> neighbors(n_slices.size());
    for(int i=0; i<n_slices.size(); ++i){
        neighbors[i] = vector<bool>(n_slices.size());
    }
    vector<ulong> neighbor_scores(neighbors.size());
    long delta_score;
    int iter = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0, 1);

    // for(int i=0; i<bin_sol.size(); ++i){
    //     if(dist(gen) > 0.5)
    //         bin_sol[i] = true;
    // }

    // start from the result of the greedy solution
    bin_sol = sol_ulong2bool(greedy_sol(n_slices, max_slices), n_slices.size());

    INITIAL_TEMP = (max_slices - compute_score(n_slices, bin_sol, max_slices)) * 2;

    for(double temp=INITIAL_TEMP; iter<=MAX_ITER && temp>=FINAL_TEMP; temp*=TEMP_DECAY){
        for(ulong cycles=0; cycles<N_CYCLES; ++cycles){

            generate_neighbors(bin_sol, neighbors, n_slices, max_slices);

            // #pragma omp parallel for
            for(int i=0; i<neighbors.size(); ++i){
                neighbor_scores[i] = compute_score(n_slices, neighbors[i], max_slices);
            }
            
            for(int i=0; i<neighbors.size(); ++i){
                delta_score = neighbor_scores[i] - score;
                
                // better solution
                if(delta_score >= 0){
                    bin_sol = neighbors[i];
                    score = neighbor_scores[i];
                }
                else{
                    double r = dist(gen);
                    if(r <= exp(delta_score / temp)){
                        bin_sol = neighbors[i];
                        score = neighbor_scores[i];
                    }
                }
            }
        }

        if(iter % 1 == 0){
            cerr << "[SA] " << fixed << right << setw(6) << iter << "  -  temp: " << temp << "  score: " << score << endl;
            
            // save solution
            auto tsol = sol_bool2ulong(bin_sol);
            ofstream out_stream(ofname);
            if(out_stream){
                out_stream << tsol.size() << endl;
                for(int i=0; i<tsol.size(); ++i)
                    out_stream << tsol[i] << " ";
                out_stream << endl;
                out_stream.close();
            }
        }
        ++iter;
    }

    auto sol = sol_bool2ulong(bin_sol);
    return sol;
}


int main(int argc, char *argv[]){
    
    // parse args
    if(argc != 3){
        cerr << "Invalid number of arguments" << endl;
        return 1;
    }
    string in_filename;
    string out_filename;

    in_filename = argv[1];
    out_filename = argv[2];

    ifstream in_stream(in_filename);

    if(!in_stream){
        cerr << "Unable to open the file " << in_filename << endl;
        return 1;
    }

    ulong M;                        // max number of pizza slices to order
    ulong N;                        // number of different types of pizza

    in_stream >> M >> N;

    vector<ulong> n_slices(N, 0);   // number of slices in each type of pizza
    
    for(int i=0; i<N; ++i)
        in_stream >> n_slices[i];
    
    in_stream.close();

    // cout << "M = " << M << ", N = " << N << endl;
    // for(auto nsl : n_slices)
    //     cout << nsl << " ";
    // cout << endl;


    // ======================== compute ========================

    vector<ulong> order;            // which types of pizza to order

    // order = greedy_sol(n_slices, M);

    order = sa_sol(n_slices, M, out_filename);

    cout << "score: " << compute_score(n_slices, order, M) << endl;


    // ========================  save   ========================
    
    ofstream out_stream(out_filename);

    if(!out_stream){
        cerr << "Unable to open the file " << out_filename << endl;
        return 1;
    }

    out_stream << order.size() << endl;
    for(int i=0; i<order.size(); ++i)
        out_stream << order[i] << " ";
    out_stream << endl;

    out_stream.close();

    cout << "saved " << out_filename << endl;

    return 0;
}