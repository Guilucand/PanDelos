#include <cstdio>
#include <string>
#include <pangene_idata.h>
#include <algorithm.h>
#include <mutex>
#include <thread>
#include <atomic>
#include <stats.h>
#include <fstream>
#include "cxxopts.hpp"

using namespace std;

void save_network_to_file(ofstream &output, map<std::pair<seq_id_t, seq_id_t>, float> const& network) {
    for (auto &element : network) {
        if (element.first.first < element.first.second) {
            output << element.first.first << "\t" << element.first.second << "\t" << element.second << endl;
        }
    }
}

int main(int argc, char **argv) {

    cxxopts::Options options(argv[0], " - PanDelos");

    bool compute_complexity = false;

    options
            .add_options()
                    ("i,input", "Input  file (.faa) to process", cxxopts::value<string>())
                    ("k,kvalue", "Length of the kmers used by the algorithm",
                     cxxopts::value<uint32_t>())
                    ("c,complexity",
                     "Compute the required number of operations without computing the network (fast)",
                     cxxopts::value<bool>(compute_complexity))
                    ("o,output", "Output file for the network", cxxopts::value<std::string>())
                    ("j,threads",
                     "Number of threads to use for the computation, defaults to # of processors",
                     cxxopts::value<int>()->default_value("8"))
                    ("h,help", "Print this help message");

    auto arguments = options.parse(argc, argv);

    if (arguments.count("help")) {
        cout << options.help({"", "Group"}) << endl;
        exit(0);
    }

    if (!arguments.count("input")) {
        cout << "Input .faa file is required! (--input flag)" << endl;
        exit(0);
    }

    uint32_t kvalue;
    if (!arguments.count("kvalue")) {
        cout << "K value is required! (--kvalue flag)" << endl;
        exit(0);
    } else {
        kvalue = arguments["kvalue"].as<uint32_t>();
    }

    pangene_idata data = pangene_idata(arguments["input"].as<string>());

    if (compute_complexity) {
        preprocess_sequences(data, kvalue, true);
        exit(0);
    }

    if (!arguments.count("output")) {
        cout << "Output network is required! (--output flag)" << endl;
        exit(0);
    }

    ofstream output(arguments["output"].as<string>());
    if (!output.is_open()) {
        cout << "Cannot open output file for writing!" << endl;
        exit(0);
    }

    size_t threads_count = 0;
    if (arguments.count("threads")) {
        threads_count = arguments["threads"].as<int>();
    }
    if (threads_count <= 0) {
        threads_count = 8;
    }


    auto start_time = chrono::steady_clock::now();
    auto algo_info = preprocess_sequences(data, kvalue, false);

    /*
     * Total complexity: |gset|**2 * ((|seq(gi)| + |seq(gj)|)**2 + ...)
     * Worst case complexity: |gset|**2 * ((|seq(gi)| + |seq(gj)|)**2 + len(seq(gi)) * (|seq(gi)| + |seq(gj)|))
     * New complexity: |gset|**2 * (|seq(gi)| * |seq(gj)| * len(avg(seq(gi), seq(g2))))
     * New complexity2: (|seq| * |seq|) + |seq| * avg(kmers_repeats)))
     *
     * */

    size_t genomes_count = data.genomes.size();
    bool multithread = threads_count > 1;

    mutex cout_mutex;
    mutex result_mutex;

    map<pair<seq_id_t, seq_id_t>, float> result_network;
    size_t done_genomes = 0;

    auto compute_genome = [&](size_t g) {

        vector<pair<pair<seq_id_t, seq_id_t>, float>> connections;

        cout << "Working on genome " << g << "/" << genomes_count << endl;
        scores scores_part = compute_scores(
                algo_info,
                algo_info.genome_sequences[g],
                multithread ? 2048 : 0xFFFFFFFF);

        cout << "Preprocessed genome " << g << "/" << genomes_count << endl;
        cout << "Filtered count: " << scores_part.non_zero.size() << endl;

        vector<float> inter_thr_sum = vector<float>(genomes_count, 0.0f);
        vector<float> inter_thr_count = vector<float>(genomes_count, 0.0f);
        vector<float> inter_max_score = vector<float>(genomes_count, 0.0f);
        vector<float> inter_min_score = vector<float>(genomes_count, 1.0f);

        vector<float> inter_max_perc = vector<float>(genomes_count);
        vector<float> inter_min_perc = vector<float>(genomes_count, 1.0f);

        // The engagged array should be different for each comparison genome
        // leading to an array of [sequencesCount][genomes_count] size.
        // this should then be used to determine if we should consider
        // the inter_min_score of the i-th genome to add this sequence to the map
        // An alternative way is to iterate two times the scores, the first time updating only
        // the inter_min/max_score arrays, then updating the relative threshold only when we set engagged to true
        vector<vector<bool>> engagged = vector<vector<bool>>(
                data.sequences.size(),
                vector<bool>(genomes_count));

//            float min_inter_max_score = 1.0f;

        vector<bool> should_add_connection = vector<bool>(scores_part.non_zero.size());

        for (size_t i = 0; i < scores_part.non_zero.size(); i++) {

            mat_cell const &non_zero = scores_part.non_zero[i];

            if (data.seq_to_genomes[non_zero.row] != data.seq_to_genomes[non_zero.col]) {
                if (non_zero.score == scores_part.max_scores[
                        scores_part.flat_map[non_zero.row] * genomes_count +
                        data.seq_to_genomes[non_zero.col]
                ] &&
                    non_zero.score == scores_part.col_max_scores[non_zero.col]) {

                    connections.push_back({{non_zero.row, non_zero.col}, non_zero.score});
                    connections.push_back({{non_zero.col, non_zero.row}, non_zero.score});

                    should_add_connection[i] = true;

                    int sg = data.seq_to_genomes[non_zero.col];
                    float score = non_zero.score;
                    float perc = non_zero.perc;
                    float otherPerc = non_zero.tr_perc;

                    engagged[non_zero.row][sg] = true;
                    inter_thr_sum[sg] += 2 * score;
                    inter_thr_count[sg] += 2;
                    if (score < 1.0 && score > inter_max_score[sg]) {
                        inter_max_score[sg] = score;
                    }

                    if (score > 0.0 && score < inter_min_score[sg]) {
                        inter_min_score[sg] = score;
                    }

                    inter_max_perc[sg] = max(inter_max_perc[sg], max(perc, otherPerc));
                    inter_min_perc[sg] = min(inter_min_perc[sg], min(perc, otherPerc));
                }
            }
        }

        for (int i = 0; i < genomes_count; i++) {

            if (i == g) continue;

            float inter_thr = inter_thr_sum[i] / inter_thr_count[i];

            cout_mutex.lock();
            {
                cout << "Comparing genome " << g << " with " << i << ":" << endl;

                cout << "Score\t" << inter_thr << "\t" << inter_min_score[i] << "\t"
                     << inter_max_score[i] << endl;
                cout << "Perc\t" << inter_min_perc[i] << "\t" << inter_max_perc[i] << endl;
            }
            cout_mutex.unlock();
        }


        vector<float> scores_row_threshold = vector<float>(data.sequences.size(), +1.0f / 0.0f);

        for (int i = 0; i < scores_part.non_zero.size(); i++) {
            if (should_add_connection[i]) {
                int row = scores_part.non_zero[i].row;
                int sg = data.seq_to_genomes[scores_part.non_zero[i].col];
                scores_row_threshold[row] = min(scores_row_threshold[row], inter_max_score[sg]);
            }
        }

//        for (int i = 0; i < genomes_count; i++) {
//            if (i == finalG) continue; // Exclude current genome from computation
//            min_inter_max_score = Math.min(min_inter_max_score, inter_max_score[i]);
//        }

        /* get inter bbh */
        /* also, calcolate threshold for inter non-bbh as the average non-null and non-bbh scores*/
        for (int i = 0; i < scores_part.non_zero.size(); i++) {

            mat_cell const &non_zero = scores_part.non_zero[i];

            size_t first_genome = data.seq_to_genomes[non_zero.row];
            size_t second_genome = data.seq_to_genomes[non_zero.col];

            if ((non_zero.row < non_zero.col) &&
//            engagged[scores_part.row[i]].get(scores_part.second_seq_genome[i]) &&
                first_genome == second_genome &&
                (non_zero.score ==
                 scores_part.max_scores[scores_part.flat_map[non_zero.row] * genomes_count +
                                        second_genome] &&
                 non_zero.score ==
                 scores_part.max_scores[scores_part.flat_map[non_zero.col] * genomes_count +
                                        second_genome] &&
                 non_zero.score >= scores_row_threshold[non_zero.row] //&&
//									non_zero.score >= min_inter_max_score
            )) {

                connections.push_back({{non_zero.row, non_zero.col}, non_zero.score});
            }
        }

        result_mutex.lock();
        {
            result_network.insert(connections.begin(), connections.end());
        }
        result_mutex.unlock();

        cout_mutex.lock();
        {
            done_genomes += 1;
            cout << "Done genome " << done_genomes << "/" << genomes_count << " => "
                 << scores_part.non_zero.size() << endl;

        }
        cout_mutex.unlock();
    };

    vector<thread> threads;

    atomic_size_t next_genome(0);

    for (int ti = 0; ti < 8; ti++) {
        threads.emplace_back([&]() {
            while (true) {
                int g = next_genome.fetch_add(1);
                if (g >= genomes_count) return;
                compute_genome(g);
            }
        });
    }

    for (thread &t : threads) {
        t.join();
    }

    cout << "Result size: " << result_network.size() << endl;

    auto end_time = chrono::steady_clock::now();

    cout << "Total time: "
         << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << "ms"
         << endl;

    cout << "----------" << endl;
    cout << "undirected degree distribution" << endl;
    print_undirected_degree_distribution(result_network);

    cout << "----------" << endl;
    cout << "directed degree distribution" << endl;
    print_directed_degree_distribution(result_network);


    cout << "----------" << endl;
    cout << "CoCo sizes" << endl;
    print_undirected_connected_components(result_network, data.sequences.size());

    cout << "----------" << endl;
    cout << "writing into " + arguments["output"].as<string>() + "" << endl;
    save_network_to_file(output, result_network);
}
