#include "algorithm.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <queue>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdint>

using namespace std;

inline rank_t update_rank(pair_info const& info, rank_t current, unsigned char next, unsigned char pop) {
    current -= info.rank_values[pop] * info.last_multiplier;
    rank_t result = current * info.rank_base + info.rank_values[next];
    return result;
}

inline rank_t update_rank_hash(pair_info const& info, rank_t current, unsigned char next, unsigned char pop) {
    rank_mul_t ncurrent = current + RABIN_MODULO - info.rank_values[pop] * info.last_multiplier;
    rank_t result = ((rank_mul_t)ncurrent * (rank_mul_t)info.rank_base + (rank_mul_t)info.rank_values[next])
            % RABIN_MODULO;
    return result;
}

void rank_init(pair_info &info, int kvalue) {

    if (kvalue <= 0) {
        printf("K value must be greater than 0.");
        exit(1);
    }

    info.kvalue = kvalue;
    int rank = 0;
    for (int i = 0; i < 256; i++) {
        if (info.rank_counters[i] > 0) info.rank_values[i] = rank++;
    }
    info.rank_base = rank;
    info.last_multiplier = 1;

    bool has_overflow = false;

    int tmp_kvalue = kvalue - 1;
    while (tmp_kvalue--) {
        rank_t ovflw_test = info.last_multiplier;
        info.last_multiplier *= info.rank_base;
        if (has_overflow) {
            info.last_multiplier %= RABIN_MODULO;
        }
        else if ((ovflw_test > info.last_multiplier) ||
                (ovflw_test * info.rank_base > info.last_multiplier * info.rank_base)) {
            has_overflow = true;
            info.hash_fallback = true;
            info.last_multiplier = ((ovflw_test % RABIN_MODULO) * info.rank_base) % RABIN_MODULO;
            cout << "Hashing fallback!" << endl;
        }
    }

    if (!has_overflow) {
        auto rank_tmp = info.last_multiplier * info.rank_base;
        info.rank_byte_order = 0;
        while (rank_tmp) {
            info.rank_byte_order++;
            rank_tmp /= 256;
        }
    }
    else {
        info.rank_byte_order = 8;
    }
}

template <class RF>
vector<rank_t> do_ranking(pair_info const& info, const char *chars, size_t len, RF ranking_function) {

    vector<rank_t> result(len - info.kvalue + 1);
    rank_t rank = 0;

    for (size_t i = 0; i < info.kvalue; i++) {
        rank = ranking_function(info, rank, chars[i], 0);
    }
    result[0] = rank;

    for (size_t i = info.kvalue; i < len; i++) {
        rank = ranking_function(info, rank, chars[i], chars[i - info.kvalue]);
        result[i - info.kvalue + 1] = rank;
    }
    return result;
}

void counting_sort(vector<rank_t> &vec, unsigned int byte_ref) {
    int counts[256 + 1] = {0};

    vector<rank_t> tmp = vec;

    for (auto val : vec) counts[((val >> (byte_ref * 8u)) & 0xFFu) + 1]++;

    for (int i = 1; i < 256; i++) {
        counts[i] += counts[i-1];
    }

    for (auto val : tmp) {
        vec[counts[(val >> (byte_ref * 8u)) & 0xFFu]++] = val;
    }
}

inline unsigned char ext_byte(rank_t value, int byteidx) {
    return (value >> (byteidx * 8u)) & 0xFFu;
}

template <class T, class F>
void counting_sort_ext(vector<T> &vec, unsigned int byte_ref, F mapping) {
    int counts[256 + 1] = {0};

    vector<T> tmp = vec;

    for (auto val : vec) counts[ext_byte(mapping(val), byte_ref) + 1]++;

    for (int i = 1; i < 256; i++) {
        counts[i] += counts[i-1];
    }

    for (auto val : tmp) {
        vec[counts[ext_byte(mapping(val), byte_ref)]++] = val;
    }
}

pair_info preprocess_sequences(const pangene_idata &data, uint32_t kvalue, bool onlyComplexity) {

    pair_info info = pair_info {};

    seq_id_t seq_count = data.sequences.size();
    info.sequences_count = seq_count;
    info.seq_gen_mapping.resize(seq_count);
    info.kseq_lengths.resize(seq_count);

    // Seq number is greater than genomes count
    info.genome_sequences.resize(seq_count);
    unsigned int genomes_count = 0;

    memset(info.rank_counters, 0, 256 * sizeof(info.rank_counters[0]));
    for (seq_id_t i = 0; i < seq_count; i++) {
        for (size_t j = 0; j < data.sequences[i].size(); j++) {
            info.rank_counters[(unsigned char)data.sequences[i][j]]++;
        }
    }

    rank_init(info, kvalue);

    auto &kmers = info.kmers;

    for (seq_id_t i = 0; i < seq_count; i++) {

        unsigned int gen_id = data.seq_to_genomes[i];

        info.seq_gen_mapping[i] = gen_id;
        genomes_count = max(genomes_count, gen_id + 1);

        info.genome_sequences[gen_id].push_back(i);

        size_t len = data.sequences[i].size();

        int64_t kseq_len = (int64_t)len - info.kvalue + 1;

        if (kseq_len > 0) {
            info.kseq_lengths[i] = kseq_len;

            for (auto kmer : do_ranking(info, data.sequences[i].data(), len,
                                        info.hash_fallback ? update_rank_hash : update_rank)) {
                kmers.emplace_back(kmer, i, 1);
            }
        }
        else {
            info.kseq_lengths[i] = 0;
        }
    }

    info.genomes_count = genomes_count;
    info.genome_sequences.resize(genomes_count);

    seq_id_t max_sortval = 1;
    for (int b = 0; max_sortval <= seq_count; b++) {
        counting_sort_ext(kmers, b, [](auto k) { return k.seq; });
        max_sortval *= 256;
    }

    for (int b = 0; b < info.rank_byte_order; b++) {
        counting_sort_ext(kmers, b, [](auto k) { return k.rank; });
    }

    for (int i = 1; i < kmers.size(); i++) {
        if (kmers[i-1].rank > kmers[i].rank) {
            printf("Not sorted ranks at pos (%d %d) vals: (%d %d)\n", i-1, i, kmers[i-1].rank, kmers[i].rank);
            exit(1);
        }
    }

    int uniqi = 0;
    for (int i = 1; i < kmers.size(); i++) {
        if (!kmers[uniqi].dupl(kmers[i])) {
            kmers[++uniqi] = kmers[i];
            kmers[uniqi].count = 1;
        }
    }
    kmers.resize(uniqi+1);

    info.kmers_ranges_start.resize(info.sequences_count);
    for (int i = 0; i < info.sequences_count; i++) {
        info.kmers_ranges_start.reserve(info.kseq_lengths[i]);
    }

    info.computation_costs.resize(info.sequences_count);
    vector<vector<seq_id_t>> &kmer_ranges = info.kmers_ranges_start;

    rank_t current_rank = kmers.begin()->rank;
    int current_rank_start = 0;
    long ranges_count = 0;
    for (int i = 0; i < kmers.size(); i++) {
        bool last = (i == kmers.size() - 1);
        if (current_rank != kmers[i].rank || last) {

            // Process current rank
            int prev_rank_begin = current_rank_start;
            int prev_rank_end = i + (last ? 1 : 0);

            if (prev_rank_end - prev_rank_begin > 1) {

                if (!onlyComplexity) {
                    // Sort kmers for memory cache optimizations
                    sort(kmers.begin() + prev_rank_begin, kmers.begin() + prev_rank_end,
                         [](auto const &a, auto const &b) {
                             return a.seq < b.seq;
                         });
                }

                for (int j = prev_rank_begin; j < prev_rank_end; j++) {
                    if (!onlyComplexity) {
                        kmer_ranges[kmers[j].seq].push_back(prev_rank_begin);
                    }
                    info.computation_costs[kmers[j].seq].total_visited += prev_rank_end - prev_rank_begin;
                }
                ranges_count++;
            }

            current_rank_start = i;
            current_rank = kmers[i].rank;
        }
    }

    long global_cost = 0;
    long total_sequences_length = 0;
    for (int i = 0; i < info.computation_costs.size(); i++) {
        global_cost += info.computation_costs[i].total_visited;
        total_sequences_length += info.kseq_lengths[i];
        info.computation_costs[i].total_seq_length = info.kseq_lengths[i];
        info.computation_costs[i].linear_ratio = ((float)info.computation_costs[i].total_visited / info.kseq_lengths[i]);
//        cout << info.computation_costs[i].first.first << " => " << info.computation_costs[i].second << endl;
    }

    cout << "------------" << endl;
    cout << "COMPUTATIONAL COSTS: " << endl;
    cout << "Total cost: " << global_cost << " lookups" << endl;
    cout << "Linear ratio: " << ((float)global_cost / (float)total_sequences_length) << endl;

    float estimated_ops_per_ms = 40505.500586716735;

    long ms = (long)((float)global_cost / estimated_ops_per_ms);

    cout << "Estimated time: ";
    if (ms < 10000) {
        cout << ms << "ms" << endl;
    }
    else if (ms < 1000 * 180) {
        cout << (ms / 1000) << "s" << endl;
    }
    else if (ms < 1000 * 60 * 180) {
        cout << (ms / 1000 / 60) << "min" << endl;
    }
    else {
        cout << (ms / 1000 / 60 / 60) << "hours" << endl;
    }

    cout << "------------" << endl << endl;

    return info;
}

static inline unsigned int get_genome_index(pair_info const& info, seq_id_t id) {
    return info.seq_gen_mapping[id];
}


static unsigned int compute_scores_index(pair_info const& info, unsigned int row, unsigned int col) {
    return row * info.genomes_count + info.seq_gen_mapping[col];
}

scores compute_scores(const pair_info &info, const vector<seq_id_t> &sequences, uint32_t step_size) {

    // Non-zero cells list
    vector<mat_cell> nonzero_cells;

    unsigned int row_seqs_count = sequences.size();

    /* Global max intra and inter scores */
    vector<float> max_scores = vector<float>(row_seqs_count * info.genomes_count, 0);
    vector<float> col_max_scores = vector<float>(info.sequences_count, 0);

    /* Single row data, overwritten on each row cycle  */
    vector<int> row_intersection_size = vector<int>(info.sequences_count, 0);
    vector<int> row_connection_perc_cnt = vector<int>(info.sequences_count, 0);
    vector<int> row_transposed_perc_cnt = vector<int>(info.sequences_count, 0);
    vector<unsigned int> reset_colors = vector<unsigned int>(info.sequences_count, 0);
    vector<int> colored_cells = vector<int>();
    colored_cells.reserve(info.sequences_count);

    vector<seq_id_t> flat_map = vector<seq_id_t >(info.sequences_count, INT32_MAX);
    int flat_idx = 0;
    for (seq_id_t row : sequences) {
        flat_map[row] = flat_idx++;
    }

    // Color cells to avoid zero-initialization
    unsigned int color = 0;

    for (seq_id_t row : sequences) {

        // Update color
        color++;
        colored_cells.clear();

//        /* Reset all row data [this step was quadratic, now optimized with colors]*/
//        memset(row_intersection_size.data(), 0, row_intersection_size.size() * sizeof(int));
//        memset(row_connection_perc_cnt.data(), 0, row_connection_perc_cnt.size() * sizeof(int));
//        memset(row_transposed_perc_cnt.data(), 0, row_transposed_perc_cnt.size() * sizeof(int));

        const kmer_rank *kmers_endptr = info.kmers.data() + info.kmers.size();

        int index = 0;
        for (auto range_start : info.kmers_ranges_start[row]) {

            auto last_it = &info.kmers[range_start];
            auto it = &info.kmers[range_start];

            count_t my_cnt;
            while (it < kmers_endptr && last_it->rank == it->rank) {
                if (it->seq == row) {
                    my_cnt = it->count;
                    break;
                }
                it++;
            }

            it = &info.kmers[range_start];
            while (it < kmers_endptr && last_it->rank == it->rank) {

                // Reinitialize the cell if color is different
                if (reset_colors[it->seq] != color) {
                    reset_colors[it->seq] = color;
                    colored_cells.push_back(it->seq);
                    row_intersection_size[it->seq] = 0;
                    row_connection_perc_cnt[it->seq] = 0;
                    row_transposed_perc_cnt[it->seq] = 0;
                }

                row_intersection_size[it->seq] += min(it->count, my_cnt);
                row_connection_perc_cnt[it->seq] += my_cnt;
                row_transposed_perc_cnt[it->seq] += it->count;
                it++;
            }
            index++;
        }

        // Clear the identity cell, it's faster than an 'if' inside the loop
        row_intersection_size[row] = 0;
        row_connection_perc_cnt[row] = 0;
        row_transposed_perc_cnt[row] = 0;

        /**
         * Faster way to compute maximums in-place,
         * with the formula |union(A, B)| = |A| + |B| - |intersection(A, B)|
         */
        for (unsigned int col_index : colored_cells) {
            int my_kcnt = info.kseq_lengths[row];
            int other_kcnt = info.kseq_lengths[col_index];
            int union_size = my_kcnt + other_kcnt - row_intersection_size[col_index];
            float perc = (float) row_connection_perc_cnt[col_index] / (float) my_kcnt;
            float tr_perc = (float) row_transposed_perc_cnt[col_index] / (float) other_kcnt;
            float threshold = 1.0f / (2.0f * (float) info.kvalue);
            bool score_valid = perc >= threshold || tr_perc >= threshold;
            float score = (float) row_intersection_size[col_index] / (float) union_size *
                          (score_valid ? 1.0f : 0.0f);

            // Add to nonzero cells only if score is > 0
            if (score > 0.0f) {
                nonzero_cells.push_back(mat_cell{
                        .score = score,
                        .perc = perc,
                        .tr_perc = tr_perc,
                        .row = row,
                        .col = col_index
                });
                float &current_max_score = max_scores[compute_scores_index(info, flat_map[row], col_index)];
                current_max_score = max(current_max_score, score);
                col_max_scores[col_index] = max(col_max_scores[col_index], score);
            }
        }
    }

    return move(scores {
            .sequences = sequences,
            .flat_map = move(flat_map),
            .non_zero = move(nonzero_cells),
            .max_scores = move(max_scores),
            .col_max_scores = move(col_max_scores)
    });
}