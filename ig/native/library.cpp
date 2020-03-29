#include "pangene_native.h"

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

typedef uint64_t rank_t;
typedef uint32_t seq_id_t;
typedef uint32_t gen_id_t;

#define RABIN_MODULO ((uint64_t)18446744073709551557ULL)

struct kmer_rank {
    rank_t rank;
    seq_id_t seq;
    unsigned int count;

    kmer_rank() {}

    kmer_rank(rank_t rank, seq_id_t seq, unsigned int count) {
        this->rank = rank;
        this->seq = seq;
        this->count = count;
    }

    bool dupl(kmer_rank &other) {
        bool dupl = this->seq == other.seq && this->rank == other.rank;
        count += dupl;
        return dupl;
    }
};

struct kmers_range {
    kmer_rank *start;
    kmer_rank *current;
    kmer_rank *end;
};

struct seq_computation_cost {

    seq_computation_cost() : total_visited(0), total_seq_length(0), linear_ratio(1.0) {}

    unsigned long total_visited;
    unsigned long total_seq_length;
    float linear_ratio;
};

struct pair_info {
    unsigned long rank_counters[256] = {0};
    unsigned char rank_values[256] = {0};
    unsigned char rank_base = 0;
    unsigned char rank_byte_order = 0;
    unsigned int kvalue = 0;
    bool hash_fallback = false;
    rank_t last_multiplier = 1;

    unsigned int sequences_count = 0;
    unsigned int genomes_count = 0;
    vector<kmer_rank> kmers;
    vector<unsigned int> kseq_lengths;
    vector<gen_id_t> seq_gen_mapping;
    vector<vector<seq_id_t>> genome_sequences;
    vector<vector<kmers_range>> kmers_ranges;
    vector<seq_computation_cost> computation_costs;
} global_info;

inline rank_t update_rank(pair_info const& info, rank_t current, jchar next, jchar pop) {
    current -= info.rank_values[pop] * info.last_multiplier;
    rank_t result = current * info.rank_base + info.rank_values[next];
    return result;
}

inline rank_t update_rank_hash(pair_info const& info, rank_t current, jchar next, jchar pop) {
    __uint128_t ncurrent = current + RABIN_MODULO - info.rank_values[pop] * info.last_multiplier;
    rank_t result = ((__uint128_t)ncurrent * (__uint128_t)info.rank_base + (__uint128_t)info.rank_values[next])
            % RABIN_MODULO;
    return result;
}

void rank_init(pair_info &info, int kvalue) {
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
        else if (ovflw_test > info.last_multiplier) {
            has_overflow = true;
            info.hash_fallback = true;
            info.last_multiplier = ((ovflw_test % RABIN_MODULO) * info.rank_base) % RABIN_MODULO;
            cout << "Hashing fallback!" << endl;
        }
    }

    if (!has_overflow) {
        auto rank_tmp = info.last_multiplier;
        info.rank_byte_order = 0;
        while (rank_tmp) {
            info.rank_byte_order++;
            rank_tmp /= 256;
        }
    }
    else {
        info.rank_byte_order = 7;
    }
}

template <class RF>
vector<rank_t> do_ranking(pair_info const& info, const jchar *chars, jsize len, RF ranking_function) {

    vector<rank_t> result(len - info.kvalue + 1);
    rank_t rank = 0;

    for (jsize i = 0; i < info.kvalue; i++) {
        rank = ranking_function(info, rank, chars[i], 0);
    }
    result[0] = rank;

    for (jsize i = info.kvalue; i < len; i++) {
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

template <class T, class F>
void counting_sort_ext(vector<T> &vec, unsigned int byte_ref, F mapping) {
    int counts[256 + 1] = {0};

    vector<T> tmp = vec;

    for (auto val : vec) counts[((mapping(val) >> (byte_ref * 8u)) & 0xFFu) + 1]++;

    for (int i = 1; i < 256; i++) {
        counts[i] += counts[i-1];
    }

    for (auto val : tmp) {
        vec[counts[(mapping(val) >> (byte_ref * 8u)) & 0xFFu]++] = val;
    }
}

void Java_infoasys_cli_pangenes_PangeneNative_preprocessSequences(JNIEnv *env, jobject obj, jobject data, jint kvalue, jboolean onlyComplexity) {

    // Reset infos
    global_info = pair_info {};

    pair_info &info = global_info;

    jclass dataClass = env->GetObjectClass(data);
    jfieldID seqId = env->GetFieldID(dataClass, "sequences", "Ljava/util/Vector;");
    jfieldID seqGen = env->GetFieldID(dataClass, "sequenceGenome", "Ljava/util/Vector;");

    jobject objSeq = env->GetObjectField(data, seqId);
    jobject objGen = env->GetObjectField(data, seqGen);

    jclass vecClass = env->GetObjectClass(objSeq);
    jmethodID elMethod = env->GetMethodID(vecClass, "get", "(I)Ljava/lang/Object;");
    jmethodID szMethod = env->GetMethodID(vecClass, "size", "()I");

    seq_id_t seq_count = env->CallIntMethod(objSeq, szMethod);
    info.sequences_count = seq_count;
    info.seq_gen_mapping.resize(seq_count);
    info.kseq_lengths.resize(seq_count);

    // Seq number is greater than genomes count
    info.genome_sequences.resize(seq_count);
    unsigned int genomes_count = 0;

    memset(info.rank_counters, 0, 256 * sizeof(info.rank_counters[0]));
    for (seq_id_t i = 0; i < seq_count; i++) {
        auto str = (jstring) env->CallObjectMethod(objSeq, elMethod, i);
        jboolean iscopy = false;
        jsize len = env->GetStringLength(str);
        const jchar *seq = env->GetStringChars(str, &iscopy);
        for (int j = 0; j < len; j++) {
            if (seq[j] < 256) {
                info.rank_counters[seq[j]]++;
            }
        }
        env->ReleaseStringChars(str, seq);
    }

    rank_init(info, kvalue);

    auto &kmers = info.kmers;

    for (seq_id_t i = 0; i < seq_count; i++) {
        auto str = (jstring) env->CallObjectMethod(objSeq, elMethod, i);
        auto gen_id_obj = (jobject)env->CallObjectMethod(objGen, elMethod, i);
        jclass intCls = env->GetObjectClass(gen_id_obj);
        jmethodID getiMethod = env->GetMethodID(intCls, "intValue", "()I");
        gen_id_t gen_id = env->CallIntMethod(gen_id_obj, getiMethod);

        info.seq_gen_mapping[i] = gen_id;
        genomes_count = max(genomes_count, gen_id + 1);

        info.genome_sequences[gen_id].push_back(i);

        jboolean iscopy = false;
        jsize len = env->GetStringLength(str);
        const jchar *seq = env->GetStringChars(str, &iscopy);

        int64_t kseq_len = (int64_t)len - info.kvalue + 1;

        if (kseq_len > 0) {
            info.kseq_lengths[i] = kseq_len;

            for (auto kmer : do_ranking(info, seq, len,
                                        info.hash_fallback ? update_rank_hash : update_rank)) {
                kmers.emplace_back(kmer, i, 0);
            }
        }
        else {
            info.kseq_lengths[i] = 0;
        }

        env->ReleaseStringChars(str, seq);
    }

    info.genomes_count = genomes_count;
    info.genome_sequences.resize(genomes_count);

    seq_id_t max_sortval = 1;
    for (int b = 0; max_sortval >= seq_count; b++) {
        counting_sort_ext(kmers, b, [](auto k) { return k.seq; });
        max_sortval *= 256;
    }

    for (int b = 0; b < info.rank_byte_order; b++) {
        counting_sort_ext(kmers, b, [](auto k) { return k.rank; });
    }

    int uniqi = 0;
    for (int i = 1; i < kmers.size(); i++) {
        if (!kmers[uniqi].dupl(kmers[i])) {
            kmers[++uniqi] = kmers[i];
            kmers[uniqi].count = 1;
        }
    }
    kmers.resize(uniqi+1);

    info.kmers_ranges.resize(info.sequences_count);
    for (int i = 0; i < info.sequences_count; i++) {
        info.kmers_ranges.reserve(info.kseq_lengths[i]);
    }

    info.computation_costs.resize(info.sequences_count);
    vector<vector<kmers_range>> &kmer_ranges = info.kmers_ranges;

    rank_t current_rank = kmers.begin()->rank;
    int current_rank_start = 0;
    long ranges_count = 0;
    for (int i = 0; i < kmers.size(); i++) {
        if (current_rank != kmers[i].rank || (i == kmers.size() - 1)) {

            // Process current rank
            int prev_rank_begin = current_rank_start;
            int prev_rank_end = i;

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
                        kmer_ranges[kmers[j].seq].push_back(
                                kmers_range{
                                        .start = kmers.data() + prev_rank_begin,
                                        .current = kmers.data() + j,
                                        .end = kmers.data() + prev_rank_end
                                });
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
    cout << "Linear ratio: " << ((float)global_cost / total_sequences_length) << endl;

    float estimated_ops_per_ms = 40505.500586716735;

    long ms = (long)(global_cost / estimated_ops_per_ms);

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
}

static inline unsigned int get_genome_index(pair_info const& info, seq_id_t id) {
    return info.seq_gen_mapping[id];
}

struct mat_cell {
    float score;
    float perc;
    float tr_perc;
    unsigned int x;
    unsigned int y;
};

jintArray create_jni_int_array(JNIEnv *env, jint *ptr, jsize len) {
    jintArray arr = env->NewIntArray(len);
    env->SetIntArrayRegion(arr, 0, len, ptr);
    return arr;
}

jfloatArray create_jni_float_array(JNIEnv *env, jfloat *ptr, jsize len) {
    jfloatArray arr = env->NewFloatArray(len);
    env->SetFloatArrayRegion(arr, 0, len, ptr);
    return arr;
}

struct scores {
    vector<seq_id_t> sequences;
    vector<seq_id_t> flat_map;
    vector<mat_cell> non_zero;
    vector<float> max_scores;
    vector<float> col_max_scores;
};

static unsigned int compute_scores_index(pair_info const& info, unsigned int row, unsigned int col) {
    return row * info.genomes_count + info.seq_gen_mapping[col];
}

static scores computeScores(pair_info const& info, vector<seq_id_t> const& sequences, jint step_size) {

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

        vector<kmer_rank*> pointers(info.kmers_ranges[row].size());
        for (size_t i = 0; i < info.kmers_ranges[row].size(); i++) {
            pointers[i] = info.kmers_ranges[row][i].start;
        }

        // Sequences steps to keep memory locality and avoid cache misses
        int sequences_step = 2048;

        for (int max_allowed_sequence = sequences_step;
            max_allowed_sequence < info.sequences_count + sequences_step;
            max_allowed_sequence += sequences_step) {

            int index = 0;
            for (auto &range : info.kmers_ranges[row]) {
                auto &it = pointers[index];
                unsigned int my_cnt = range.current->count;
                while (it < range.end && it->seq <= max_allowed_sequence) {

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
                        .x = row,
                        .y = col_index
                });
                jfloat &current_max_score = max_scores[compute_scores_index(info, flat_map[row], col_index)];
                current_max_score = max(current_max_score, score);
                col_max_scores[col_index] = max(col_max_scores[col_index], score);
            }
        }
    }

    return move(scores {
        .sequences = move(sequences),
        .flat_map = move(flat_map),
        .non_zero = move(nonzero_cells),
        .max_scores = move(max_scores),
        .col_max_scores = move(col_max_scores)
    });
}

void Java_infoasys_cli_pangenes_PangeneNative_computeScores(JNIEnv *env, jobject obj, jint genome_id, jobject out_scores, jint step_size) {

    pair_info const& info = global_info;

    unsigned int rows_count = info.genome_sequences[genome_id].size();

    cout << "Genome " << genome_id << " cost = " << std::accumulate(
            info.genome_sequences[genome_id].begin(),
            info.genome_sequences[genome_id].end(),
            0UL, [&info](auto a, auto b) { return a + info.computation_costs[b].total_visited; }) << endl;

    scores current_scores = computeScores(info, info.genome_sequences[genome_id], step_size);

    jclass scores_out_cl = env->GetObjectClass(out_scores);

    env->SetIntField(out_scores, env->GetFieldID(scores_out_cl, "scoresCount", "I"), current_scores.non_zero.size());
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "scoresMaxMappings", "[I"),
            create_jni_int_array(env, (jint*)current_scores.flat_map.data(), current_scores.flat_map.size()));
    assert(sizeof(jint) == sizeof(seq_id_t));

    unsigned int cells_count = current_scores.non_zero.size();

    vector<jint> tmp_integer(cells_count);
    vector<jfloat> tmp_floating(cells_count);

    for (int i = 0; i < cells_count; i++) tmp_floating[i] = current_scores.non_zero[i].score;
    jfloatArray cells_scores_arr = create_jni_float_array(env, tmp_floating.data(), tmp_floating.size());

    for (int i = 0; i < cells_count; i++) tmp_floating[i] = current_scores.non_zero[i].perc;
    jfloatArray cells_percs_arr = create_jni_float_array(env, tmp_floating.data(), tmp_floating.size());

    for (int i = 0; i < cells_count; i++) tmp_floating[i] = current_scores.non_zero[i].tr_perc;
    jfloatArray cells_trpercs_arr = create_jni_float_array(env, tmp_floating.data(), tmp_floating.size());


    for (int i = 0; i < cells_count; i++) {  tmp_integer[i] = current_scores.non_zero[i].x; }
    jintArray cells_row_arr = create_jni_int_array(env, tmp_integer.data(), tmp_integer.size());

    for (int i = 0; i < cells_count; i++) {  tmp_integer[i] = current_scores.non_zero[i].y; }
    jintArray cells_column_arr = create_jni_int_array(env, tmp_integer.data(), tmp_integer.size());


    for (int i = 0; i < cells_count; i++) { tmp_integer[i] = get_genome_index(info, current_scores.non_zero[i].x); }
    jintArray first_genome_index_arr = create_jni_int_array(env, tmp_integer.data(), tmp_integer.size());

    for (int i = 0; i < cells_count; i++) { tmp_integer[i] = get_genome_index(info, current_scores.non_zero[i].y); }
    jintArray second_genome_index_arr = create_jni_int_array(env, tmp_integer.data(), tmp_integer.size());

    jobjectArray max_scores_arr = env->NewObjectArray(rows_count, env->FindClass("[F"),
                                                      nullptr);

    vector<jfloat> tmp_floating2 = vector<jfloat>(info.sequences_count);
    for (int i = 0; i < info.sequences_count; i++) tmp_floating2[i] = current_scores.col_max_scores[i];
    jfloatArray col_max_scores_arr = create_jni_float_array(env, tmp_floating2.data(), tmp_floating2.size());

    for (unsigned int row = 0; row < rows_count; row++) {
        env->SetObjectArrayElement(max_scores_arr, row,
                                   create_jni_float_array(env,
                                                          current_scores.max_scores.data() +
                                                          row * info.genomes_count,
                                                          info.genomes_count));
    }

    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "scores", "[F"), cells_scores_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "percs", "[F"), cells_percs_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "tr_percs", "[F"), cells_trpercs_arr);

    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "row", "[I"), cells_row_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "column", "[I"), cells_column_arr);

    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "first_seq_genome", "[I"), first_genome_index_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "second_seq_genome", "[I"), second_genome_index_arr);

    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "max_genome_score", "[[F"), max_scores_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "max_genome_score_col", "[F"), col_max_scores_arr);
}