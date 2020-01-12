#include "pangene_native.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <queue>
#include <cassert>
#include <cstring>

using namespace std;

template <class T>
using min_heap = priority_queue<T, vector<T>, greater<T>>;

typedef unsigned long rank_t;
typedef unsigned int seq_id_t;

vector<vector<rank_t>> kmers;

unsigned char rank_values[256];
unsigned char rank_base;
unsigned char rank_byte_order;
unsigned int kvalue;
rank_t max_rank = 1;

rank_t update_rank(rank_t current, jchar next) {
    return (current * rank_base + rank_values[next]) % max_rank;
}

void rank_init(int kvalue) {
    ::kvalue = kvalue;
    for (int i = 'A'; i <= 'Z'; i++) {
        rank_values[i] = i - 'A';
    }
    rank_base = 'Z' - 'A' + 1;
    max_rank = 1;
    while (kvalue--) max_rank *= rank_base;
    auto rank_tmp = max_rank;

    rank_byte_order = 0;
    while (rank_tmp) {
        rank_byte_order++;
        rank_tmp /= 256;
    }
}

vector<rank_t> do_ranking(const jchar *chars, jsize len) {

    vector<rank_t> result(len - kvalue + 1);
    rank_t rank = 0;

    for (jsize i = 0; i < kvalue; i++) {
        rank = update_rank(rank, chars[i]);
    }
    result[0] = rank;

    for (jsize i = kvalue; i < len; i++) {
        rank = update_rank(rank, chars[i]);
        result[i - kvalue + 1] = rank;
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

void Java_infoasys_cli_pangenes_PangeneNative_initialize(JNIEnv *env, jobject obj, jint kvalue) {
    rank_init(kvalue);
}

void Java_infoasys_cli_pangenes_PangeneNative_preprocessSequences(JNIEnv *env, jobject obj, jobject data) {
    jclass dataClass = env->GetObjectClass(data);
    jfieldID seqId = env->GetFieldID(dataClass, "sequences", "Ljava/util/Vector;");

    jobject objSeq = env->GetObjectField(data, seqId);

    jclass vecClass = env->GetObjectClass(objSeq);
    jmethodID elMethod = env->GetMethodID(vecClass, "get", "(I)Ljava/lang/Object;");
    jmethodID szMethod = env->GetMethodID(vecClass, "size", "()I");

    int seq_number = env->CallIntMethod(objSeq, szMethod);

    for (int i = 0; i < seq_number; i++) {
        auto str = (jstring) env->CallObjectMethod(objSeq, elMethod, i);

        jboolean iscopy = false;
        jsize len = env->GetStringLength(str);
        const jchar *seq = env->GetStringChars(str, &iscopy);

        kmers.push_back(do_ranking(seq, len));
        for (int b = 0; b < rank_byte_order; b++) {
            counting_sort(kmers.back(), b);
        }

        env->ReleaseStringChars(str, seq);
    }
}

struct kmers_range {
    pair<pair<rank_t, seq_id_t>, int> *start;
    pair<pair<rank_t, seq_id_t>, int> *current;
    pair<pair<rank_t, seq_id_t>, int> *end;
};

struct pair_info {
    vector<int> flat_map;
    vector<int> rev_flat_map;
    int first_genome_seq_count;
    vector<pair<pair<rank_t, seq_id_t>, int>> kmers_list;
    vector<vector<kmers_range>> kmers_ranges;
};

jlong Java_infoasys_cli_pangenes_PangeneNative_computePair(JNIEnv *env, jobject obj, jintArray first, jintArray second) {

    pair_info *info = new pair_info;

    jboolean is_copy = false;
    jsize first_sz = env->GetArrayLength(first);
    jsize second_sz = env->GetArrayLength(second);

    jint *first_el = env->GetIntArrayElements(first, &is_copy);
    jint *second_el = env->GetIntArrayElements(second, &is_copy);

    vector<jint> g1_seqs(first_el, first_el + first_sz);
    vector<jint> g2_seqs(second_el, second_el + second_sz);

    info->first_genome_seq_count = first_sz;
    info->flat_map.resize(kmers.size());
    vector<jint> &flat_map = info->flat_map;

    env->ReleaseIntArrayElements(first, first_el, JNI_ABORT);
    env->ReleaseIntArrayElements(second, second_el, JNI_ABORT);

    vector<jint> all_seqs(g1_seqs);
    all_seqs.insert(all_seqs.end(), g2_seqs.begin(), g2_seqs.end());
// TODO: Ensure that all sequences are unique, else collapse them
//    sort(all_seqs.begin(), all_seqs.end());
//    all_seqs.resize(distance(all_seqs.begin(), unique(all_seqs.begin(), all_seqs.end())));


    vector<int> positions(kmers.size(), 0);
    min_heap<pair<rank_t, seq_id_t>> mheap;

    info->rev_flat_map.resize(all_seqs.size());

    for (auto seq_index : all_seqs) {
        if (kmers[seq_index].size() > 0) {
            flat_map[seq_index] = mheap.size();
            info->rev_flat_map[mheap.size()] = seq_index;
            mheap.push({kmers[seq_index][0], seq_index});
        }
    }

    vector<pair<pair<rank_t, seq_id_t>, int>> &kmers_cross = info->kmers_list;
    int allowed_seq_count = mheap.size();

    while (!mheap.empty()) {
        auto top = mheap.top(); mheap.pop();
        if (kmers_cross.size() && kmers_cross.back().first == top) {
            kmers_cross.back().second++;
        }
        else {
            kmers_cross.push_back({top, 1});
        }
        rank_t rank = top.first;
        seq_id_t seq_id = top.second;
        positions[seq_id]++;
        if (positions[seq_id] < kmers[seq_id].size()) {
            mheap.push({ kmers[seq_id][positions[seq_id]], seq_id });
        }
    }

    info->kmers_ranges.resize(allowed_seq_count);
    vector<vector<kmers_range>> &kmer_ranges = info->kmers_ranges;

    rank_t current_rank = kmers_cross.begin()->first.first;
    int current_rank_start = 0;
    for (int i = 0; i < kmers_cross.size(); i++) {
        if (current_rank != kmers_cross[i].first.first || (i == kmers_cross.size() - 1)) {

            // Process current rank
            int prev_rank_begin = current_rank_start;
            int prev_rank_end = i;

            for (int j = prev_rank_begin; j < prev_rank_end; j++) {
                kmer_ranges[flat_map[kmers_cross[j].first.second]].push_back(
                        kmers_range {
                            .start = kmers_cross.data() + prev_rank_begin,
                            .current = kmers_cross.data() + j,
                            .end = kmers_cross.data() + prev_rank_end
                        });
            }

            current_rank_start = i;
            current_rank = kmers_cross[i].first.first;
        }
    }
    //    accumulate(first_el, first_el + first_sz, 0, [](int current, jint next) { return current + ; });
    return (long)info;
}

static unsigned int get_genome_index(pair_info *info, int flat_index) {
    return flat_index >= info->first_genome_seq_count;
}

struct mat_cell {
    float score;
    float perc;
    float tr_perc;
    int x;
    int y;
};

jbooleanArray create_jni_boolean_array(JNIEnv *env, jboolean *ptr, jsize len) {
    jbooleanArray arr = env->NewBooleanArray(len);
    env->SetBooleanArrayRegion(arr, 0, len, ptr);
    return arr;
}

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

void
Java_infoasys_cli_pangenes_PangeneNative_computeSequenceScores(JNIEnv *env, jobject obj, jlong addr_ref, jint row_start, jint row_end, jobject out_scores) {
    pair_info *info = (pair_info *) addr_ref;
    int seqs_count = info->rev_flat_map.size();

    // Non-zero cells list
    vector<mat_cell> nonzero_cells;

    /* Global max intra and inter scores */
    vector<float> max_intra_scores = vector<float>(seqs_count, 0);
    vector<float> max_inter_scores = vector<float>(seqs_count, 0);

    /* Single row data, overwritten on each row cycle  */
    vector<int> row_intersection_size = vector<int>(seqs_count, 0);
    vector<int> row_connection_perc_cnt = vector<int>(seqs_count, 0);
    vector<int> row_transposed_perc_cnt = vector<int>(seqs_count, 0);

    for (int row = row_start; row < row_end; row++) {

        // Avoid computing duplicate results
        int col_limit = row;

        /* Reset all row data */
        memset(row_intersection_size.data(), 0, row_intersection_size.size() * sizeof(int));
        memset(row_connection_perc_cnt.data(), 0, row_connection_perc_cnt.size() * sizeof(int));
        memset(row_transposed_perc_cnt.data(), 0, row_transposed_perc_cnt.size() * sizeof(int));

        unsigned int my_gindex = get_genome_index(info, row);

        for (auto &range : info->kmers_ranges[row]) {
            auto it = range.start;
            int my_cnt = range.current->second;
            while (it < range.end) {
                row_intersection_size[info->flat_map[it->first.second]] += min(it->second, my_cnt);
                row_connection_perc_cnt[info->flat_map[it->first.second]] += my_cnt;
                row_transposed_perc_cnt[info->flat_map[it->first.second]] += it->second;
                it++;
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
        for (int i = 0; i < seqs_count; i++) {
            int my_kcnt = kmers[info->rev_flat_map[row]].size();
            int other_kcnt = kmers[info->rev_flat_map[i]].size();
            int union_size = my_kcnt + other_kcnt - row_intersection_size[i];
            float perc = (float) row_connection_perc_cnt[i] / (float) my_kcnt;
            float tr_perc = (float) row_transposed_perc_cnt[i] / (float) other_kcnt;
            float threshold = 1.0f / (2.0f * (float) kvalue);
            bool score_valid = perc >= threshold || tr_perc >= threshold;
            float score = (float) row_intersection_size[i] / (float) union_size *
                            (score_valid ? 1.0f : 0.0f);

            // Add to nonzero cells only if score is > 0
            if (score > 0.0f) {

                if (i < col_limit) {
                    nonzero_cells.push_back(mat_cell{
                            .score = score,
                            .perc = perc,
                            .tr_perc = tr_perc,
                            .x = row,
                            .y = i
                    });
                }

                unsigned int other_gindex = get_genome_index(info, i);

                jfloat *max_scores[2] = {
                        &max_intra_scores[row],
                        &max_inter_scores[row]
                };

                // If other and my gindex are equal, update intra score (xor of equal values is always 0)
                // else update the inter score (xor of different values is 1 if one is 0 and the other is 1)
                jfloat &current_score = *max_scores[other_gindex ^ my_gindex];
                current_score = max(current_score, score);
            }
        }
    }

    jfloatArray max_intra_scores_arr = create_jni_float_array(env, max_intra_scores.data(), max_intra_scores.size());
    jfloatArray max_inter_scores_arr = create_jni_float_array(env, max_inter_scores.data(), max_inter_scores.size());


    vector<jboolean> is_same_genome(nonzero_cells.size());
    for (int i = 0; i < is_same_genome.size(); i++) {
        is_same_genome[i] =
                get_genome_index(info, nonzero_cells[i].x) == get_genome_index(info, nonzero_cells[i].y);
    }
    jbooleanArray is_same_genome_arr = create_jni_boolean_array(env, is_same_genome.data(), is_same_genome.size());

    vector<jint> tmpibuffer(nonzero_cells.size());

    for (int i = 0; i < tmpibuffer.size(); i++) tmpibuffer[i] = info->rev_flat_map[nonzero_cells[i].x];
    jintArray cells_firstseq_arr = create_jni_int_array(env, tmpibuffer.data(), tmpibuffer.size());

    for (int i = 0; i < tmpibuffer.size(); i++) tmpibuffer[i] = nonzero_cells[i].x;
    jintArray cells_row_arr = create_jni_int_array(env, tmpibuffer.data(), tmpibuffer.size());

    for (int i = 0; i < tmpibuffer.size(); i++) tmpibuffer[i] = info->rev_flat_map[nonzero_cells[i].y];
    jintArray cells_secondseq_arr = create_jni_int_array(env, tmpibuffer.data(), tmpibuffer.size());

    for (int i = 0; i < tmpibuffer.size(); i++) tmpibuffer[i] = nonzero_cells[i].y;
    jintArray cells_column_arr = create_jni_int_array(env, tmpibuffer.data(), tmpibuffer.size());

    vector<jfloat> tmpfbuffer(nonzero_cells.size());

    for (int i = 0; i < tmpibuffer.size(); i++) tmpfbuffer[i] = nonzero_cells[i].score;
    jfloatArray cells_scores_arr = create_jni_float_array(env, tmpfbuffer.data(), tmpfbuffer.size());

    for (int i = 0; i < tmpibuffer.size(); i++) tmpfbuffer[i] = nonzero_cells[i].perc;
    jfloatArray cells_percs_arr = create_jni_float_array(env, tmpfbuffer.data(), tmpfbuffer.size());

    for (int i = 0; i < tmpibuffer.size(); i++) tmpfbuffer[i] = nonzero_cells[i].tr_perc;
    jfloatArray cells_trpercs_arr = create_jni_float_array(env, tmpfbuffer.data(), tmpfbuffer.size());

    jclass scores_out_cl = env->GetObjectClass(out_scores);

    env->SetIntField(out_scores, env->GetFieldID(scores_out_cl, "scoresCount", "I"), nonzero_cells.size());

    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "sameGenome", "[Z"), is_same_genome_arr);

    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "max_intra_score", "[F"), max_intra_scores_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "max_inter_score", "[F"), max_inter_scores_arr);

    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "scores", "[F"), cells_scores_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "percs", "[F"), cells_percs_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "tr_percs", "[F"), cells_trpercs_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "firstSeqIndex", "[I"), cells_firstseq_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "secondSeqIndex", "[I"), cells_secondseq_arr);

    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "row", "[I"), cells_row_arr);
    env->SetObjectField(out_scores, env->GetFieldID(scores_out_cl, "column", "[I"), cells_column_arr);
}

void Java_infoasys_cli_pangenes_PangeneNative_freePairStruct(JNIEnv *env, jobject obj, jlong addr) {
    pair_info *info = (pair_info*)addr;
    delete info;
}
