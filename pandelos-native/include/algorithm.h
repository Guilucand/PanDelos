//
// Created by andrea on 09/08/20.
//

#ifndef NATIVE_ALGORITHM_H
#define NATIVE_ALGORITHM_H

#include <pangene_idata.h>

typedef uint64_t rank_t;
typedef uint32_t seq_id_t;
typedef uint32_t gen_id_t;

struct mat_cell {
    float score;
    float perc;
    float tr_perc;
    unsigned int row;
    unsigned int col;
};

struct scores {
    std::vector<seq_id_t> sequences;
    std::vector<seq_id_t> flat_map;
    std::vector<mat_cell> non_zero;
    std::vector<float> max_scores;
    std::vector<float> col_max_scores;
};

struct kmer_rank {
    rank_t rank;
    seq_id_t seq;
    unsigned int count;

    inline kmer_rank() {}

    inline kmer_rank(rank_t rank, seq_id_t seq, unsigned int count) {
        this->rank = rank;
        this->seq = seq;
        this->count = count;
    }

    inline bool dupl(kmer_rank &other) {
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

    inline seq_computation_cost() : total_visited(0), total_seq_length(0), linear_ratio(1.0) {}

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
    std::vector<kmer_rank> kmers;
    std::vector<unsigned int> kseq_lengths;
    std::vector<gen_id_t> seq_gen_mapping;
    std::vector<std::vector<seq_id_t>> genome_sequences;
    std::vector<std::vector<kmers_range>> kmers_ranges;
    std::vector<seq_computation_cost> computation_costs;
};

scores compute_scores(pair_info const& info, std::vector<seq_id_t> const& sequences, uint32_t step_size);

pair_info preprocess_sequences(pangene_idata const& data, uint32_t kvalue, bool onlyComplexity);


#endif //NATIVE_ALGORITHM_H
