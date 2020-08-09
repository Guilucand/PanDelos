//
// Created by andrea on 09/08/20.
//

#ifndef NATIVE_STATS_H
#define NATIVE_STATS_H

#include <map>
#include "pangene_idata.h"
#include "algorithm.h"

void print_undirected_degree_distribution(std::map<std::pair<seq_id_t, seq_id_t>, float> const& network);
void print_directed_degree_distribution(std::map<std::pair<seq_id_t, seq_id_t>, float> const& network);
void print_undirected_connected_components(std::map<std::pair<seq_id_t, seq_id_t>, float> const& network, size_t max_seq);

#endif //NATIVE_STATS_H
