//
// Created by andrea on 09/08/20.
//

#include "stats.h"
#include <iostream>

using namespace std;

template<class T>
static void print_freqs(T &collection) {
    map<size_t, size_t> freqs;
    for (auto &item : collection) {
        freqs[item.second]++;
    }

    for (auto &freq : freqs) {
        cout << freq.first << "\t" << freq.second << endl;
    }
}

void print_undirected_degree_distribution(std::map<std::pair<seq_id_t, seq_id_t>, float> const& network) {

    map<seq_id_t, size_t> links;

    for (auto &element : network) {
        seq_id_t first = element.first.first;
        seq_id_t second = element.first.second;

        if (second > first && network.count({second, first}) > 0) {
            continue;
        }

        links[first]++;
        links[second]++;
    }

    print_freqs(links);
}

void print_directed_degree_distribution(std::map<std::pair<seq_id_t, seq_id_t>, float> const& network) {
    map<seq_id_t, size_t> links;

    for (auto &element : network) {
        seq_id_t first = element.first.first;
        links[first]++;
    }

    print_freqs(links);
}

static vector<seq_id_t> P;

static seq_id_t find(seq_id_t x) {
    if (P[x] == (seq_id_t)-1) return x;
    return P[x] = find(P[x]);
}

static void _union(seq_id_t a, seq_id_t b) {
    a = find(a);
    b = find(b);
    if (a != b) {
        P[a] = b;
    }
}

void print_undirected_connected_components(const map<std::pair<seq_id_t, seq_id_t>, float> &network, size_t max_seq) {
    P.clear();
    P.resize(max_seq, -1);

    for (auto &element : network) {
        _union(element.first.first, element.first.second);
    }

    map<seq_id_t, size_t> coco_sizes;

    for (seq_id_t seq = 0; seq < max_seq; seq++) {
        coco_sizes[find(seq)]++;
    }

    print_freqs(coco_sizes);
}
