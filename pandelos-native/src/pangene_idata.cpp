//
// Created by andrea on 09/08/20.
//

#include "pangene_idata.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_map>

using namespace std;

// Taken from: https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

pangene_idata::pangene_idata(std::string const& input) {
    ifstream input_file(input);

    if (!input_file.is_open()) {
        cout << "Cannot open file '" << input << "'" << endl;
    }

    auto genomes_map = unordered_map<string, size_t>();

    while (!input_file.eof()) {
        string gen_name;
        string seq_name;
        string product;
        string sequence;

        getline(input_file, gen_name, '\t');
        getline(input_file, seq_name, '\t');
        getline(input_file, product);

        getline(input_file, sequence);

        trim(sequence);

        if (sequence.length() > 0) {

            if (genomes_map.count(gen_name) == 0) {
                genomes_map.insert({ gen_name, this->genomes.size() });
                this->genomes.push_back(gen_name);
            }

            this->sequences_names.push_back(seq_name);
            this->sequences_desc.push_back(product);
            this->seq_to_genomes.push_back(genomes_map[gen_name]);
            this->sequences.push_back(sequence);
        }
    }
}
