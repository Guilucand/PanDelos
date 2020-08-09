//
// Created by andrea on 09/08/20.
//

#ifndef NATIVE_PANGENE_IDATA_H
#define NATIVE_PANGENE_IDATA_H

#include <vector>
#include <string>

struct pangene_idata {
    std::vector<std::string> sequences_names;
    std::vector<std::string> sequences_desc;
    std::vector<std::string> sequences;

    std::vector<std::string> genomes;
    std::vector<size_t> seq_to_genomes;

    pangene_idata(std::string const& input);
};



#endif //NATIVE_PANGENE_IDATA_H
