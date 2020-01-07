//
// Created by bs674 on 8/23/19.
//


#ifndef AND_CNS_FASTAWITHINDEX_H
#define AND_CNS_FASTAWITHINDEX_H

#include <fstream>//required for file operation:
#include <iostream>//required for input/output:
#include <string>

#include "FastaIndexEntryImpl.h"
#include "fasta.h"
class FastaWithIndex {
    private:
        FastaIndexEntryImpl fastaIndexEntryImpl;
        std::ifstream inFile;
    public:
        FastaIndexEntryImpl &getFastaIndexEntryImpl();
        void setFastaIndexEntryImpl(const FastaIndexEntryImpl &fastaIndexEntryImpl);
        FastaWithIndex(const std::string & fastaFileLocation);
        std::string getSubSequence(std::string & chr, int32_t start1, int32_t end1, int8_t strand) ;
        ~FastaWithIndex();
};


#endif //AND_CNS_FASTAWITHINDEX_H
