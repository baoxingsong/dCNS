//
// Created by bs674 on 8/23/19.
//

#ifndef AND_CNS_FASTAINDEXENTRYIMPL_H
#define AND_CNS_FASTAINDEXENTRYIMPL_H

#include <string>
#include <vector>
#include <map>
#include <regex>
#include "../model/model.h"

class FastaIndexEntryImpl {
    private:
        std::vector<std::string> names; // the vector could keep the order
        std::map<std::string, FastaIndexEntry> entries;
    public:
        void readFastaIndexFile( std::string & fastaIndexFileLocation );
//        void createFastaIndexFile( std::string & fastaFileLocation );

        std::vector<std::string> &getNames();
        void setNames(const std::vector<std::string> &names);
        std::map<std::string, FastaIndexEntry> &getEntries();
        void setEntries(const std::map<std::string, FastaIndexEntry> &entries);
};


#endif //AND_CNS_FASTAINDEXENTRYIMPL_H
