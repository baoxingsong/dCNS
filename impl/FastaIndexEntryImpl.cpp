//
// Created by bs674 on 8/23/19.
//

#include "FastaIndexEntryImpl.h"


void FastaIndexEntryImpl::readFastaIndexFile( std::string & fastaIndexFileLocation ){
    std::ifstream infile(fastaIndexFileLocation);
    if( ! infile.good()){
        std::cerr << "error in opening fasta index file " << fastaIndexFileLocation << "." << std::endl << " Please index your fasta file using samtools samtools faidx command." << std::endl;
        exit (1);
    }


    std::regex fastaIndexRegex ("^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");
    std::smatch match;
    std::string line;
    while (std::getline(infile, line)){
        regex_search(line, match, fastaIndexRegex);
        if (!match.empty()) {
            //std::cout << "match" << line << std::endl;
            FastaIndexEntry fastaIndexEntry (match[1],  std::stoi(match[2]), std::stol(match[3]), std::stoi(match[4]), std::stoi(match[5]));
            this->entries[match[1]] = fastaIndexEntry;
        }
    }
    infile.close();
}




std::vector<std::string> &FastaIndexEntryImpl::getNames()  {
    return names;
}

void FastaIndexEntryImpl::setNames(const std::vector<std::string> &names) {
    FastaIndexEntryImpl::names = names;
}

std::map<std::string, FastaIndexEntry> &FastaIndexEntryImpl::getEntries()  {
    return entries;
}

void FastaIndexEntryImpl::setEntries(const std::map<std::string, FastaIndexEntry> &entries) {
    FastaIndexEntryImpl::entries = entries;
}
