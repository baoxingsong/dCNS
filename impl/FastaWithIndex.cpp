//
// Created by bs674 on 8/23/19.
//


/**
 * do not try to creat fasta index file using this software, since we maybe need to read the fasta with a lot of process, that would cause problem
 * */
#include "FastaWithIndex.h"

FastaWithIndex::FastaWithIndex ( const std::string & fastaFileLocation) {
    std::string fastaIndexFilePath = fastaFileLocation + ".fai";
    this->fastaIndexEntryImpl.readFastaIndexFile(  fastaIndexFilePath );
    inFile.open(fastaFileLocation, std::ios::in);
}

std::string FastaWithIndex::getSubSequence(std::string & chr, int32_t start1, int32_t end1,
        int8_t strand) { // it should be OK without synchronized
    int start;
    int end;
    if (start1 < end1) {
        start = start1;
        end = end1;
    } else {
        start = end1;
        end = start1;
    }
//    for( std::map<std::string, FastaIndexEntry>::iterator i = fastaIndexEntryImpl.getEntries().begin();
//        i != fastaIndexEntryImpl.getEntries().end(); ++i ){
//        std::cout << i->first << std::endl;
//    }

    std::string subSequence;
    if (fastaIndexEntryImpl.getEntries().find(chr) != fastaIndexEntryImpl.getEntries().end()) {
        FastaIndexEntry entry = fastaIndexEntryImpl.getEntries()[chr];
        start--;
        if (start < 0) {
            start = 0;
        }

        if (entry.getLength() < end) {
            end = entry.getLength();
        }
        if (entry.getLength() < start) {
            start = entry.getLength();
        }
        if (end < 0) {
            end = 0;
        }

        int length = (end - start);

        int newlines_before = start > 0 ? (start - 1) / entry.getLineBlen() : 0;
        int newlines_by_end = (start + length - 1) / entry.getLineBlen();
        int newlines_inside = newlines_by_end - newlines_before;
        int seqlen = length + newlines_inside;

        inFile.seekg(entry.getOffset() + newlines_before + start, std::ios::beg);

        char text;
        for( int32_t  i = 0; i<seqlen && !inFile.eof();++i ){
            inFile.get(text); //read data:
            subSequence += text;
        }
        subSequence.erase(remove_if(subSequence.begin(), subSequence.end(), isspace), subSequence.end());
    }
    if( strand == 0 ){
        return getReverseComplementary(subSequence);
    }
    return subSequence;
}


FastaWithIndex::~FastaWithIndex(){
    inFile.close();
}

FastaIndexEntryImpl &FastaWithIndex::getFastaIndexEntryImpl() {
    return fastaIndexEntryImpl;
}

void FastaWithIndex::setFastaIndexEntryImpl(const FastaIndexEntryImpl &fastaIndexEntryImpl) {
    FastaWithIndex::fastaIndexEntryImpl = fastaIndexEntryImpl;
}
