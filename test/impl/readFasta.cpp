//
// Created by Baoxing song on 2019-01-03.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>

#include "../../InputParser.h"
#include "../../model/model.h"
#include "../../controlLayer.h"

#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>

// most general one
TEST(readFasta, c1){ // just to make sure that every line has been analysed

    std::string refFasta = "/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa";
    std::string chr = "1";
    std::cout << std::endl;


    for( int i=0; i<1000; ++i ) {
        std::map<std::string, std::string> sequences0;
        std::vector<std::string> seqNames0;
        readFastaFile( refFasta, sequences0, seqNames0);
        std::string seq1 = getSubSequence( sequences0[chr], 101,  203,  1);
        std::cout << "seq1:" << seq1 << std::endl;
    }
    /*
    for( int i=0; i<1000; ++i ){
        FastaWithIndex fastaWithIndex(refFasta);
        std::string seq2 = fastaWithIndex.getSubSequence( chr, 101,  203,  1);
        std::cout << "seq2:" << seq2 << std::endl;
    }
*/
}

