//
// Created by Baoxing song on 2019-01-03.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>

TEST(SequenceCharToUInt8, c1){ // just to make sure that every line has been analysed
    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq = "ATGCTGTGGCTGTCGATGCTGACTAE";
    int8_t * seq_int8 = sequenceCharToUInt8.seq_to_int8(seq);
    std::string seq_recover = sequenceCharToUInt8.int8_to_seq(seq_int8, seq.length());
    std::cout << "seq         " << seq << std::endl;
    std::cout << "seq_recover " << seq_recover << std::endl;
    int8_t * seq_int8_reCom = sequenceCharToUInt8.rev_comp(seq_int8, seq.length());
    std::string seq_recCom = sequenceCharToUInt8.int8_to_seq(seq_int8_reCom, seq.length());
    std::cout << "seq_recCom  " << seq_recCom << std::endl;
    ASSERT_EQ(0, 0);
}
