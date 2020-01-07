//
// Created by Baoxing Song on 2019-08-07.
//

#ifndef AND_CNS_MASKGENOME_H
#define AND_CNS_MASKGENOME_H


#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>

void maskGenome( std::map<std::string, std::string>& genome, char & alterChar, std::string & kmerProfileFile,
                 int32_t & maskKmerFrequency, int32_t & outputLineWidth, std::map<std::string, std::string> & ifCds,
                 std::string & outputFile, std::vector<std::string> & seqNames, bool & ifSoftmask );

#endif //AND_CNS_MASKGENOME_H
