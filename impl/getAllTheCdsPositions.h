//
// Created by Baoxing Song on 2019-08-07.
//

#ifndef AND_CNS_GETALLTHECDSPOSITIONS_H
#define AND_CNS_GETALLTHECDSPOSITIONS_H

#include "gffToCategory.h"
void gffToMask (const std::string& filePath, std::map<std::string, std::string>& genome,
                std::map<std::string, std::string> & ifCds);
void gffToMaskGene (const std::string& filePath, std::map<std::string, std::string>& genome,
                    std::map<std::string, std::string> & ifCds);

void samToMask (const std::string& filePath, std::map<std::string, std::string>& genome,
                std::map<std::string, std::string> & ifCds, double & similarity, std::map<std::string, std::string> & cdsSequences);
#endif //AND_CNS_GETALLTHECDSPOSITIONS_H
