//
// Created by bs674 on 8/9/19.
//

#ifndef AND_CNS_SAMFILETOSEEDS_H
#define AND_CNS_SAMFILETOSEEDS_H

#include <string>
#include <vector>
#include <map>
#include <regex>
#include "../model/model.h"

void samFileToSeeds ( std::string & samFile, std::map<std::string, std::map<std::string, std::vector<Seed>>> & positiveSeeds,
                      std::map<std::string, std::map<std::string, std::vector<Seed>>> & negativeSeeds,
                      std::map<std::string, std::string> & refGenome, std::map<std::string, std::string> & queryGenome);


#endif //AND_CNS_SAMFILETOSEEDS_H
