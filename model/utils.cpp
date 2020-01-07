//
// Created by bs674 on 8/23/19.
//

#include "utils.h"


// split a string by specific seperator, and generate a vector
void split(const std::string &s, char& delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (item.length() > 0) {
            elems.push_back(item);
        }
    }
}

