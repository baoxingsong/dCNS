/*
 * =====================================================================================
 *
 *       Filename:  gffToCategory.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/04/2019 19:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song, songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************
generate weigth vector by reading gff file



 ************************************************************************/


#ifndef SSW_GFFTOCATEGORY_H
#define SSW_GFFTOCATEGORY_H


#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <regex>

// and gff records
struct element{
    int32_t start;
    int32_t end;
};


void split(const std::string &s, char& delim, std::vector<std::string> &elems);

void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, int32_t> & chrSize, const std::string & outputFile);
void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, int32_t> & chrSize, std::map<std::string, int16_t *> & categories);
void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, int32_t> & chrSize, std::map<std::string, int16_t *> & categories, std::map<std::string, int16_t *> & categories_rev);
#endif //SSW_GFFTOCATEGORY_H
