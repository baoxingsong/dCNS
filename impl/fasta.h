//
// Created by Baoxing song on 2019-01-02.
//

#ifndef SONG_CNS_READFASTA_H
#define SONG_CNS_READFASTA_H

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <vector>
#include <algorithm>
#include <regex>



class FastaMeta{
private:
    std::string species;
    std::string chr;
    int32_t start;
    int32_t end;
    int strand;
public:
    FastaMeta();
    FastaMeta(const std::string &species, const std::string &chr, int32_t & start, int32_t & end, int & strand);
    const std::string &getSpecies() const;
    void setSpecies(const std::string &species);
    const std::string &getChr() const;
    void setChr(const std::string &chr);
    int32_t getStart() const;
    void setStart(int32_t start);
    int32_t getEnd() const;
    void setEnd(int32_t end);
    int getStrand() const;
    void setStrand(int strand);

};


void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences, std::vector<std::string> & seqNames);
void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences);
void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences,
                    std::map<std::string, FastaMeta>& metaInformations, std::vector<std::string> & seqNames );
std::string getReverseComplementary(const std::string& sequence);
std::string getSubSequence( std::string sequence, const int& _start, const int& _end);
std::string getSubSequence( std::string & sequence, const int& _start, const int& _end, const int& strand);


#endif //SONG_CNS_READFASTA_H
