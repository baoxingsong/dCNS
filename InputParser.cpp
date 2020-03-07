/*
 * =====================================================================================
 *
 *       Filename:  InputParser.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 01:13:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include "InputParser.h"




bool str2bool( std::string p, bool defaultValue ){
    std::transform(p.begin(), p.end(), p.begin(),::toupper);
    if( p.compare("YES")==1 || p.compare("TRUE")==1 || p.compare("Y")==1 || p.compare("T")==1 || p.compare("1")==1){
        return true;
    }else if ( p.compare("NO")==1 || p.compare("FALSE")==1 || p.compare("F")==1 || p.compare("N")==1 || p.compare("0")==1){
        return true;
    }
    return defaultValue;
}

InputParser::InputParser (int &argc, char **argv){
    for (int i=1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

std::string InputParser::getCmdOption( std::string &option) {
    std::vector<std::string>::iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }
    return "";
}

std::vector<std::string> InputParser::getCmdOptionMultipleParameters( const char* o) {
    std::string option = o;
    std::vector<std::string> sp;

    std::vector<std::string>::iterator itr;
    itr = std::find(this->tokens.begin(), this->tokens.end(), option);
    while (itr != this->tokens.end() && ++itr != this->tokens.end() && (*itr)[0] != '-' && (*itr)[0] != '>' ) {
        sp.push_back(*itr);
    }
    return sp;
}

std::string InputParser::getCmdOption( const char* o) {
    std::string option = o;
    return getCmdOption(option);
}

bool InputParser::cmdOptionExists(std::string &option) {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
           != this->tokens.end();
}

bool InputParser::cmdOptionExists( const char* o){
    std::string option = o;
    return cmdOptionExists(option);
}

void usage( ){
    std::string progName = "dCNS";
    std::cout << "Program: " << progName << " version: " << VERSION << std::endl <<
    "Usage: "<<progName<<" <command> [options]"<< std::endl <<
    "Commands:"<< std::endl <<
    "    slideWindow         use a sliding window to align against the reference sequence" << std::endl <<
    "    pairCnsXExtend      sliding window seed, x-drop extension" << std::endl <<
//    "    pairCns2Gaps        sliding window seed, x-drop extension, 2-piece gap cost extension" << std::endl <<
//    "    weighted1Gap        sliding window seed, x-drop extension, weighted x-drop re-alignment and extension" << std::endl <<
//    "    weighted2Gaps       sliding window seed, x-drop extension, weighted x-drop 2-piece gap cost re-alignment and extension" << std::endl <<
    "    cut1Gap             masking genome, sliding window seed, x-drop extension" << std::endl <<
//    "    cut2Gaps            masking genome sliding window seed, x-drop extension, 2-piece gap cost extension" << std::endl <<
    "    cut2Gaps2           extend the cut1Gap result using cut2Gaps approach" << std::endl <<
    "    multCns             detect conserved fragment among a group of sequences" << std::endl <<
    std::endl <<
    "Complementary commands:"<< std::endl <<
    "    maskGenome          masking genome for high repeat regions" << std::endl <<
    "    ranSco              get the smith-waterman scores using random sequence fragments" << std::endl <<
    "                        for parameters fitting and lambda and k value calculating aim" << std::endl;


    //"    pair-cns        detect conserved fragment between a pair of sequences" << std::endl <<

}

