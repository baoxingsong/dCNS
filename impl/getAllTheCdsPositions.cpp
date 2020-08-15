//
// Created by Baoxing Song on 2019-08-07.
//

#include "getAllTheCdsPositions.h"


// read gff file tell where is CDS, this is for genome masking purpose
void gffToMask (const std::string& filePath, std::map<std::string, std::string>& genome,
                 std::map<std::string, std::string> & ifCds){

    ifCds.clear();
    for( std::map<std::string, std::string>::iterator it = genome.begin(); it!=genome.end(); ++it ){
        ifCds[it->first] = std::string(it->second.size(), '0');
    }

    std::set<std::string> cdssNickNames;
    cdssNickNames.insert("CDS");
    cdssNickNames.insert("start_codon");
    cdssNickNames.insert("stop_codon");

    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit (1);
    }
    int32_t start;
    int32_t end;
    int32_t temp;
    int32_t length;
    char splim='\t';
    std::string line;
    std::string chromosomeName;
    while (std::getline(infile, line)){
        if(line[0]=='#' || line.size()<9){
        }else {
            std::vector<std::string> elemetns;
            split(line, splim, elemetns);
            if(elemetns.size()==9){
                start = stoi(elemetns[3]);
                end = stoi(elemetns[4]);
                if (start > end) {
                    temp = start;
                    start = end;
                    end = temp;
                }
                --start; // change the coordinates to 0 based
                length = end - start;
                chromosomeName = elemetns[0];
                if ( cdssNickNames.find(elemetns[2]) != cdssNickNames.end() ) {
                    ifCds[chromosomeName].replace(start, length, length, '1');
                }
            }
        }
    }
}


void gffToMaskGene (const std::string& filePath, std::map<std::string, std::string>& genome,
               std::map<std::string, std::string> & ifCds){

    ifCds.clear();
    for( std::map<std::string, std::string>::iterator it = genome.begin(); it!=genome.end(); ++it ){
        ifCds[it->first] = std::string(it->second.size(), '0');
    }

    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit (1);
    }
    int32_t start;
    int32_t end;
    int32_t temp;
    int32_t length;
    char splim='\t';
    std::string line;
    std::string chromosomeName;
    while (std::getline(infile, line)){
        if(line[0]=='#' || line.size()<9){
        }else {
            std::vector<std::string> elemetns;
            split(line, splim, elemetns);
            if(elemetns.size()==9){
                start = stoi(elemetns[3]);
                end = stoi(elemetns[4]);
                if (start > end) {
                    temp = start;
                    start = end;
                    end = temp;
                }
                --start; // change the coordinates to 0 based
                length = end - start;
                chromosomeName = elemetns[0];
                if ( elemetns[2].compare("gene") == 0 && elemetns[8].find("biotype=protein_coding") != std::string::npos ) {
                    ifCds[chromosomeName].replace(start, length, length, '1');
                }
            }
        }
    }
}



void samToMask (const std::string& filePath, std::map<std::string, std::string>& genome,
                std::map<std::string, std::string> & ifCds, double & similarity, std::map<std::string, std::string> & cdsSequences){
    bool recreatIfCds = false;
    if (ifCds.size() ==0){
        recreatIfCds = true;
    }
    for( std::map<std::string, std::string>::iterator it = genome.begin(); it!=genome.end(); ++it ){
        if(ifCds[it->first].size() != it->second.size()){
            recreatIfCds = true;
        }
    }

    if (recreatIfCds){
        ifCds.clear();
        for( std::map<std::string, std::string>::iterator it = genome.begin(); it!=genome.end(); ++it ){
            ifCds[it->first] = std::string(it->second.size(), '0');
        }
        std::cout << "GFF was not used for genome masking" << std::endl;
    }
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit (1);
    }
    std::regex samRegex("^(.*?)\\s+(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\w+)\\t.*NM:i:(\\d+).*cs:(\\S+)");
    std::smatch samMatch;
    std::regex cigarRegex ("(\\d+)([MIDNSH])");
    std::smatch cigarMatch;

    std::regex scRegex (":(\\d+)");
    std::smatch scMatch;

    int32_t start;
    int32_t length;
    std::string line;
    std::string cigarString;
    std::string cigarChar;
    std::string chromosomeName;
    while (std::getline(infile, line)){
        if(line[0]=='@' || line.size()<9){
        }else {
            regex_search(line, samMatch, samRegex);
            if (!samMatch.empty()) {
                double match = 0;
                std::string sc = samMatch[8];
                while ( std::regex_search(sc, scMatch, scRegex) ){
                    match += stoi(scMatch[1]);
                    sc = scMatch.suffix().str();
                }
                if( match >= similarity * cdsSequences[samMatch[1]].size() ){
                    chromosomeName = samMatch[3];
                    start = stoi(samMatch[4]);
                    cigarString = samMatch[6];
                    while ( std::regex_search(cigarString, cigarMatch, cigarRegex) ){
                        cigarChar = cigarMatch[2];
                        length = stoi(cigarMatch[1]);
                        if( cigarChar[0] == 'M' || cigarChar[0] == 'D' ){
                            ifCds[chromosomeName].replace(start-1, length, length, '1');
                            start += length;
                        }else if ( cigarChar[0] == 'N'){
                            start += length;
                        }
                        cigarString = cigarMatch.suffix().str();
                    }
                }
            }
        }
    }
}
