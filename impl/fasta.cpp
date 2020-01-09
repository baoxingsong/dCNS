//
// Created by Baoxing song on 2019-01-02.
//

#include "fasta.h"

void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences, std::vector<std::string> & seqNames ){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit (1);
    }
    std::string name;
    std::stringstream sequencestream;
    std::string line="";
    while (std::getline(infile, line)){
        if( line[0] == '>'  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            // this method is faster than regular expression
            if( name.size()>0 ){
                std::string sequence = sequencestream.str();
                 (sequence.begin(), sequence.end(), sequence.begin(),::toupper);
                sequences[name]=sequence;
                seqNames.push_back(name);
            }
            name=line.substr(1, line.find(" ", 0)-1);
            if( name[0] == '>' ){
                name=name.substr(1, name.find("\t", 0)-1);
            }else{
                name=name.substr(0, name.find("\t", 0)-1);
            }
            sequencestream.str(std::string());
        }else{
            sequencestream << line;
        }
    }
    if( name.size()>0 ){
        std::string sequence = sequencestream.str();
        std::transform(sequence.begin(), sequence.end(), sequence.begin(),::toupper);
        if( sequence.size()>0 ){
            sequences[name]=sequence;
            seqNames.push_back(name);
        }
    }
}

void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences,
        std::map<std::string, FastaMeta>& metaInformations, std::vector<std::string> & seqNames ){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit (1);
    }

    std::regex fastaRegex (">(\\S+)\\s");
    std::regex metaRegex ("species:(\\S+)\\s+chr:(\\S+)\\s+strand:(\\S+)\\s+start:(\\S+)\\s+end:(\\S+)");
    std::smatch match;
    std::string name="";
    std::string meta_line="";
    std::stringstream sequencestream;
    std::string line="";
    while (std::getline(infile, line)){
        if( line[0] == '>'  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            // this method is faster than regular expression
            if( name.size()>0 ){
                std::string sequence = sequencestream.str();
                (sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
                sequences[name]=sequence;
//                std::cout << std::endl << "line 68" << name << " " << sequence << std::endl;
                seqNames.push_back(name);
                regex_search(meta_line, match, metaRegex);
//                std::cout << "line 70:" << meta_line << std::endl;
                if (!match.empty()) {
//                    std::cout << "line 72:" << meta_line << std::endl;
                    int strand = std::stoi(match[3]);
                    int32_t start = std::stoi(match[4]);
                    int32_t end = std::stoi(match[5]);
                    FastaMeta fastaMeta(match[1], match[2], start, end, strand);
                    metaInformations[name]=fastaMeta;
                }
            }
            meta_line = line;
            name="";
            regex_search(line, match, fastaRegex);
            if (!match.empty()) {
                name = match[1];
            }
            sequencestream.str(std::string());
        }else{
            sequencestream << line;
        }
    }
    if( name.size()>0 ){
        std::string sequence = sequencestream.str();
        std::transform(sequence.begin(), sequence.end(), sequence.begin(),::toupper);
        //if( sequence.size()>0 ){
            sequences[name]=sequence;
            seqNames.push_back(name);
//            std::cout << std::endl << "line 98" << name << " " << sequence << std::endl;
            regex_search(meta_line, match, metaRegex);
            if (!match.empty()) {
//                std::cout << "line 72:" << meta_line << std::endl;
                int strand = std::stoi(match[3]);
                int32_t start = std::stoi(match[4]);
                int32_t end = std::stoi(match[5]);
                FastaMeta fastaMeta(match[1], match[2], start, end, strand);
                metaInformations[name]=fastaMeta;
            }
        //}
    }
    infile.close();
}

void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences ){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit (1);
    }
    std::string name;
    std::stringstream sequencestream;
    std::string line="";
    while (std::getline(infile, line)){
        if( line[0] == '>'  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            // this method is faster than regular expression
            if( name.size()>0 ){
                std::string sequence = sequencestream.str();
                (sequence.begin(), sequence.end(), sequence.begin(),::toupper);
                sequences[name]=sequence;
            }
            name=line.substr(1, line.find(" ", 0)-1);
            if( name[0] == '>' ){
                name=name.substr(1, name.find("\t", 0)-1);
            }else{
                name=name.substr(0, name.find("\t", 0)-1);
            }
            sequencestream.str(std::string());
        }else{
            sequencestream << line;
        }
    }
    if( name.size()>0 ){
        std::string sequence = sequencestream.str();
        std::transform(sequence.begin(), sequence.end(), sequence.begin(),::toupper);
        if( sequence.size()>0 ){
            sequences[name]=sequence;
        }
    }
    infile.close();
}


std::string getReverseComplementary(const std::string& sequence) {
    std::string reserveString;
    reserveString.reserve(sequence.size());
    std::stringstream reversecomplementary(reserveString);
    for (int i = sequence.length() - 1; i >= 0; i--) {
        char c = sequence[i];
        if ('A' == c) {
            c = 'T';
        } else if ('T' == c) {
            c = 'A';
        } else if ('U' == c) {
            c = 'A';
        } else if ('C' == c) {
            c = 'G';
        } else if ('G' == c) {
            c = 'C';
        } else if ('R' == c) {
            c = 'Y';
        } else if ('Y' == c) {
            c = 'R';
        } else if ('K' == c) {
            c = 'M';
        } else if ('M' == c) {
            c = 'K';
        } else if ('B' == c) {
            c = 'V';
        } else if ('V' == c) {
            c = 'B';
        } else if ('D' == c) {
            c = 'H';
        } else if ('H' == c) {
            c = 'D';
        }
        reversecomplementary<< c;
    }
    return reversecomplementary.str();
}

std::string getSubSequence( std::string sequence, const int& _start, const int& _end){
    size_t start = _start;
    size_t end = _end;
    if( start > end ){
        size_t temp = start;
        start=end;
        end=temp;
    }
    if( start < 1 ){
        start = 1;
    }
    if( start > sequence.size() ){
        start = sequence.size();
        end = start;
    }else if( end > sequence.size() ){
        end = sequence.size();
    }
    return sequence.substr(start-1, end-start+1);
}

std::string getSubSequence( std::string & sequence, const int& _start, const int& _end, const int& strand){
    if( strand == 1 ){ //1 is positive
        return  getSubSequence( sequence, _start, _end);
    }else{
        std::string seq = getSubSequence( sequence,  _start, _end);
        return getReverseComplementary(seq);
    }
}


FastaMeta::FastaMeta(){

}
FastaMeta::FastaMeta(const std::string &species, const std::string &chr, int32_t & start, int32_t & end, int & strand)
        : species(species), chr(chr), start(start), end(end), strand(strand) {}
const std::string &FastaMeta::getSpecies() const {
    return species;
}

void FastaMeta::setSpecies(const std::string &species) {
    FastaMeta::species = species;
}

const std::string &FastaMeta::getChr() const {
    return chr;
}

void FastaMeta::setChr(const std::string &chr) {
    FastaMeta::chr = chr;
}

int32_t FastaMeta::getStart() const {
    return start;
}

void FastaMeta::setStart(int32_t start) {
    FastaMeta::start = start;
}

int32_t FastaMeta::getEnd() const {
    return end;
}

void FastaMeta::setEnd(int32_t end) {
    FastaMeta::end = end;
}

int FastaMeta::getStrand() const {
    return strand;
}

void FastaMeta::setStrand(int strand) {
    FastaMeta::strand = strand;
}
