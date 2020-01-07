//
// Created by bs674 on 8/9/19.
//

#include "samFileToSeeds.h"


int32_t numberOfNInAInterval( std::string & seq, int32_t start, int32_t end ){
    int32_t numberOfNs = 0;
    for( int32_t i= start; i<=end; ++i){
        if( seq[i]=='n' ){
            ++numberOfNs;
        }
    }
    return numberOfNs;
}


void samFileToSeeds ( std::string & samFile, std::map<std::string, std::map<std::string, std::vector<Seed>>> & positiveSeeds,
                      std::map<std::string, std::map<std::string, std::vector<Seed>>> & negativeSeeds,
                      std::map<std::string, std::string> & refGenome, std::map<std::string, std::string> & queryGenome){
    std::ifstream infile(samFile);
    if( ! infile.good()){
        std::cerr << "error in opening sam file " << samFile << std::endl;
        exit (1);
    }
    std::regex samRegex("^(.*?)\\s+(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\w+)");
    std::smatch samMatch;
    std::regex cigarRegex ("(\\d+)([MIDNSH])");
    std::smatch cigarMatch;

    std::regex cigarRegex1 ("^(\\d+)([SH])");
    std::smatch cigarMatch1;

    int32_t refStart;
    int32_t refEnd;
    int32_t queryStart;
    int32_t queryEnd;
    int32_t length;
    std::string line;
    std::string cigarString;
    std::string cigarChar;
    std::string chromosomeName;
    int32_t samLineNumber  = 0;
    while (std::getline(infile, line)){
        if( 0 == samLineNumber%1000 ){
            std::cout << "reading samfile file line " << samLineNumber << std::endl;
        }
        ++samLineNumber;
//        std::cout << line << std::endl;
        if(line[0]=='@' || line.size()<9){
        }else {
            regex_search(line, samMatch, samRegex);
            if (!samMatch.empty()) {
                cigarString = samMatch[6];
                regex_search(cigarString, cigarMatch1, cigarRegex1);
                if ( !cigarMatch1.empty() ) {
                    std::string referenceChromosomeName = samMatch[3];
                    std::string queryChromosomeName = samMatch[1];
                    refStart = stoi(samMatch[4]);
                    queryStart = stoi(cigarMatch1[1]);
                    refEnd = refStart;
                    queryEnd = queryStart;
                    while (std::regex_search(cigarString, cigarMatch, cigarRegex)) {
                        cigarChar = cigarMatch[2];
                        length = stoi(cigarMatch[1]);
                        if (cigarChar[0] == 'M'){
                            refEnd += length;
                            queryEnd += length;
                        }else if( cigarChar[0] == 'D') {
                            refEnd += length;
                        }else if (cigarChar[0] == 'N') {
                            refEnd += length;
                        }else if (cigarChar[0] == 'I') {
                            queryEnd += length;
                        }
                        cigarString = cigarMatch.suffix().str();
                    }
                    int32_t cigarFlag = std::stoi(samMatch[2]);

                    refStart = refStart - 1;
                    refEnd = refEnd - 2;
                    queryStart = queryStart - 1;
                    queryEnd = queryEnd - 2;


                    int32_t numberOfNs1 = numberOfNInAInterval( refGenome[referenceChromosomeName], 0, refStart-1);
                    int32_t numberOfNs2 = numberOfNInAInterval( refGenome[referenceChromosomeName], refStart, refEnd);
                    refStart -= numberOfNs1;
                    refEnd -= numberOfNs1;
                    refEnd -= numberOfNs2;

                    if ( 0 == cigarFlag ){
                        if ( positiveSeeds.find(referenceChromosomeName) == positiveSeeds.end() ){
                            std::map<std::string, std::vector<Seed>> temp;
                            positiveSeeds[referenceChromosomeName] = temp;
                            std::vector<Seed> temp2;
                            positiveSeeds[referenceChromosomeName][queryChromosomeName] = temp2;
                        }else if( positiveSeeds[referenceChromosomeName].find(queryChromosomeName) == positiveSeeds[referenceChromosomeName].end() ){
                            std::vector<Seed> temp2;
                            positiveSeeds[referenceChromosomeName][queryChromosomeName] = temp2;
                        }

                        numberOfNs1 = numberOfNInAInterval( queryGenome[queryChromosomeName], 0, queryStart-1);
                        numberOfNs2 = numberOfNInAInterval( queryGenome[queryChromosomeName], queryStart, queryEnd);
                        queryStart -= numberOfNs1;
                        queryEnd -= numberOfNs1;
                        queryEnd -= numberOfNs2;

                        Seed seed(refStart, refEnd, queryStart, queryEnd);
                        positiveSeeds[referenceChromosomeName][queryChromosomeName].push_back(seed);
//                        std::cerr << line << std::endl;
//                        std::cerr << "x_seed.start1:" << seed.getStart1() << "x_seed.end1:" << seed.getEnd1() << "x_seed.start2:" << seed.getStart2() << "x_seed.end2:" << seed.getEnd1() << std::endl;
                    }else{

                        numberOfNs1 = numberOfNInAInterval( queryGenome[queryChromosomeName], queryStart, queryEnd);
                        numberOfNs2 = numberOfNInAInterval( queryGenome[queryChromosomeName], queryEnd+1, queryGenome[queryChromosomeName].length()-1 );

                        queryEnd = queryGenome[queryChromosomeName].length()-1-queryStart;
                        queryStart = queryGenome[queryChromosomeName].length()-1-queryEnd;

                        queryStart -= numberOfNs2;
                        queryEnd -= numberOfNs2;
                        queryEnd -= numberOfNs1;

                        if ( negativeSeeds.find(referenceChromosomeName) == negativeSeeds.end() ){
                            std::map<std::string, std::vector<Seed>> temp;
                            negativeSeeds[referenceChromosomeName] = temp;
                            std::vector<Seed> temp2;
                            negativeSeeds[referenceChromosomeName][queryChromosomeName] = temp2;
                        }else if( negativeSeeds[referenceChromosomeName].find(queryChromosomeName) == negativeSeeds[referenceChromosomeName].end() ){
                            std::vector<Seed> temp2;
                            negativeSeeds[referenceChromosomeName][queryChromosomeName] = temp2;
                        }
                        Seed seed(refStart, refEnd, queryStart, queryEnd);
                        negativeSeeds[referenceChromosomeName][queryChromosomeName].push_back(seed);
                    }
                }
            }
        }
    }
}
