
#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>

TEST(kmerreading, c1){
    int32_t  maskKmerFrequency = 100;
    std::string kmerProfileFile = "/media/bs674/2t/testPan_cns/realStory/maskGenomeForGenomeAlignment/sorghum_k20_count_dumps.fa";
    std::ifstream infile(kmerProfileFile);
    if( ! infile.good()){
        std::cerr << "error in opening kmer Profile file " << kmerProfileFile << std::endl;
        exit (1);
    }
    std::set<std::string> kmers;
    int32_t count=0;
    std::stringstream sequencestream;
    std::string line="";
    while (std::getline(infile, line)){
        if( line[0] == '>'  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            // this method is faster than regular expression
            if( count > maskKmerFrequency ){
                std::string sequence = sequencestream.str();
                std::transform(sequence.begin(), sequence.end(), sequence.begin(),::toupper);
                kmers.insert(sequence);
            }
            count = std::stoi(line.substr(1, line.find(" ", 0)-1));
            sequencestream.str(std::string());
        }else{
            sequencestream << line;
        }
    }
    std::string sequence = sequencestream.str();
    std::transform(sequence.begin(), sequence.end(), sequence.begin(),::toupper);
    int32_t kmerSize=sequence.size();
    if( count > maskKmerFrequency ){
        kmers.insert(sequence);

    }
    std::cout << "kmerSize:" << kmerSize << std::endl;

    std::string maskKmer = std::string(kmerSize, 'n');
    std::cout << "maskKmer:" << maskKmer << std::endl;
//    std::map<std::string, std::string> ifMasks;
//    for( std::map<std::string, std::string>::iterator it = genome.begin(); it!=genome.end(); ++it ){
//        std::string ifMask = std::string(it->second.size(), '0');
//        for( int32_t i=0; i<= (it->second).size() - kmerSize; ++i ){
//            if ( kmers.find(it->second.substr(i, kmerSize) )!= kmers.end() ){
//                ifMask.replace(i, kmerSize, kmerSize, '1');
//            }
//        }
//        ifMasks[it->first] = ifMask;
//    }

    ASSERT_EQ(0, 0);
}
