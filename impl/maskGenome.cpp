//
// Created by Baoxing Song on 2019-08-07.
//

#include "maskGenome.h"


void maskGenome( std::map<std::string, std::string>& genome, char & alterChar, std::string & kmerProfileFile,
        int32_t & maskKmerFrequency, int32_t & outputLineWidth, std::map<std::string, std::string> & ifCds,
        std::string & outputFile, std::vector<std::string> & seqNames, bool & ifSoftmask ){
    //assert(maskKmerFrequency > 0);



    std::set<std::string> kmers;
    int32_t kmerSize=10;
    if( maskKmerFrequency > 0 ){
        std::ifstream infile(kmerProfileFile);
        if( ! infile.good()){
            std::cerr << "error in opening kmer Profile file " << kmerProfileFile << std::endl;
            exit (1);
        }
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
        kmerSize=sequence.size();
        if( count > maskKmerFrequency ){
            kmers.insert(sequence);
        }
        infile.close();
    }

    std::string maskKmer = std::string(kmerSize, alterChar);
//    std::cout << "maskKmer:" << maskKmer << std::endl;
    std::map<std::string, std::string> ifMasks;
    for( std::map<std::string, std::string>::iterator it = genome.begin(); it!=genome.end(); ++it ){
        std::string ifMask = std::string(it->second.size(), '0');
        for( int32_t i=0; i<= (it->second).size() - kmerSize; ++i ){
            if ( kmers.find(it->second.substr(i, kmerSize) )!= kmers.end() ){
                ifMask.replace(i, kmerSize, kmerSize, '1');
            }
        }
        ifMasks[it->first] = ifMask;
    }

    std::ofstream ofile;
    ofile.open(outputFile);
    std::locale loc;
    for( std::string chr : seqNames ){
        ofile << ">" << chr << std::endl;
        std::string seq = genome[chr];
        std::transform(seq.begin(), seq.end(), seq.begin(),::toupper);
        if (ifSoftmask){
            if( ifCds.find(chr) != ifCds.end() ){
                for( int32_t i=0; i < seq.size(); ++i ){
                    if ( ifMasks[chr][i] == '0' && ifCds[chr][i] == '0' ){
                        ofile << seq[i];
                    }else{
                        ofile << std::tolower(seq[i], loc);
                    }
                    if ( i % outputLineWidth == (outputLineWidth-1) ){
                        ofile << std::endl;
                    }
                }
            }else{
                for( int32_t i=0; i < seq.size(); ++i ){
                    if ( ifMasks[chr][i] == '0' ){
                        ofile << seq[i];
                    }else{
                        ofile << std::tolower(seq[i], loc);
                    }
                    if ( i % outputLineWidth == (outputLineWidth-1) ){
                        ofile << std::endl;
                    }
                }
            }
        }else{
            if( ifCds.find(chr) != ifCds.end() ){
                for( int32_t i=0; i < seq.size(); ++i ){
                    if ( ifMasks[chr][i] == '0' && ifCds[chr][i] == '0' ){
                        ofile << seq[i];
                    }else{
                        ofile << alterChar;
                    }
                    if ( i % outputLineWidth == (outputLineWidth-1) ){
                        ofile << std::endl;
                    }
                }
            }else{
                for( int32_t i=0; i < seq.size(); ++i ){
                    if ( ifMasks[chr][i] == '0' ){
                        ofile << seq[i];
                    }else{
                        ofile << alterChar;
                    }
                    if ( i % outputLineWidth == (outputLineWidth-1) ){
                        ofile << std::endl;
                    }
                }
            }
        }

        ofile << std::endl;
    }
    ofile.close();
}
