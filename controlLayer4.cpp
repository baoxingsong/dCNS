/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 09:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************




 ************************************************************************/

#include "controlLayer.h"


std::string softwareName = "dCNS";
int32_t seed_window_size = 38;
int32_t mini_cns_score = 40;
int32_t step_size = 8;
int32_t matrix_boundary_distance = 0;

int32_t matchingScore = 2;
int32_t mismatchingPenalty = -3;
int32_t openGapPenalty1 = -4;
int32_t extendGapPenalty1 = -2;

int32_t openGapPenalty2 = -45;
int32_t extendGapPenalty2 = 0;

int32_t zDrop = 50; // this one should be slightly larger than openGapPenalty2
int32_t bandwidth = 20000; // using very large value

double lambda = 0.382291; // for 100kb *100kb sequence alignment, score 53-54 gives a p-value around 0.1
double kValue = 0.006662;

bool onlySyntenic = false;
double pvalues = 0.1;

int32_t w = 10;  //this is the band width for band sequence alignments
int32_t xDrop = 20;



void outputNoMasking(std::string _output, std::vector<PairedSimilarFragment> pairedSimilarFragments0, int32_t refStrand,
        int32_t  queStrand, int32_t  refStart, int32_t  queStart, std::string refChr, std::string queChr,
                     std::string & seq1_string, std::string & seq2_string){
    std::ofstream ofile;
    ofile.open(_output+".sam", std::ios_base::app);
    std::ofstream ofile2;
    ofile2.open(_output+"o", std::ios_base::app);
//    std::cout << "line 63, size:" << pairedSimilarFragments0.size() <<"refStrand:" << refStrand << "queStrand:" << queStrand<< std::endl;
//    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    ofile2 << "#" << refChr << ":" << refStrand << ":" << refStart << "-" << queChr << ":" << queStrand << ":" << queStart << std::endl;
    for( int32_t i=0; i<pairedSimilarFragments0.size(); ++i ){
        if( refStrand == 1 && 1==queStrand ){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int32_t j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << "\t*" << std::endl;
        }else if ( refStrand == 0 && 1==queStrand){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + seq1_string.size() -1 - pairedSimilarFragments0[i].getEnd1()+1;
            ofile << queChr << "\t16\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int32_t j=pairedSimilarFragments0[i].getCigar().size()-1; j>=0; --j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2(), 0) << "\t*" << std::endl;
        }else if( refStrand == 1 && 0==queStrand ){
//            std::cout << "line 86" << std::endl;
            int32_t quePos = queStart + seq2_string.size()-1-pairedSimilarFragments0[i].getEnd2()+1; // this is the position on the positive strand
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            ofile << queChr << "\t16\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int32_t j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << "\t*" << std::endl;
        }else  if( refStrand == 0 && refStrand==queStrand ){
            int32_t quePos = queStart + seq2_string.size()-1-pairedSimilarFragments0[i].getEnd2()+1;
            int32_t refPos = refStart + seq1_string.size() -1 - pairedSimilarFragments0[i].getEnd1()+1;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int32_t j=pairedSimilarFragments0[i].getCigar().size()-1; j>=0; --j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2(), 0) << "\t*" << std::endl;
        }

        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
               << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

    }
    ofile.close();
    ofile2 << "#done" << std::endl;
    ofile2.close();
}



void outputNoMasking(std::string _output, PairedSimilarFragment pairedSimilarFragment, int32_t refStrand,
                     int32_t  queStrand, int32_t  refStart, int32_t  queStart, std::string refChr, std::string queChr,
                     std::string & seq1_string, std::string & seq2_string){
    std::ofstream ofile;
    ofile.open(_output+".sam", std::ios_base::app);
    std::ofstream ofile2;
    ofile2.open(_output+"o", std::ios_base::app);

    ofile2 << "#" << refChr << ":" << refStrand << ":" << refStart << "-" << queChr << ":" << queStrand << ":" << queStart << std::endl;
    if( refStrand == 1 && 1==queStrand ){
        int32_t quePos = queStart + pairedSimilarFragment.getStart2() - 1;
        int32_t refPos = refStart + pairedSimilarFragment.getStart1() - 1;
        ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragment.getScore() << "\t" << quePos << "H";
        for( int32_t j=0; j<pairedSimilarFragment.getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragment.getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragment.getCigar()[j]&0xf;
            ofile << cigarLength << "MID"[cigarType];
        }
        ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragment.getStart2(), pairedSimilarFragment.getEnd2()) << "\t*" << std::endl;
    }else if ( refStrand == 0 && 1==queStrand){
        int32_t quePos = queStart + pairedSimilarFragment.getStart2() - 1;
        int32_t refPos = refStart + seq1_string.size() -1 - pairedSimilarFragment.getEnd1()+1;
        ofile << queChr << "\t16\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragment.getScore() << "\t" << quePos << "H";
        for( int32_t j=pairedSimilarFragment.getCigar().size()-1; j>=0; --j ){
            uint32_t cigarLength = pairedSimilarFragment.getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragment.getCigar()[j]&0xf;
            ofile << cigarLength << "MID"[cigarType];
        }
        ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragment.getStart2(), pairedSimilarFragment.getEnd2(), 0) << "\t*" << std::endl;
    }else if( refStrand == 1 && 0==queStrand ){
//            std::cout << "line 86" << std::endl;
        int32_t quePos = queStart + seq2_string.size()-1-pairedSimilarFragment.getEnd2()+1; // this is the position on the positive strand
        int32_t refPos = refStart + pairedSimilarFragment.getStart1() - 1;
        ofile << queChr << "\t16\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragment.getScore() << "\t" << quePos << "H";
        for( int32_t j=0; j<pairedSimilarFragment.getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragment.getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragment.getCigar()[j]&0xf;
            ofile << cigarLength << "MID"[cigarType];
        }
        ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragment.getStart2(), pairedSimilarFragment.getEnd2()) << "\t*" << std::endl;
    }else  if( refStrand == 0 && refStrand==queStrand ){
        int32_t quePos = queStart + seq2_string.size()-1-pairedSimilarFragment.getEnd2()+1;
        int32_t refPos = refStart + seq1_string.size() -1 - pairedSimilarFragment.getEnd1()+1;
        ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragment.getScore() << "\t" << quePos << "H";
        for( int32_t j=pairedSimilarFragment.getCigar().size()-1; j>=0; --j ){
            uint32_t cigarLength = pairedSimilarFragment.getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragment.getCigar()[j]&0xf;
            ofile << cigarLength << "MID"[cigarType];
        }
        ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragment.getStart2(), pairedSimilarFragment.getEnd2(), 0) << "\t*" << std::endl;
    }

    ofile2 << pairedSimilarFragment.getStart1() << " " << pairedSimilarFragment.getEnd1()
           << " " << pairedSimilarFragment.getStart2() << " " << pairedSimilarFragment.getEnd2() << " " << pairedSimilarFragment.getScore() << " " << pairedSimilarFragment.getPValue() << std::endl;


    ofile.close();
    ofile2 << "#done" << std::endl;
    ofile2.close();
}


void pairCnsXExtend(std::string & _input,  std::string & _reference, std::string & _output,
        int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty,
        int32_t & _extendGapPenalty, int32_t & _seed_window_size, int32_t & _mini_cns_score,
        int32_t & _step_size, int32_t & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
        const double & _kValue, const int32_t & _w, const int32_t & _xDrop, double & _pvalues){

    Scorei m(_matchingScore, _mismatchingPenalty);

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    std::map<std::string, FastaMeta> metaInformations;
    readFastaFile( _input, sequences, metaInformations, seqNames );

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_string = sequences[_reference];
    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int32_t length1 = seq1_string.length();
    //sequences.clear();
    std::cout << "reference:" << seq1_string << std::endl;
    for( std::string _query : seqNames ){
        if( _query.compare(_reference) != 0) {
            std::string seq2_string = sequences[_query];
            //seq2_string = getReverseComplementary(seq2_string);
            int8_t *seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
            int8_t *seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
            int32_t length2 = seq2_string.length();
            std::cout << "_query:" << _query << ": " << seq2_string << std::endl;
            std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
                    findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                          seq2, seq2_rev_com, length1, length2, _seed_window_size,
                                                          _mini_cns_score, _matrix_boundary_distance, _openGapPenalty,
                                                          _extendGapPenalty, _matchingScore,
                                                          _mismatchingPenalty, m, _step_size, seq1_string, seq2_string,
                                                          _pvalues, _lambda, _kValue, _w, _xDrop);

            if (_onlySyntenic) {
                std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                pairedSimilarFragments0 = pairedSimilarFragments;
            }

            outputNoMasking(_output+metaInformations[_query].getSpecies(), pairedSimilarFragments0,
                    metaInformations[_reference].getStrand(), metaInformations[_query].getStrand(),
                    metaInformations[_reference].getStart(), metaInformations[_query].getStart(),
                    metaInformations[_reference].getChr(), metaInformations[_query].getChr(),
                    seq1_string, seq2_string);
            delete seq2;
            delete seq2_rev_com;
        }
    }
    delete seq1;
    delete seq1_rev_com;
}


int pairCnsXExtend(int argc, char** argv){

    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" pairCnsXExtend -i input -r reference -q query -o output" << std::endl<<
        "Options" << std::endl <<
        " -i FILE    input file in fasta format" << std::endl <<
        " -r STRING  reference sequence entry name in input fasta file" << std::endl <<
//        " -q STRING  query sequence entry name in input fasta file" << std::endl <<
        " -o FILE    output file, exists file will be overwrote" << std::endl <<
        " -y bool    only output syntenic result (default: false)" << std::endl <<
        " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
        " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
        " -O INT     open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
        " -E INT     extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
        " -u INT     xextend alignment band width (default: " << w << ")" << std::endl <<
        " -w INT     the windows size used to run the smith-waterman algorithm to get the alignment seed (default: "<<seed_window_size<<")" << std::endl <<
        //" -m INT     minimum seeds size to trigger a alignment extension and the minimum conserved sequence to report (default: " << mini_cns_seed_size << ")" << std::endl <<
        " -c INT     minimum seeds score to trigger a alignment extension (default: " << mini_cns_score << ")" << std::endl <<
        " -s INT     step size for sliding the smith-waterman seeds alignment window (default: " << step_size << ")" << std::endl <<
        " -p DOUBLE  pvalue for significant alignment output (default: " << pvalues << ")" << std::endl <<
        " -b INT     minimum number of base-pairs between the maximum smith-waterman score with the score matrix boundary (default: " << matrix_boundary_distance << ") , it is not important" << std::endl;
    InputParser inputParser (argc, argv);

    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O") );
    }
    if( inputParser.cmdOptionExists("-E") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E") );
    }
    if( inputParser.cmdOptionExists("-u") ){
        w = std::stoi( inputParser.getCmdOption("-u") );
    }
    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }
    if( inputParser.cmdOptionExists("-w") ){
        seed_window_size = std::stoi( inputParser.getCmdOption("-w") );
    }
    if( inputParser.cmdOptionExists("-c") ){
        mini_cns_score = std::stoi( inputParser.getCmdOption("-c") );
    }
    if( inputParser.cmdOptionExists("-s") ){
        step_size = std::stoi( inputParser.getCmdOption("-s") );
    }
    if( inputParser.cmdOptionExists("-b") ){
        matrix_boundary_distance = std::stoi( inputParser.getCmdOption("-b") );
    }
    if( inputParser.cmdOptionExists("-p") ){
        pvalues = std::stod( inputParser.getCmdOption("-p") );
    }

    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && /*inputParser.cmdOptionExists("-q") &&*/ inputParser.cmdOptionExists("-o")  ){
        std::string input = inputParser.getCmdOption("-i");
        std::string reference = inputParser.getCmdOption("-r");
        //std::string query = inputParser.getCmdOption("-q");
        std::string output = inputParser.getCmdOption("-o");

        pairCnsXExtend(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                      extendGapPenalty1, seed_window_size, mini_cns_score,
                      step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, xDrop, pvalues);

        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}



void pairCns2Gaps(std::string & _input,  std::string & _reference, std::string & _output,
                    int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty1,
                    int32_t & _extendGapPenalty1, int32_t & _openGapPenalty2,
                    int32_t & _extendGapPenalty2, int32_t & _seed_window_size, int32_t & _mini_cns_score,
                    int32_t & _step_size, int32_t & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
                    const double & _kValue, const int32_t & _w, const int32_t & _bandwidth, const int32_t & _xDrop,
                    const int32_t & _zDrop, double & _pvalues){

    std::map<std::string, FastaMeta> metaInformations;

    Scorei m(_matchingScore, _mismatchingPenalty);


    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( _input, sequences, metaInformations, seqNames );

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_string = sequences[_reference];
    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int32_t length1 = seq1_string.length();
    sequences.clear();
    for( std::string _query : seqNames ){
        if( _query.compare(_reference) != 0) {
            std::string seq2_string = sequences[_query];
            int8_t *seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
            int8_t *seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
            int32_t length2 = seq2_string.length();

            std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
                    findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                           _seed_window_size, _mini_cns_score,
                                                           _matrix_boundary_distance, _openGapPenalty1, _extendGapPenalty1,
                                                           _openGapPenalty2, _extendGapPenalty2, _matchingScore,
                                                           _mismatchingPenalty, m, _step_size, seq1_string, seq2_string, _pvalues,
                                                           _lambda, _kValue, _zDrop, _bandwidth, _w, _xDrop);

            if (_onlySyntenic) {
                std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                pairedSimilarFragments0 = pairedSimilarFragments;
            }

            outputNoMasking(_output+metaInformations[_query].getSpecies(), pairedSimilarFragments0, metaInformations[_reference].getStrand(),
                            metaInformations[_query].getStrand(), metaInformations[_reference].getStart(), metaInformations[_query].getStart(), metaInformations[_reference].getChr(), metaInformations[_query].getChr(),
                            seq1_string, seq2_string);
            delete seq2;
            delete seq2_rev_com;
        }
    }
    delete seq1;
    delete seq1_rev_com;
}

int pairCns2Gaps(int argc, char** argv){
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" pairCns2Gaps -i input -r reference -q query -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    input file in fasta format" << std::endl <<
          " -r STRING  reference sequence entry name in input fasta file" << std::endl <<
          //        " -q STRING  query sequence entry name in input fasta file" << std::endl <<
          " -o FILE    output file, exists file will be overwrote" << std::endl <<
          " -y bool    only output syntenic result (default: false)" << std::endl <<
          " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT    open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT    extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -u INT     xextend alignment band width (default: " << w << ")" << std::endl <<
          " -O2 INT    open gap penalty (default: " << openGapPenalty2 << ")" << std::endl <<
          " -E2 INT    extend gap penalty (default: " << extendGapPenalty2 << ")" << std::endl <<
          " -u2 INT    z drop alignment band width (default: " << bandwidth << ")" << std::endl <<
          //" - INT     xextend alignment band width (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -w INT     the windows size used to run the smith-waterman algorithm to get the alignment seed (default: "<<seed_window_size<<")" << std::endl <<
          //" -m INT     minimum seeds size to trigger a alignment extension and the minimum conserved sequence to report (default: " << mini_cns_seed_size << ")" << std::endl <<
          " -c INT     minimum seeds score to trigger a alignment extension (default: " << mini_cns_score << ")" << std::endl <<
          " -s INT     step size for sliding the smith-waterman seeds alignment window (default: " << step_size << ")" << std::endl <<
          " -p DOUBLE  pvalue for significant alignment output (default: " << pvalues << ")" << std::endl <<
          " -b INT     minimum number of base-pairs between the maximum smith-waterman score with the score matrix boundary (default: " << matrix_boundary_distance << ") , it is not important" << std::endl;
    InputParser inputParser (argc, argv);


    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O1") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
    }
    if( inputParser.cmdOptionExists("-E1") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
    }
    if( inputParser.cmdOptionExists("-O2") ){
        openGapPenalty2 = std::stoi( inputParser.getCmdOption("-O2") );
    }
    if( inputParser.cmdOptionExists("-E2") ){
        extendGapPenalty2 = std::stoi( inputParser.getCmdOption("-E2") );
    }
    if( inputParser.cmdOptionExists("-u") ){
        w = std::stoi( inputParser.getCmdOption("-u") );
    }
    if( inputParser.cmdOptionExists("-u2") ){
        bandwidth = std::stoi( inputParser.getCmdOption("-u2") );
    }
    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }
    if( inputParser.cmdOptionExists("-w") ){
        seed_window_size = std::stoi( inputParser.getCmdOption("-w") );
    }
//    if( inputParser.cmdOptionExists("-m") ){
//        mini_cns_seed_size = std::stoi( inputParser.getCmdOption("-m") );
//    }
    if( inputParser.cmdOptionExists("-c") ){
        mini_cns_score = std::stoi( inputParser.getCmdOption("-c") );
    }
    if( inputParser.cmdOptionExists("-s") ){
        step_size = std::stoi( inputParser.getCmdOption("-s") );
    }
    if( inputParser.cmdOptionExists("-b") ){
        matrix_boundary_distance = std::stoi( inputParser.getCmdOption("-b") );
    }
    if( inputParser.cmdOptionExists("-p") ){
        pvalues = std::stod( inputParser.getCmdOption("-p") );
    }


    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && /*inputParser.cmdOptionExists("-q") &&*/ inputParser.cmdOptionExists("-o")  ){
        std::string input = inputParser.getCmdOption("-i");
        std::string reference = inputParser.getCmdOption("-r");
        //std::string query = inputParser.getCmdOption("-q");
        std::string output = inputParser.getCmdOption("-o");

        pairCns2Gaps(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                       extendGapPenalty1, openGapPenalty2, extendGapPenalty2, seed_window_size, mini_cns_score,
                       step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, bandwidth, xDrop, zDrop, pvalues);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}




void weighted1Gap(std::string & _input,  std::string & _reference, std::string & _output,
                  int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty1,
                  int32_t & _extendGapPenalty1,  int32_t & _seed_window_size, int32_t & _mini_cns_score,
                  int32_t & _step_size, int32_t & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
                  const double & _kValue, int32_t & _bandwidth, int32_t & _w, const int32_t & _xDrop, const int32_t & _zDrop, const std::string & scoreFoler,
                  const std::string & refFasta, const std::string & gffFile, double & _pvalues){

    std::map<std::string, FastaMeta> metaInformations;

    Scorei m(_matchingScore, _mismatchingPenalty);

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( _input, sequences, metaInformations, seqNames );

    Score score(scoreFoler);
    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    //readFastaFile( "/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa", sequences0, seqNames0);
    readFastaFile( refFasta, sequences0, seqNames0);
    std::map<std::string, int32_t>  chrSize;
    for( std::map<std::string, std::string>::iterator it = sequences0.begin(); it!=sequences0.end(); ++it ){
        chrSize[it->first] = it->second.size();
    }
    std::map<std::string, int16_t *> weight;
    std::map<std::string, int16_t *> weight_rev;
    //readGffFileWithEveryThing("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv3.31.gff3", chrSize, weight, weight_rev);
    readGffFileWithEveryThing(gffFile, chrSize, weight, weight_rev);

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string refChr = metaInformations[_reference].getChr();
    std::string seq1_string = sequences[_reference];
    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int32_t length1 = seq1_string.length();
    sequences0.clear();
    for( std::string _query : seqNames ){
        if( _query.compare(_reference) != 0) {
            std::string seq2_string = sequences[_query];
            int8_t *seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
            int8_t *seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
            int32_t length2 = seq2_string.length();

            std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
                    findSimilarFragmentsForPairedSequence_wighted_1gap ( seq1, seq1_rev_com,
                                     seq2, seq2_rev_com, length1, length2, _seed_window_size,
                                     _mini_cns_score, _matrix_boundary_distance,
                                     _openGapPenalty1, _extendGapPenalty1,  _matchingScore, _mismatchingPenalty,
                                     m, _step_size, seq1_string, seq2_string,
                                     _pvalues, _lambda, _kValue, _zDrop,
                                     _bandwidth, _w, _xDrop, score,
                                     weight[refChr]+metaInformations[_reference].getStart()-1,
                                     weight_rev[refChr] + (chrSize[refChr]-(metaInformations[_reference].getStart()+length1-1)));



            if (_onlySyntenic) {
                std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                pairedSimilarFragments0 = pairedSimilarFragments;
            }

            outputNoMasking(_output+metaInformations[_query].getSpecies(), pairedSimilarFragments0, metaInformations[_reference].getStrand(),
                            metaInformations[_query].getStrand(), metaInformations[_reference].getStart(),
                            metaInformations[_query].getStart(), metaInformations[_reference].getChr(), metaInformations[_query].getChr(),
                            seq1_string, seq2_string);
            delete seq2;
            delete seq2_rev_com;
        }
    }
    delete seq1;
    delete seq1_rev_com;
}

int weighted1Gap(int argc, char** argv){
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" weighted1Gap -i input -r reference -q query -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    input file in fasta format" << std::endl <<
          " -r STRING  reference sequence entry name in input fasta file" << std::endl <<
          " -f STRING  gff3 format reference genome annotation file" << std::endl <<
          " -ra STRING masked reference genome sequence in fasta format" << std::endl <<
          " -sp STRING the folder where score parameters located" << std::endl <<
          //        " -q STRING  query sequence entry name in input fasta file" << std::endl <<
          " -o FILE    output file, exists file will be overwrote" << std::endl <<
          " -y bool    only output syntenic result (default: false)" << std::endl <<
          " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT     open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT     extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
//          " -O2 INT     open gap penalty (default: " << openGapPenalty2 << ")" << std::endl <<
//          " -E2 INT     extend gap penalty (default: " << extendGapPenalty2 << ")" << std::endl <<
          " -u INT     xextend alignment band width (default: " << w << ")" << std::endl <<
          //" - INT     xextend alignment band width (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -w INT     the windows size used to run the smith-waterman algorithm to get the alignment seed (default: "<<seed_window_size<<")" << std::endl <<
          //" -m INT     minimum seeds size to trigger a alignment extension and the minimum conserved sequence to report (default: " << mini_cns_seed_size << ")" << std::endl <<
          " -c INT     minimum seeds score to trigger a alignment extension (default: " << mini_cns_score << ")" << std::endl <<
          " -s INT     step size for sliding the smith-waterman seeds alignment window (default: " << step_size << ")" << std::endl <<
          " -p DOUBLE  pvalue for significant alignment output (default: " << pvalues << ")" << std::endl <<
          " -b INT     minimum number of base-pairs between the maximum smith-waterman score with the score matrix boundary (default: " << matrix_boundary_distance << ") , it is not important" << std::endl;
    InputParser inputParser (argc, argv);


    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O1") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
    }
    if( inputParser.cmdOptionExists("-E1") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
    }
    if( inputParser.cmdOptionExists("-u") ){
        w = std::stoi( inputParser.getCmdOption("-u") );
    }
    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }
    if( inputParser.cmdOptionExists("-w") ){
        seed_window_size = std::stoi( inputParser.getCmdOption("-w") );
    }
//    if( inputParser.cmdOptionExists("-m") ){
//        mini_cns_seed_size = std::stoi( inputParser.getCmdOption("-m") );
//    }
    if( inputParser.cmdOptionExists("-c") ){
        mini_cns_score = std::stoi( inputParser.getCmdOption("-c") );
    }
    if( inputParser.cmdOptionExists("-s") ){
        step_size = std::stoi( inputParser.getCmdOption("-s") );
    }
    if( inputParser.cmdOptionExists("-b") ){
        matrix_boundary_distance = std::stoi( inputParser.getCmdOption("-b") );
    }
    if( inputParser.cmdOptionExists("-p") ){
        pvalues = std::stod( inputParser.getCmdOption("-p") );
    }


    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && /*inputParser.cmdOptionExists("-q") &&*/
        inputParser.cmdOptionExists("-o") && inputParser.cmdOptionExists("-ra") && inputParser.cmdOptionExists("-sp")
                                                                                   && inputParser.cmdOptionExists("-f")){
        std::string input = inputParser.getCmdOption("-i");
        std::string reference = inputParser.getCmdOption("-r");
        //std::string query = inputParser.getCmdOption("-q");
        std::string output = inputParser.getCmdOption("-o");
        std::string refFasta = inputParser.getCmdOption("-ra");
        std::string scoreFoler = inputParser.getCmdOption("-sp");
        std::string gffFile = inputParser.getCmdOption("-f");

        weighted1Gap(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                     extendGapPenalty1,seed_window_size, mini_cns_score,
                     step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, bandwidth, w, xDrop, zDrop, scoreFoler,
                     refFasta, gffFile, pvalues);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}



void weighted2Gaps(std::string & _input,  std::string & _reference, std::string & _output,
                  int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty1,
                  int32_t & _extendGapPenalty1, int32_t & _seed_window_size, int32_t & _mini_cns_score,
                  int32_t & _step_size, int32_t & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
                  const double & _kValue, int32_t & _w, const int32_t & _bandwidth, const int32_t & _xDrop, const int32_t & _zDrop, const std::string & scoreFoler,
                  const std::string & refFasta,  const std::string & gffFile, double & _pvalues){

    std::map<std::string, FastaMeta> metaInformations;

    Scorei m(_matchingScore, _mismatchingPenalty);


    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( _input, sequences, metaInformations, seqNames );

    Score score(scoreFoler);
    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    //readFastaFile( "/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa", sequences0, seqNames0);
    readFastaFile( refFasta, sequences0, seqNames0);
    std::map<std::string, int32_t>  chrSize;
    for( std::map<std::string, std::string>::iterator it = sequences0.begin(); it!=sequences0.end(); ++it ){
        chrSize[it->first] = it->second.size();
    }
    std::map<std::string, int16_t *> weight;
    std::map<std::string, int16_t *> weight_rev;
    //readGffFileWithEveryThing("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv3.31.gff3", chrSize, weight, weight_rev);
    readGffFileWithEveryThing(gffFile, chrSize, weight, weight_rev);


    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string refChr = metaInformations[_reference].getChr();
    std::string seq1_string = sequences[_reference];
    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int32_t length1 = seq1_string.length();
    sequences0.clear();
    for( std::string _query : seqNames ){
        if( _query.compare(_reference) != 0) {
            std::string seq2_string = sequences[_query];
            int8_t *seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
            int8_t *seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
            int32_t length2 = seq2_string.length();

            std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
                    findSimilarFragmentsForPairedSequence_wighted ( seq1, seq1_rev_com,
                            seq2, seq2_rev_com, length1, length2, _seed_window_size,
                            _mini_cns_score, _matrix_boundary_distance,
                            _openGapPenalty1, _extendGapPenalty1, _matchingScore, _mismatchingPenalty,
                            m, _step_size, seq1_string, seq2_string,
                            _pvalues, _lambda, _kValue, _zDrop,
                            _bandwidth, _w, _xDrop, score,
                            weight[refChr]+metaInformations[_reference].getStart()-1,
                            weight_rev[refChr] + ((chrSize[refChr]-(metaInformations[_reference].getStart()+length1-1))) );


            if (_onlySyntenic) {
                std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                pairedSimilarFragments0 = pairedSimilarFragments;
            }

            outputNoMasking(_output+metaInformations[_query].getSpecies(), pairedSimilarFragments0, metaInformations[_reference].getStrand(),
                            metaInformations[_query].getStrand(), metaInformations[_reference].getStart(),
                            metaInformations[_query].getStart(), metaInformations[_reference].getChr(), metaInformations[_query].getChr(),
                            seq1_string, seq2_string);
            delete seq2;
            delete seq2_rev_com;
        }
    }
    delete seq1;
    delete seq1_rev_com;
}

int weighted2Gaps(int argc, char** argv){
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" weighted2Gaps -i input -r reference -q query -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    input file in fasta format" << std::endl <<
          " -r STRING  reference sequence entry name in input fasta file" << std::endl <<
          " -f STRING  gff3 format reference genome annotation file" << std::endl <<
          " -ra STRING masked reference genome sequence in fasta format" << std::endl <<
          " -sp STRING the folder where score parameters located" << std::endl <<
          //        " -q STRING  query sequence entry name in input fasta file" << std::endl <<
          " -o FILE    output file, exists file will be overwrote" << std::endl <<
          " -y bool    only output syntenic result (default: false)" << std::endl <<
          " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT    open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT    extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -u INT     xextend alignment band width (default: " << w << ")" << std::endl <<
          " -u2 INT    z drop alignment band width (default: " << bandwidth << ")" << std::endl <<
          //" - INT     xextend alignment band width (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -w INT     the windows size used to run the smith-waterman algorithm to get the alignment seed (default: "<<seed_window_size<<")" << std::endl <<
          //" -m INT     minimum seeds size to trigger a alignment extension and the minimum conserved sequence to report (default: " << mini_cns_seed_size << ")" << std::endl <<
          " -c INT     minimum seeds score to trigger a alignment extension (default: " << mini_cns_score << ")" << std::endl <<
          " -s INT     step size for sliding the smith-waterman seeds alignment window (default: " << step_size << ")" << std::endl <<
          " -p DOUBLE  pvalue for significant alignment output (default: " << pvalues << ")" << std::endl <<
          " -b INT     minimum number of base-pairs between the maximum smith-waterman score with the score matrix boundary (default: " << matrix_boundary_distance << ") , it is not important" << std::endl;
    InputParser inputParser (argc, argv);


    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O1") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
    }
    if( inputParser.cmdOptionExists("-E1") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
    }
    if( inputParser.cmdOptionExists("-u") ){
        w = std::stoi( inputParser.getCmdOption("-u") );
    }
    if( inputParser.cmdOptionExists("-u2") ){
        bandwidth = std::stoi( inputParser.getCmdOption("-u2") );
    }
    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }
    if( inputParser.cmdOptionExists("-w") ){
        seed_window_size = std::stoi( inputParser.getCmdOption("-w") );
    }
//    if( inputParser.cmdOptionExists("-m") ){
//        mini_cns_seed_size = std::stoi( inputParser.getCmdOption("-m") );
//    }
    if( inputParser.cmdOptionExists("-c") ){
        mini_cns_score = std::stoi( inputParser.getCmdOption("-c") );
    }
    if( inputParser.cmdOptionExists("-s") ){
        step_size = std::stoi( inputParser.getCmdOption("-s") );
    }
    if( inputParser.cmdOptionExists("-b") ){
        matrix_boundary_distance = std::stoi( inputParser.getCmdOption("-b") );
    }
    if( inputParser.cmdOptionExists("-p") ){
        pvalues = std::stod( inputParser.getCmdOption("-p") );
    }


    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r")
            && /*inputParser.cmdOptionExists("-q") &&*/ inputParser.cmdOptionExists("-o")
            && inputParser.cmdOptionExists("-ra") && inputParser.cmdOptionExists("-sp")
                                                     && inputParser.cmdOptionExists("-f")){
        std::string input = inputParser.getCmdOption("-i");
        std::string reference = inputParser.getCmdOption("-r");
        //std::string query = inputParser.getCmdOption("-q");
        std::string output = inputParser.getCmdOption("-o");
        std::string scoreFoler = inputParser.getCmdOption("-sp");
        std::string refFasta = inputParser.getCmdOption("-ra");
        std::string gffFile = inputParser.getCmdOption("-f");


        weighted2Gaps(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                     extendGapPenalty1,  seed_window_size, mini_cns_score,
                     step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, bandwidth, xDrop, zDrop, scoreFoler,
                      refFasta, gffFile, pvalues);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}




void outputMasking(std::string _output, std::vector<PairedSimilarFragment> & pairedSimilarFragments0, int32_t refStrand,
                     int32_t queStrand, int32_t  refStart, int32_t queStart, std::string refChr, std::string queChr,
                     std::string & seq1_string, std::string & seq2_string, std::string & seq1_0_string, std::string & seq2_0_string){
    std::ofstream ofile;
    ofile.open(_output+".sam", std::ios_base::app);
    std::ofstream ofile2;
    ofile2.open(_output+"o", std::ios_base::app);
//    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    ofile2 << "#" << refChr << ":" << refStrand << ":" << refStart << "-" << queChr << ":" << queStrand << ":" << queStart << std::endl;
    for( int32_t i=0; i<pairedSimilarFragments0.size(); ++i ){

        if( refStrand==1 && queStrand ==1){

            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;

            int32_t start = 1;
            int32_t end = 1;
            int32_t number_of_chr=0;
            int32_t number_of_n=0;
            int32_t started=0;

            for ( int32_t j=0; j<seq2_0_string.size(); ++j ){
                if ( seq2_0_string[j] == 'n' ){
                    ++number_of_n;
                }else{
                    ++number_of_chr;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart2() && started==0){
                    start += j;
                    quePos += number_of_n;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd2() ){
                    end += j;
                    break;
                }
            }
            number_of_chr=0;
            number_of_n=0;
            int32_t ref_start = 1;
            started=0;
            for ( int32_t j=0; j<seq1_0_string.size(); ++j ){
                if ( seq1_0_string[j] != 'n' ){
                    ++number_of_chr;
                }else{
                    ++number_of_n;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart1() && started==0){
                    ref_start += j;
                    refPos += number_of_n;
                    break;
                }
            }
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";

            int32_t que_i = start-1; // the index of reference
            int32_t ref_i = ref_start-1; // the index of query

            std::vector<uint32_t> newCigar;
            for( int32_t j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                uint32_t op = 0;
                uint32_t length = 1;
                for ( int32_t h=0; h<cigarLength; ++h ){
                    op = cigarType;
                    if( op==0 ){
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }

                        ++ref_i;
                        ++que_i;
                    }else if(op==1) {
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else{
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                    }
                    op = cigarType;
                    if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                        newCigar.push_back(length << 4 | op);
                    }else{
                        newCigar[newCigar.size()-1] += length<<4;
                    }
                }
            }

            for( int32_t j=0; j<newCigar.size(); ++j ){
                uint32_t cigarLength = newCigar[j]>>4;
                uint32_t cigarType = newCigar[j]&0xf;
                ofile << cigarLength << "MIDDSHI=XB"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_0_string, start, end) << "\t*" << std::endl;


        }else if ( refStrand==1 && queStrand ==0) { // this one should be correct

            int32_t quePos = queStart + seq2_0_string.size()-1-pairedSimilarFragments0[i].getEnd2()+1; // this is the position on the positive strand
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;




//            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
//            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;


            int32_t start = 1;
            int32_t end = 1;
            int32_t number_of_chr=0;
            int32_t number_of_n=0;
            int32_t started=0;

            for ( int32_t j=0; j<seq2_0_string.size(); ++j ){
                if ( seq2_0_string[j] == 'n' ){
                    ++number_of_n;
                }else{
                    ++number_of_chr;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart2() && started==0){
                    start += j;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd2() ){
                    end += j;
                    quePos -= number_of_n;
                    break;
                }
            }


            number_of_chr=0;
            number_of_n=0;
            int32_t ref_start = 1;
            for ( int32_t j=0; j<seq1_0_string.size(); ++j ){
                if ( seq1_0_string[j] != 'n' ){
                    ++number_of_chr;
                }else{
                    ++number_of_n;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart1() ){
                    ref_start += j;
                    refPos += number_of_n;
                    break;
                }
            }


//            std::cout << "number_of_chr:" << number_of_chr << " pairedSimilarFragments0[i].getStart1(): " << pairedSimilarFragments0[i].getStart1() << std::endl;
            ofile << queChr << "\t16\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
//            std::cout << "line 696 quePos:" << quePos<< std::endl;
//            std::cout << "line 697 refPos:" << refPos<< std::endl;
//            std::cout << "line 698 start:" << start<< std::endl;
//            std::cout << "line 699 end:" << end<< std::endl;

            int32_t que_i = start-1; // the index of reference
            int32_t ref_i = ref_start-1; // the index of query


            std::vector<uint32_t> newCigar;
            for( int32_t j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                uint32_t op = 0;
                uint32_t length = 1;
                for ( int32_t h=0; h<cigarLength; ++h ){
                    op = cigarType;
                    if( op==0 ){
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else if(op==1) {
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else{
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                    }
                    op = cigarType;
                    if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                        newCigar.push_back(length << 4 | op);
                    }else{
                        newCigar[newCigar.size()-1] += length<<4;
                    }
                }
            }
            for( int32_t j=0; j<newCigar.size(); ++j ){
                uint32_t cigarLength = newCigar[j]>>4;
                uint32_t cigarType = newCigar[j]&0xf;
                ofile << cigarLength << "MIDDSHI=XB"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_0_string, start, end) << "\t*" << std::endl;


        }else if ( refStrand == 0 && queStrand ==1){


            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + seq1_0_string.size() -1 - pairedSimilarFragments0[i].getEnd1()+1;

            int32_t start = 1;
            int32_t end = 1;
            int32_t number_of_chr=0;
            int32_t number_of_n=0;
            int32_t started=0;

            for ( int32_t j=0; j<seq2_0_string.size(); ++j ){
                if ( seq2_0_string[j] == 'n' ){
                    ++number_of_n;
                }else{
                    ++number_of_chr;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart2() && started==0){
                    start += j;
                    quePos += number_of_n;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd2() ){
                    end += j;
                    break;
                }
            }
            number_of_chr=0;
            number_of_n=0;
            int32_t ref_start = 1;
            started=0;
            for ( int32_t j=0; j<seq1_0_string.size(); ++j ){
                if ( seq1_0_string[j] != 'n' ){
                    ++number_of_chr;
                }else{
                    ++number_of_n;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart1() && started==0){
                    ref_start += j;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd1() ){
                    refPos -= number_of_n;
                    break;
                }
            }


//            std::cout << "number_of_chr:" << number_of_chr << " pairedSimilarFragments0[i].getStart1(): " << pairedSimilarFragments0[i].getStart1() << std::endl;
            ofile << queChr << "\t16\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
//            std::cout << "line 696 quePos:" << quePos<< std::endl;
//            std::cout << "line 697 refPos:" << refPos<< std::endl;
//            std::cout << "line 698 start:" << start<< std::endl;
//            std::cout << "line 699 end:" << end<< std::endl;

            int32_t que_i = start-1; // the index of reference
            int32_t ref_i = ref_start-1; // the index of query
/*
            for( int32_t j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ) {
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j] >> 4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j] & 0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile << "\t";
*/
            std::vector<uint32_t> newCigar;
            for( int32_t j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                uint32_t op = 0;
                uint32_t length = 1;
                for ( int32_t h=0; h<cigarLength; ++h ){
                    op = cigarType;
                    if( op==0 ){
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else if(op==1) {
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else{
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                    }
                    op = cigarType;
                    if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                        newCigar.push_back(length << 4 | op);
                    }else{
                        newCigar[newCigar.size()-1] += length<<4;
                    }
                }
            }


            for( int32_t j=newCigar.size()-1; j>=0; --j ){
                uint32_t cigarLength = newCigar[j]>>4;
                uint32_t cigarType = newCigar[j]&0xf;
                ofile << cigarLength << "MIDDSHI=XB"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_0_string, start, end, 0) << "\t*" << std::endl;
        }else if ( refStrand == 0 && queStrand ==0){   // it seems this one is not correct

            int32_t quePos = queStart + seq2_0_string.size()-1-pairedSimilarFragments0[i].getEnd2()+1;
            int32_t refPos = refStart + seq1_0_string.size() -1 - pairedSimilarFragments0[i].getEnd1()+1;



            int32_t start = 1;
            int32_t end = 1;
            int32_t number_of_chr=0;
            int32_t number_of_n=0;
            int32_t started=0;

            for ( int32_t j=0; j<seq2_0_string.size(); ++j ){
                if ( seq2_0_string[j] == 'n' ){
                    ++number_of_n;
                }else{
                    ++number_of_chr;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart2() && started==0){
                    start += j;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd2() ){
                    end += j;
                    quePos -= number_of_n;
                    break;
                }
            }

            number_of_chr=0;
            number_of_n=0;
            int32_t ref_start = 1;
            started=0;
            for ( int32_t j=0; j<seq1_0_string.size(); ++j ){
                if ( seq1_0_string[j] != 'n' ){
                    ++number_of_chr;
                }else{
                    ++number_of_n;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart1() && started==0){
                    ref_start += j;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd1() ){
                    refPos -= number_of_n;
                    break;
                }
            }

            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";

            int32_t que_i = start-1; // the index of reference
            int32_t ref_i = ref_start-1; // the index of query

            std::vector<uint32_t> newCigar;
            for( int32_t j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                uint32_t op = 0;
                uint32_t length = 1;
                for ( int32_t h=0; h<cigarLength; ++h ){
                    op = cigarType;
                    if( op==0 ){
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else if(op==1) {
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else{
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                    }
                    op = cigarType;
                    if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                        newCigar.push_back(length << 4 | op);
                    }else{
                        newCigar[newCigar.size()-1] += length<<4;
                    }
                }
            }

            for( int32_t j=newCigar.size()-1; j>=0; --j ){
                uint32_t cigarLength = newCigar[j]>>4;
                uint32_t cigarType = newCigar[j]&0xf;
                ofile << cigarLength << "MIDDSHI=XB"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_0_string, start, end, 0) << "\t*" << std::endl;
        }
        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
               << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

    }
    ofile.close();
    ofile2.close();
}


void cut1Gap(std::string & _input,  std::string & _reference, std::string & _output,
                   int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty1,
                   int32_t & _extendGapPenalty1,  int32_t & _seed_window_size, int32_t & _mini_cns_score,
                   int32_t & _step_size, int32_t & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
                   const double & _kValue, int32_t & _w, const int32_t & _xDrop,
                   const std::string & refFasta, const std::string & queryFasta, double & _pvalues){

    Scorei m(_matchingScore, _mismatchingPenalty);

    FastaWithIndex refFastaWithIndex(refFasta);
    FastaWithIndex queryFastaWithIndex(queryFasta);

    std::map<std::string, FastaMeta> metaInformations;
    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( _input, sequences, metaInformations, seqNames );

    std::map<std::string, int32_t>  chrSize;
    for( std::map<std::string, FastaIndexEntry>::iterator it = refFastaWithIndex.getFastaIndexEntryImpl().getEntries().begin(); it!=refFastaWithIndex.getFastaIndexEntryImpl().getEntries().end(); ++it ){
        chrSize[it->first] = it->second.getLength();
    }
    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string refChr = metaInformations[_reference].getChr();

    std::string seq1_0_string = refFastaWithIndex.getSubSequence( refChr, metaInformations[_reference].getStart(),  metaInformations[_reference].getEnd(),  metaInformations[_reference].getStrand());  //getSubSequence( sequences0[refChr], metaInformations[_reference].getStart(),  metaInformations[_reference].getEnd(),  metaInformations[_reference].getStrand());
    std::string seq1_string=(seq1_0_string);
    seq1_string.erase(std::remove(seq1_string.begin(), seq1_string.end(), 'n'), seq1_string.end());
    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int32_t length1 = seq1_string.length();

    for( std::string _query : seqNames ){
        if( _query.compare(_reference) != 0) {
            std::ofstream ofile;
            ofile.open(_output+metaInformations[_query].getSpecies()+".sam", std::ios_base::app);
            std::ofstream ofile2;
            ofile2.open(_output+metaInformations[_query].getSpecies()+"o", std::ios_base::app);
            ofile.close();
            ofile2.close();

            if( sequences[_query].size() > 0 ){
//            std::string seq2_0_string = getSubSequence( sequences1[metaInformations[_query].getChr()], metaInformations[_query].getStart(),  metaInformations[_query].getEnd(),  metaInformations[_query].getStrand());
                std::string queryChr = metaInformations[_query].getChr();
                std::string seq2_0_string = queryFastaWithIndex.getSubSequence(queryChr, metaInformations[_query].getStart(),  metaInformations[_query].getEnd(),  metaInformations[_query].getStrand() );
                std::string seq2_string=std::string(seq2_0_string);
                seq2_string.erase(std::remove(seq2_string.begin(), seq2_string.end(), 'n'), seq2_string.end());
                int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
                int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
                int32_t length2 = seq2_string.length();

                std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                      seq2, seq2_rev_com, length1, length2, _seed_window_size,
                                                      _mini_cns_score, _matrix_boundary_distance, _openGapPenalty1,
                                                      _extendGapPenalty1, _matchingScore,
                                                      _mismatchingPenalty, m, _step_size, seq1_string, seq2_string,
                                                      _pvalues, _lambda, _kValue, _w, _xDrop);

                if (_onlySyntenic) {
                    std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                    pairedSimilarFragments0 = pairedSimilarFragments;
                }

                outputMasking(_output+metaInformations[_query].getSpecies(), pairedSimilarFragments0, metaInformations[_reference].getStrand(),
                        metaInformations[_query].getStrand(), metaInformations[_reference].getStart(),
                        metaInformations[_query].getStart(), metaInformations[_reference].getChr(), metaInformations[_query].getChr(),
                        seq1_string, seq2_string, seq1_0_string, seq2_0_string);
                delete seq2;
                delete seq2_rev_com;
            }

            ofile2.open(_output+metaInformations[_query].getSpecies()+"o", std::ios_base::app);
            ofile2 << "#done" << std::endl;
            ofile2.close();

        }
    }
    delete seq1;
    delete seq1_rev_com;
}

int cut1Gap(int argc, char** argv){
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" cut1Gap -i input -r reference -q query -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    input file in fasta format" << std::endl <<
          " -r STRING  reference sequence entry name in input fasta file" << std::endl <<
          " -ra STRING masked reference genome sequence in fasta format" << std::endl <<
          " -qa STRING masked query genome sequence in fasta format" << std::endl <<
          //        " -q STRING  query sequence entry name in input fasta file" << std::endl <<
          " -o FILE    output file, exists file will be overwrote" << std::endl <<
          " -y bool    only output syntenic result (default: false)" << std::endl <<
          " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT    open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT    extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -u INT     xextend alignment band width (default: " << w << ")" << std::endl <<
          //" - INT     xextend alignment band width (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -w INT     the windows size used to run the smith-waterman algorithm to get the alignment seed (default: "<<seed_window_size<<")" << std::endl <<
          //" -m INT     minimum seeds size to trigger a alignment extension and the minimum conserved sequence to report (default: " << mini_cns_seed_size << ")" << std::endl <<
          " -c INT     minimum seeds score to trigger a alignment extension (default: " << mini_cns_score << ")" << std::endl <<
          " -s INT     step size for sliding the smith-waterman seeds alignment window (default: " << step_size << ")" << std::endl <<
          " -p DOUBLE  pvalue for significant alignment output (default: " << pvalues << ")" << std::endl <<
/*          " -b INT     minimum number of base-pairs between the maximum smith-waterman score with the score matrix boundary (default: " << matrix_boundary_distance << ") , it is not important" << */std::endl;
    InputParser inputParser (argc, argv);

    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O1") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
    }
    if( inputParser.cmdOptionExists("-E1") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
    }
    if( inputParser.cmdOptionExists("-u") ){
        w = std::stoi( inputParser.getCmdOption("-u") );
    }
    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }
    if( inputParser.cmdOptionExists("-w") ){
        seed_window_size = std::stoi( inputParser.getCmdOption("-w") );
    }
    if( inputParser.cmdOptionExists("-c") ){
        mini_cns_score = std::stoi( inputParser.getCmdOption("-c") );
    }
    if( inputParser.cmdOptionExists("-s") ){
        step_size = std::stoi( inputParser.getCmdOption("-s") );
    }
    if( inputParser.cmdOptionExists("-b") ){
        matrix_boundary_distance = std::stoi( inputParser.getCmdOption("-b") );
    }
    if( inputParser.cmdOptionExists("-p") ){
        pvalues = std::stod( inputParser.getCmdOption("-p") );
    }

    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") &&
        /*inputParser.cmdOptionExists("-q") &&*/ inputParser.cmdOptionExists("-o") &&
            inputParser.cmdOptionExists("-ra") && inputParser.cmdOptionExists("-qa") ){
        std::string input = inputParser.getCmdOption("-i");
        std::string reference = inputParser.getCmdOption("-r");
        //std::string query = inputParser.getCmdOption("-q");
        std::string output = inputParser.getCmdOption("-o");
        std::string refFasta = inputParser.getCmdOption("-ra");
        std::string queryFasta = inputParser.getCmdOption("-qa");

        cut1Gap(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                      extendGapPenalty1, seed_window_size, mini_cns_score,
                      step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, xDrop,
                      refFasta, queryFasta, pvalues);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}




void cut2Gaps(std::string & _input,  std::string & _reference, std::string & _output,
             int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty1,
             int32_t & _extendGapPenalty1, int32_t & _openGapPenalty2,
             int32_t & _extendGapPenalty2, int32_t & _seed_window_size, int32_t & _mini_cns_score,
             int32_t & _step_size, int32_t & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
             const double & _kValue, int32_t & _w, const int32_t & _bandwidth, const int32_t & _xDrop, const int32_t & _zDrop,
             const std::string & refFasta, const std::string & queryFasta, double & _pvalues){

    Scorei m(_matchingScore, _mismatchingPenalty);


    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    readFastaFile( refFasta, sequences0, seqNames0);

    std::map<std::string, std::string> sequences1;
    std::vector<std::string> seqNames1;
    readFastaFile( queryFasta, sequences1, seqNames1);

    std::map<std::string, FastaMeta> metaInformations;
    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( _input, sequences, metaInformations, seqNames );

    std::map<std::string, int32_t>  chrSize;
    for( std::map<std::string, std::string>::iterator it = sequences0.begin(); it!=sequences0.end(); ++it ){
        chrSize[it->first] = it->second.size();
    }
    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string refChr = metaInformations[_reference].getChr();

    std::string seq1_0_string = getSubSequence( sequences0[refChr], metaInformations[_reference].getStart(),  metaInformations[_reference].getEnd(),  metaInformations[_reference].getStrand());
    std::string seq1_string=seq1_0_string;
    seq1_string.erase(std::remove(seq1_string.begin(), seq1_string.end(), 'n'), seq1_string.end());
    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int32_t length1 = seq1_string.length();
    sequences0.clear();

    for( std::string _query : seqNames ){
        if( _query.compare(_reference) != 0) {

            std::string seq2_0_string = getSubSequence( sequences1[metaInformations[_query].getChr()], metaInformations[_query].getStart(),  metaInformations[_query].getEnd(),  metaInformations[_query].getStrand());
            std::string seq2_string=seq2_0_string;
            seq2_string.erase(std::remove(seq2_string.begin(), seq2_string.end(), 'n'), seq2_string.end());
            int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
            int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
            int32_t length2 = seq2_string.length();

            std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
                    findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                           _seed_window_size, _mini_cns_score,
                                                           _matrix_boundary_distance, _openGapPenalty1, _extendGapPenalty1,
                                                           _openGapPenalty2, _extendGapPenalty2, _matchingScore,
                                                           _mismatchingPenalty, m, _step_size, seq1_string, seq2_string, _pvalues,
                                                           _lambda, _kValue, _zDrop, _bandwidth, _w, _xDrop);

            if (_onlySyntenic) {
                std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                pairedSimilarFragments0 = pairedSimilarFragments;
            }

            outputMasking(_output+metaInformations[_query].getSpecies(), pairedSimilarFragments0, metaInformations[_reference].getStrand(),
                          metaInformations[_query].getStrand(), metaInformations[_reference].getStart(),
                          metaInformations[_query].getStart(), metaInformations[_reference].getChr(), metaInformations[_query].getChr(),
                          seq1_string, seq2_string, seq1_0_string, seq2_0_string);
            delete seq2;
            delete seq2_rev_com;
        }
    }
    delete seq1;
    delete seq1_rev_com;
}

int cut2Gaps(int argc, char** argv){
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" cut2Gaps -i input -r reference -q query -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    input file in fasta format" << std::endl <<
          " -r STRING  reference sequence entry name in input fasta file" << std::endl <<
          " -ra STRING masked reference genome sequence in fasta format" << std::endl <<
          " -qa STRING masked query genome sequence in fasta format" << std::endl <<
          //        " -q STRING  query sequence entry name in input fasta file" << std::endl <<
          " -o FILE    output file, exists file will be overwrote" << std::endl <<
          " -y bool    only output syntenic result (default: false)" << std::endl <<
          " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT    open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT    extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -u INT     xextend alignment band width (default: " << w << ")" << std::endl <<
          " -O2 INT    open gap penalty (default: " << openGapPenalty2 << ")" << std::endl <<
          " -E2 INT    extend gap penalty (default: " << extendGapPenalty2 << ")" << std::endl <<
          " -u2 INT    z drop alignment band width (default: " << bandwidth << ")" << std::endl <<
          //" - INT     xextend alignment band width (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -w INT     the windows size used to run the smith-waterman algorithm to get the alignment seed (default: "<<seed_window_size<<")" << std::endl <<
          //" -m INT     minimum seeds size to trigger a alignment extension and the minimum conserved sequence to report (default: " << mini_cns_seed_size << ")" << std::endl <<
          " -c INT     minimum seeds score to trigger a alignment extension (default: " << mini_cns_score << ")" << std::endl <<
          " -s INT     step size for sliding the smith-waterman seeds alignment window (default: " << step_size << ")" << std::endl <<
          " -p DOUBLE  pvalue for significant alignment output (default: " << pvalues << ")" << std::endl <<
          " -b INT     minimum number of base-pairs between the maximum smith-waterman score with the score matrix boundary (default: " << matrix_boundary_distance << ") , it is not important" << std::endl;
    InputParser inputParser (argc, argv);

    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O1") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
    }
    if( inputParser.cmdOptionExists("-E1") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
    }
    if( inputParser.cmdOptionExists("-O2") ){
        openGapPenalty2 = std::stoi( inputParser.getCmdOption("-O2") );
    }
    if( inputParser.cmdOptionExists("-E2") ){
        extendGapPenalty2 = std::stoi( inputParser.getCmdOption("-E2") );
    }
    if( inputParser.cmdOptionExists("-u") ){
        w = std::stoi( inputParser.getCmdOption("-u") );
    }
    if( inputParser.cmdOptionExists("-u2") ){
        bandwidth = std::stoi( inputParser.getCmdOption("-u2") );
    }
    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }
    if( inputParser.cmdOptionExists("-w") ){
        seed_window_size = std::stoi( inputParser.getCmdOption("-w") );
    }
    if( inputParser.cmdOptionExists("-c") ){
        mini_cns_score = std::stoi( inputParser.getCmdOption("-c") );
    }
    if( inputParser.cmdOptionExists("-s") ){
        step_size = std::stoi( inputParser.getCmdOption("-s") );
    }
    if( inputParser.cmdOptionExists("-b") ){
        matrix_boundary_distance = std::stoi( inputParser.getCmdOption("-b") );
    }
    if( inputParser.cmdOptionExists("-p") ){
        pvalues = std::stod( inputParser.getCmdOption("-p") );
    }

    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") &&
              /*inputParser.cmdOptionExists("-q") &&*/ inputParser.cmdOptionExists("-o") &&
              inputParser.cmdOptionExists("-ra") && inputParser.cmdOptionExists("-qa") ){
        std::string input = inputParser.getCmdOption("-i");
        std::string reference = inputParser.getCmdOption("-r");
        //std::string query = inputParser.getCmdOption("-q");
        std::string output = inputParser.getCmdOption("-o");
        std::string refFasta = inputParser.getCmdOption("-ra");
        std::string queryFasta = inputParser.getCmdOption("-qa");


        cut2Gaps(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                extendGapPenalty1, openGapPenalty2, extendGapPenalty2, seed_window_size, mini_cns_score,
                step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, bandwidth,
                xDrop, zDrop, refFasta, queryFasta, pvalues);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}




void mapCNSToGenome(std::string & _input,  std::string & _reference, std::string & _output,
              int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty1,
              int32_t & _extendGapPenalty1, int32_t & _openGapPenalty2,
              int32_t & _extendGapPenalty2){

    int32_t refStrand = 1;
    int32_t queStrand = 1;

    Scorei m(_matchingScore, _mismatchingPenalty);

    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    readFastaFile( _reference, sequences0, seqNames0);

    std::map<std::string, std::string> sequences1;
    std::vector<std::string> seqNames1;
    readFastaFile( _input, sequences1, seqNames1);
    SequenceCharToUInt8 sequenceCharToUInt8;

    for( std::string _ref : seqNames0 ){
        std::string seq1_string = sequences0[_ref];
        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
        int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
        int32_t length1 = seq1_string.length();
        for( std::map<std::string, std::string>::iterator it=sequences1.begin(); it!=sequences1.end(); ++it ){
            std::string _query = it->first;
            std::string seq2_string = sequences1[_query];
            int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
            int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
            int32_t length2 = seq2_string.length();
            std::cout << _ref << " " << _query << std::endl;
            PairedSimilarFragment pairedSimilarFragments = mapCNSToGenome (  seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                                _openGapPenalty1, _extendGapPenalty1,
                                                                _openGapPenalty2, _extendGapPenalty2,
                                                                matchingScore, mismatchingPenalty, m);

            outputNoMasking( _output, pairedSimilarFragments, refStrand, queStrand, 1, 1, _ref, _query,
                                 seq1_string, seq2_string);

            delete seq2;
            delete seq2_rev_com;
        }

        delete seq1;
        delete seq1_rev_com;
    }
}

int mapCNSToGenome(int argc, char** argv){
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" cut2Gaps -i input -r reference -q query -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    CNS sequence file in fasta format" << std::endl <<
          " -r STRING masked reference genome sequence in fasta format" << std::endl <<
          " -o FILE    output file, exists file will be overwrote" << std::endl <<
          " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT    open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT    extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -O2 INT    open gap penalty (default: " << openGapPenalty2 << ")" << std::endl <<
          " -E2 INT    extend gap penalty (default: " << extendGapPenalty2 << ")" << std::endl << std::endl;
    InputParser inputParser (argc, argv);

    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O1") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
    }
    if( inputParser.cmdOptionExists("-E1") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
    }
    if( inputParser.cmdOptionExists("-O2") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O2") );
    }
    if( inputParser.cmdOptionExists("-E2") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E2") );
    }

    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") &&
              /*inputParser.cmdOptionExists("-q") &&*/ inputParser.cmdOptionExists("-o") ){
        std::string input = inputParser.getCmdOption("-i");
        std::string reference = inputParser.getCmdOption("-r");
        std::string output = inputParser.getCmdOption("-o");


        mapCNSToGenome(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                 extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}


int cut2Gaps2(int argc, char** argv){
    int32_t maximumAlignLength = 140000;
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" cut2Gaps2 -i input -r reference -q query -o output" << std::endl<<
          "          the default parameters use about 20GB RAM" << std::endl <<
          "Options" << std::endl <<
          " -ra STRING masked reference genome sequence in fasta format" << std::endl <<
          " -qa STRING masked query genome sequence in fasta format" << std::endl <<
          " -sa STRING sam file from the cut1Gap " << std::endl <<
          " -o  FILE   output file, exists file will be overwrote" << std::endl <<
          " -y  bool   only output syntenic result (default: false)" << std::endl <<
          " -A  INT    Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B  INT    Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT    open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT    extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -O2 INT    open gap penalty (default: " << openGapPenalty2 << ")" << std::endl <<
          " -E2 INT    extend gap penalty (default: " << extendGapPenalty2 << ")" << std::endl <<
          " -u2 INT    z drop alignment band width (default: " << bandwidth << ")" << std::endl <<
          " -m  INT    z drop alignment band width (default: " << maximumAlignLength << ")" << std::endl <<
          " -w  INT    the windows size used to run the smith-waterman algorithm to get the alignment seed (default: "<<seed_window_size<<")" << std::endl;
    InputParser inputParser (argc, argv);

    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O1") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
    }
    if( inputParser.cmdOptionExists("-E1") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
    }
    if( inputParser.cmdOptionExists("-O2") ){
        openGapPenalty2 = std::stoi( inputParser.getCmdOption("-O2") );
    }
    if( inputParser.cmdOptionExists("-E2") ){
        extendGapPenalty2 = std::stoi( inputParser.getCmdOption("-E2") );
    }
    if( inputParser.cmdOptionExists("-m") ){
        maximumAlignLength = std::stoi( inputParser.getCmdOption("-m") );
    }
    if( inputParser.cmdOptionExists("-u2") ){
        bandwidth = std::stoi( inputParser.getCmdOption("-u2") );
    }
    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }
    if( inputParser.cmdOptionExists("-w") ){
        seed_window_size = std::stoi( inputParser.getCmdOption("-w") );
    }
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-sa") && inputParser.cmdOptionExists("-o") &&
              inputParser.cmdOptionExists("-ra") && inputParser.cmdOptionExists("-qa") ){
        std::string samFile = inputParser.getCmdOption("-sa");
        std::string output = inputParser.getCmdOption("-o");
        std::string refFasta = inputParser.getCmdOption("-ra");
        std::string queryFasta = inputParser.getCmdOption("-qa");

        std::map<std::string, std::string> query_genome;
        std::vector<std::string> query_seqNames;
        readFastaFile( queryFasta, query_genome, query_seqNames);
        std::cout << "query genome reading done" << std::endl;

        std::map<std::string, std::string> ref_genome;
        std::vector<std::string> ref_seqNames;
        readFastaFile( refFasta, ref_genome, ref_seqNames);
        std::cout << "reference genome reading done" << std::endl;

        std::map<std::string, std::map<std::string, std::vector<Seed>>> positiveSeeds;
        std::map<std::string, std::map<std::string, std::vector<Seed>>> negativeSeeds;

        SequenceCharToUInt8 sequenceCharToUInt8;
        samFileToSeeds ( samFile, positiveSeeds, negativeSeeds, ref_genome, query_genome);
        std::cout << "sam file reading done" << std::endl;
        Scorei m(matchingScore, mismatchingPenalty);

        Matrix T(maximumAlignLength+1, maximumAlignLength + 1);
        for ( std::map<std::string, std::map<std::string, std::vector<Seed>>>::iterator it = positiveSeeds.begin(); it!=positiveSeeds.end(); ++it ){

            std::string seq1_0_string = ref_genome[it->first];
            std::string seq1_string=seq1_0_string;
            seq1_string.erase(std::remove(seq1_string.begin(), seq1_string.end(), 'n'), seq1_string.end());
            int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
            int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
            int32_t length1 = seq1_string.length();

            for ( std::map<std::string, std::vector<Seed>>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
                std::cout << "position maize:" << it->first << " sorghum: " << it2->first << std::endl;
                std::string seq2_0_string = query_genome[it2->first];
                std::string seq2_string=seq2_0_string;
                seq2_string.erase(std::remove(seq2_string.begin(), seq2_string.end(), 'n'), seq2_string.end());
                int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
                int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
                int32_t length2 = seq2_string.length();

                std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence (
                        seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,  seed_window_size,
                        openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, matchingScore,
                         mismatchingPenalty, m, seq1_string, seq2_string, pvalues,
                         lambda, kValue,zDrop,  bandwidth,positiveSeeds[it->first][it2->first], T, maximumAlignLength);
                if (onlySyntenic) {
                    std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                    pairedSimilarFragments0 = pairedSimilarFragments;
                }

                outputMasking( output, pairedSimilarFragments0, 1, 1, 1, 1, it->first, it2->first,
                              seq1_string, seq2_string, seq1_0_string, seq2_0_string);
                delete seq2;
                delete seq2_rev_com;
            }
            delete seq1;
            delete seq1_rev_com;
        }

        for ( std::map<std::string, std::map<std::string, std::vector<Seed>>>::iterator it = negativeSeeds.begin(); it!=negativeSeeds.end(); ++it ){

            std::string seq1_0_string = ref_genome[it->first];
            std::string seq1_string=seq1_0_string;
            seq1_string.erase(std::remove(seq1_string.begin(), seq1_string.end(), 'n'), seq1_string.end());
            int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
            int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
            int32_t length1 = seq1_string.length();

            for ( std::map<std::string, std::vector<Seed>>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
                std::cout << "negative maize:" << it->first << " sorghum: " << it2->first << std::endl;
                std::string seq2_0_string = query_genome[it2->first];
                seq2_0_string = getReverseComplementary(seq2_0_string);

                std::string seq2_string=seq2_0_string;
                seq2_string.erase(std::remove(seq2_string.begin(), seq2_string.end(), 'n'), seq2_string.end());
                int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
                int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());
                int32_t length2 = seq2_string.length();

                std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
                        findSimilarFragmentsForPairedSequence ( seq1, seq1_rev_com,
                                        seq2, seq2_rev_com, length1, length2,  seed_window_size,
                                        openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, matchingScore,
                                        mismatchingPenalty, m, seq1_string, seq2_string, pvalues,
                                        lambda, kValue,zDrop,  bandwidth,positiveSeeds[it->first][it2->first], T, maximumAlignLength);
                if (onlySyntenic) {
                    std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                    pairedSimilarFragments0 = pairedSimilarFragments;
                }

                outputMasking( output, pairedSimilarFragments0, 1, 1, 1, 1, it->first, it2->first,
                               seq1_string, seq2_string, seq1_0_string, seq2_0_string);
                delete seq2;
                delete seq2_rev_com;
            }
            delete seq1;
            delete seq1_rev_com;
        }
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}




int maskGenome(int argc, char** argv){
    int32_t maskKmerFrequency = 50;
    int32_t outputLineWidth = 60;
    double similarity = 0.8;
    char alterChar = 'n';
    bool ifSoftmask = false;
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" maskGenome -i input -r reference -q query -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    input file in fasta format" << std::endl <<
          " -o FILE    output file in fasta format" << std::endl <<
          " -g FILE    genome annotation file in GFF3 format" << std::endl <<
          " -z         mask the coding gene, default mask CDS (only useful if -g is used)" << std::endl <<
          " -s FILE    sam file (conflict with -g)" << std::endl <<
          " -c FILE    cds sequence file (must present if -s is used)" << std::endl <<
          " -p FLOAT   similarity of CDS alignment for masking (only useful if -s is used)" << std::endl <<
          " -k FILE    k-mer count file from jellyfish dump command" << std::endl <<
          " -l CHAR    alter character for of masked positions (default:" << alterChar << ")" << std::endl <<
          " -d         soft mask. Option -l has no effect when this option is in use." << std::endl <<
          " -f INT     k-mer frequency threshold for masking (default:" << maskKmerFrequency << ")" << std::endl <<
          " -w INT     output fasta file line width (default:" << outputLineWidth << ")" << std::endl;
    InputParser inputParser (argc, argv);

    if( inputParser.cmdOptionExists("-f") ){
        maskKmerFrequency = std::stoi( inputParser.getCmdOption("-f") );
    }
    if( inputParser.cmdOptionExists("-w") ){
        outputLineWidth = std::stoi( inputParser.getCmdOption("-w") );
    }
    if( inputParser.cmdOptionExists("-l") ){
        alterChar = inputParser.getCmdOption("-l")[0];
    }

    if( inputParser.cmdOptionExists("-d") ){
        ifSoftmask = true;
    }
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if ( inputParser.cmdOptionExists("-s") && inputParser.cmdOptionExists("-g") ){
        std::cerr << "please give one and only one of gff3 file or bam file" << std::endl;
    }else if(
            inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-o")){

        std::string input = inputParser.getCmdOption("-i");
        std::string outputFile = inputParser.getCmdOption("-o");
        std::string kmerProfileFile;
        if(inputParser.cmdOptionExists("-k")){
            kmerProfileFile = inputParser.getCmdOption("-k");
        }else{
            maskKmerFrequency = -1;
        }
        std::map<std::string, std::string>genome;
        std::vector<std::string> seqNames;
        readFastaFile( input, genome, seqNames);

        std::map<std::string, std::string> ifCds;
        if( inputParser.cmdOptionExists("-g") ){
            std::cout << "masking genome using GFF3 file and sam file" << std::endl;
            std::string gff = inputParser.getCmdOption("-g");
            if ( inputParser.cmdOptionExists("-z") ){
                gffToMaskGene ( gff, genome, ifCds);
            }else{
                gffToMask ( gff, genome, ifCds);
            }
        }else if (inputParser.cmdOptionExists("-s") and inputParser.cmdOptionExists("-c") ){
            std::cout << "masking genome using kmer and sam file" << std::endl;
            std::string sam = inputParser.getCmdOption("-s");

            if ( inputParser.cmdOptionExists("-p")) {
                similarity = std::stod( inputParser.getCmdOption("-p") );
            }
            std::string cds = inputParser.getCmdOption("-c");
            std::map<std::string, std::string> cdsSequences;
            std::vector<std::string> cdsNames;
            readFastaFile( cds, cdsSequences, cdsNames);
            samToMask (sam, genome, ifCds, similarity, cdsSequences);

        }

        maskGenome( genome, alterChar,  kmerProfileFile, maskKmerFrequency, outputLineWidth,
                ifCds, outputFile, seqNames, ifSoftmask);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}



int smitherWaterManScoreOfRandomFragments(int argc, char** argv){
    int32_t length = 1000;
    int32_t permutationTimes = 10000;
    int32_t seed = 1;
    std::stringstream usage;
    bool removen = true;
    usage << "Usage: "<<softwareName<<" ranSco -i input -r reference -q query -o output" << std::endl<<
          "Options" << std::endl <<
          " -r FILE    reference genome file in fasta format" << std::endl <<
          " -i FILE    input file in fasta format" << std::endl <<
          " -l INT     random sequence fragment length (default: " << length << ")" << std::endl <<
          " -p INT     random time (default: " << permutationTimes << ")" << std::endl <<
          " -s INT     random seed (default: " << seed << ")" << std::endl <<
          " -n bool    remove n from the genome sequences (default: true)" << std::endl <<
          " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O INT     open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E INT     extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl;


    InputParser inputParser (argc, argv);

    if( inputParser.cmdOptionExists("-n") ){
        removen = str2bool( inputParser.getCmdOption("-n"), onlySyntenic );
    }

    if( inputParser.cmdOptionExists("-l") ){
        length = std::stoi( inputParser.getCmdOption("-l") );
    }
    if( inputParser.cmdOptionExists("-p") ){
        permutationTimes = std::stoi( inputParser.getCmdOption("-p") );
    }
    if( inputParser.cmdOptionExists("-s") ){
        seed = std::stoi( inputParser.getCmdOption("-s") );
    }
    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O") );
    }
    if( inputParser.cmdOptionExists("-E") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E") );
    }

    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-r") && inputParser.cmdOptionExists("-i") ){
        std::vector<std::string> queryFastas = inputParser.getCmdOptionMultipleParameters("-i");

        std::string referenceFasta = inputParser.getCmdOption("-r");


        permutationLsqLambdaK( referenceFasta, queryFastas, openGapPenalty1,  extendGapPenalty1, matchingScore,
                               mismatchingPenalty, length, permutationTimes, seed, removen);

        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}


int multCns( int argc, char** argv ){
    int32_t minimumNumberOfSpecies = 3;
    uint32_t mini_cns_size = 7;
    double outputWithMinimumLengthPercentage = 0.8;
    std::stringstream usage;

    usage << "Usage: "<<softwareName<<" mult-cns -i input -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    reference genome file in fasta format" << std::endl <<
          //" -g FILE    query genome files in fasta format" << std::endl <<
          " -s FILEs   sam files output of any pairwise approach" << std::endl <<
          " -o FILE    output file in fasta format" << std::endl <<
          " -m FILE    minimum output size (default:" << mini_cns_size << ")" << std::endl <<
          //" -c INT     minimum number of species" << std::endl <<
          " -e DOUBLE  minimum out put of the full out put fragment length" << std::endl <<
          " -p INT     the minimum number of non-reference species that should contain similar fragment for output (default:"<< minimumNumberOfSpecies<<")" << std::endl <<
          " -n bool    use only single sequence overlapping with reference (default: true)" << std::endl <<
          "            if setting it as false, it might cost a lot of RAM and CPU time" << std::endl <<
          " -y bool    only output syntenic result (default: false)" << std::endl;
    InputParser inputParser (argc, argv);

    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }

    if( inputParser.cmdOptionExists("-m") ){
        mini_cns_size = std::stoi( inputParser.getCmdOption("-m") );
    }

    if( inputParser.cmdOptionExists("-e") ){
        outputWithMinimumLengthPercentage = std::stod( inputParser.getCmdOption("-e") );
    }

    if( inputParser.cmdOptionExists("-p") ){
        minimumNumberOfSpecies = std::stoi( inputParser.getCmdOption("-p") );
    }

    bool onlyPickOneSequenceForEachSamForMSA  = true;
    if( inputParser.cmdOptionExists("-n") ){
        onlyPickOneSequenceForEachSamForMSA = str2bool( inputParser.getCmdOption("-n"), onlySyntenic );
    }

    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-o") && inputParser.cmdOptionExists("-s")  ){
        std::string referenceGenomeFile = inputParser.getCmdOption("-i");
        std::string output = inputParser.getCmdOption("-o");

        std::vector<std::string> files = inputParser.getCmdOptionMultipleParameters("-s");
        std::map<std::string, std::string> samFiles;
        for( std::string file : files ){
            std::string filename = file;

            const size_t last_slash_idx = filename.find_last_of("\\/");
            if (std::string::npos != last_slash_idx){
                filename.erase(0, last_slash_idx + 1);
            }

            const size_t period_idx = filename.rfind('.');
            if (std::string::npos != period_idx){
                filename.erase(period_idx);
            }
            if( samFiles.find(filename) != samFiles.end() ){
                std::cout << "please give each input sam file a different prefix" << std::endl;
                return 1;
            }
            samFiles[filename] = file;
        }
        getCnsForMultipleSpecies (  onlySyntenic, output,
                //std::map<std::string, std::string> & sequences, /*species, fastaFile*/
                samFiles,/*species, sameFile*/
                referenceGenomeFile, minimumNumberOfSpecies, mini_cns_size,
                outputWithMinimumLengthPercentage, onlyPickOneSequenceForEachSamForMSA);

        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}




int slideWindow( int argc, char** argv ){
    matchingScore = 2;
    mismatchingPenalty = -3;
    openGapPenalty1 = -4;
    extendGapPenalty1 = -2;

    int32_t windowsSize = 150;
    step_size = 50;


    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" slideWindow -i input -o output" << std::endl<<
          "Options" << std::endl <<
          " -i FILE    sequence file in fasta format" << std::endl <<
          " -o FILE    output file in fasta format" << std::endl <<
          " -A INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O INT     open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E INT     extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -w INT     the windows size (default: "<<seed_window_size<<")" << std::endl <<
          " -s INT     step size for sliding the smith-waterman alignment window (default: " << step_size << ")" << std::endl <<
          " -m FILE    minimum output sites" << std::endl << std::endl;
    InputParser inputParser (argc, argv);


    if( inputParser.cmdOptionExists("-y") ){
        onlySyntenic = str2bool( inputParser.getCmdOption("-y"), onlySyntenic );
    }

    if( inputParser.cmdOptionExists("-A") ){
        matchingScore = std::stoi( inputParser.getCmdOption("-A") );
    }
    if( inputParser.cmdOptionExists("-B") ){
        mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
    }
    if( inputParser.cmdOptionExists("-O") ){
        openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O") );
    }
    if( inputParser.cmdOptionExists("-E") ){
        extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E") );
    }
    if( inputParser.cmdOptionExists("-w") ){
        windowsSize = std::stoi( inputParser.getCmdOption("-w") );
    }
    if( inputParser.cmdOptionExists("-s") ){
        step_size = std::stoi( inputParser.getCmdOption("-s") );
    }

    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-o")  ){
        std::string filePath = inputParser.getCmdOption("-i");
        std::string output = inputParser.getCmdOption("-o");

        std::map<std::string, std::string> sequences;
        std::vector<std::string> seqNames;
        readFastaFile( filePath, sequences, seqNames );
        int8_t ** seqs = new int8_t * [seqNames.size()];
        SequenceCharToUInt8 sequenceCharToUInt8;
        std::vector<int32_t> lengths;
        for( int i=0; i <seqNames.size(); ++i ){
            seqs[i] = sequenceCharToUInt8.seq_to_int8(sequences[seqNames[i]]);
            lengths.push_back(sequences[seqNames[i]].length());
        }


        Scorei m(matchingScore, mismatchingPenalty);
        slidingWindowAlignment ( seqs, lengths, seqNames, windowsSize, openGapPenalty1, extendGapPenalty1, matchingScore,
                                 mismatchingPenalty, output, step_size, m);

        for ( int i=0; i<seqNames.size(); ++i ){
            delete[] seqs[i];
        }
        delete[] seqs;

        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}
