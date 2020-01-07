//
// Created by bs674 on 6/8/19.
//

#include "permutationLsqLambdaK.h"

size_t randomeNumber( size_t range ){
//    std::srand(time(NULL));
//    return (std::rand() % (range+1)); // return a random value with the range [0, range]
    return (std::rand() % (range+1));
}


std::string randomReferenceFragments( std::map<std::string, std::string> & sequences, std::vector<std::string> & seqNames, const size_t & length){
    std::string seqName = seqNames[randomeNumber(seqNames.size()-1)];
    std::string seq = "";
    if( sequences[seqName].size()>length ){
        size_t start = randomeNumber(sequences[seqName].size()-length) + 1;
        size_t end = start + length -1;
        seq = getSubSequence(sequences[seqName], start, end);
    }
    return seq;
}

std::string randomReferenceFragments2( std::map<std::string, std::map<std::string, std::string>> & sequences,
        std::map<std::string, std::vector<std::string>> & seqNames, std::vector<std::string> & queryFastas,
        const size_t & length){
    std::string spe = queryFastas[randomeNumber(queryFastas.size()-1)];
    return randomReferenceFragments(sequences[spe], seqNames[spe], length);
}

void permutationLsqLambdaK( std::string & referenceFasta, std::vector<std::string> & queryFastas,
        const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
        const int & mismatchingPenalty, int32_t & length, int & permutationTimes, int32_t & seed, bool & removen){

    std::srand(seed);

    std::map<std::string, std::string> referenceSequences;
    std::vector<std::string> referenceSeqNames;
    readFastaFile( referenceFasta, referenceSequences, referenceSeqNames );

    std::map<std::string, std::map<std::string, std::string>> querySequences;
    std::map<std::string, std::vector<std::string>> quyerySeqNames;

    for ( std::string queryFasta :  queryFastas){
        std::map<std::string, std::string> querySequence;
        std::vector<std::string> querySeqName;
        readFastaFile( queryFasta, querySequence, querySeqName );
        querySequences[queryFasta]=querySequence;
        quyerySeqNames[queryFasta]=querySeqName;
    }

    if( removen ){
        for( std::map<std::string, std::string>::iterator it=referenceSequences.begin(); it!=referenceSequences.end(); ++it ){
            std::string seq_string=(it->second);
            seq_string.erase(std::remove(seq_string.begin(), seq_string.end(), 'n'), seq_string.end());
            referenceSequences[it->first] = seq_string;
        }
        for ( std::string queryFasta :  queryFastas){
//            std::cout << "queryFasta:" << queryFasta << " size:" << querySequences[queryFasta].size() << std::endl;
            for( std::map<std::string, std::string>::iterator it=querySequences[queryFasta].begin(); it!=querySequences[queryFasta].end(); ++it ){
                std::string seq_string=(it->second);
                seq_string.erase(std::remove(seq_string.begin(), seq_string.end(), 'n'), seq_string.end());
                querySequences[queryFasta][it->first] = seq_string;
//                std::cout << "queryFasta:" << queryFasta << " size:" << querySequences[queryFasta].size() << " " << it->first << " " << seq_string.size() << std::endl;
            }
        }
    }


    Scorei m(matchingScore, mismatchingPenalty);
    SequenceCharToUInt8 sequenceCharToUInt8;


    for ( int i=0; i<permutationTimes; ++i ){
//        std::cout << "line 71" << std::endl;
        std::string referenceSeq = randomReferenceFragments( referenceSequences, referenceSeqNames, length);
        while ( (referenceSeq.length() <  length) || referenceSeq.find('n')!=std::string::npos || referenceSeq.find('N')!=std::string::npos || referenceSeq.find('-')!=std::string::npos ){
            referenceSeq = randomReferenceFragments( referenceSequences, referenceSeqNames, length);
        }
//        std::cout << "line 76" << std::endl;
        std::string querySeq = randomReferenceFragments2( querySequences, quyerySeqNames, queryFastas, length);
        while( querySeq.length() < length || querySeq.find('n')!=std::string::npos || querySeq.find('N')!=std::string::npos || querySeq.find('-')!=std::string::npos){
            querySeq = randomReferenceFragments2( querySequences, quyerySeqNames, queryFastas, length);
        }
//        std::cout << "line 81" << std::endl;
        int32_t length1 = referenceSeq.length();
        int32_t length2 = querySeq.length();
//        std::cout << "line 84" << std::endl;
        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(referenceSeq);
        int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);

        int32_t maxScore;

        int32_t endPosition1, endPosition2;
        bool returnPosition =false;
        SmithWaterman(seq1, seq2, length1, length2, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1, endPosition2, m, returnPosition);
        delete seq1;
        delete seq2;
        std::cout << maxScore << std::endl;
    }
}


// this one implemented the same function as last one
// but it is slower, since calling the full version of smith-waterman algorithm
void permutationLsqLambdaKslow( std::string & referenceFasta, std::vector<std::string> & queryFastas,
                            const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                            const int & mismatchingPenalty, int32_t & length, int & permutationTimes, int32_t & seed){

    std::srand(seed);

    std::map<std::string, std::string> referenceSequences;
    std::vector<std::string> referenceSeqNames;
    readFastaFile( referenceFasta, referenceSequences, referenceSeqNames );

    std::map<std::string, std::map<std::string, std::string>> querySequences;
    std::map<std::string, std::vector<std::string>> quyerySeqNames;

    for ( std::string queryFasta :  queryFastas){
        std::map<std::string, std::string> querySequence;
        std::vector<std::string> querySeqName;
        readFastaFile( queryFasta, querySequence, querySeqName );
        querySequences[queryFasta]=querySequence;
        quyerySeqNames[queryFasta]=querySeqName;
    }

    Scorei m(matchingScore, mismatchingPenalty);
    std::vector<int64_t> maximumScores;
    SequenceCharToUInt8 sequenceCharToUInt8;



    for ( int i=0; i<permutationTimes; ++i ){

        std::string referenceSeq = randomReferenceFragments( referenceSequences, referenceSeqNames, length);
        while ( (referenceSeq.length() <  length) || referenceSeq.find('n')!=std::string::npos || referenceSeq.find('N')!=std::string::npos || referenceSeq.find('-')!=std::string::npos ){
            referenceSeq = randomReferenceFragments( referenceSequences, referenceSeqNames, length);
        }

        std::string querySeq = randomReferenceFragments2( querySequences, quyerySeqNames, queryFastas, length);
        while( querySeq.length() < length || querySeq.find('n')!=std::string::npos || querySeq.find('N')!=std::string::npos || querySeq.find('-')!=std::string::npos){
            querySeq = randomReferenceFragments2( querySequences, quyerySeqNames, queryFastas, length);
        }

        int32_t length1 = referenceSeq.length();
        int32_t length2 = querySeq.length();
        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(referenceSeq);
        int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);

        int32_t maxScore;
        int32_t endPosition1;
        int32_t endPosition2;
        std::vector<uint32_t > cigar = SmithWaterman(seq1, seq2, length1, length2, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1, endPosition2, m, false, false);
        maximumScores.push_back(maxScore);
        std::cout << maxScore << std::endl;
    }
}
