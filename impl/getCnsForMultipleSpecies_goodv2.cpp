//
// Created by Baoxing song on 2019-01-09.
//

#include "getCnsForMultipleSpecies.h"

/**
 * The problem with this one is that the MSA output must have the reference sequence there
 * The alignment without reference would not be used
 *
 * the second problem is that the alignment between reference sequence and the first query sequence is used for coordinating
 * Even some alignments without overlap with any alignment between the reference sequence and the first query sequence could be outputted
 * */


/*
 * take to PairedSimilarFragment and check if they overlap with each other in the terms of reference range
 * return a bool value
 * **/

bool pairedSimilarFragmentsOverlap( const PairedSimilarFragment & pairedSimilarFragment1, const PairedSimilarFragment & pairedSimilarFragment2){
    if( pairedSimilarFragment1.getStart1()>=pairedSimilarFragment2.getStart1() && pairedSimilarFragment1.getStart1()<pairedSimilarFragment2.getEnd1() ){
        return true;
    }else if( pairedSimilarFragment2.getStart1()>=pairedSimilarFragment1.getStart1() && pairedSimilarFragment2.getStart1()<pairedSimilarFragment1.getEnd1() ){
        return true;
    }
    return false;
}


//the return is a matrix which lists all the combinations, like [[2,3,5], [4,5,7]]   2,3,5 is the index from allPairedSimilarFragments[1], allPairedSimilarFragments[2] and allPairedSimilarFragments[3]
// allPairedSimilarFragments[0] is the reference one, so do not use it here
std::vector< std::vector<int>> generateAllCombinations(std::vector<std::vector<bool>> & overlapped){
    std::vector< std::vector<int>> indexs;
    for ( int i=0; i < overlapped.size(); ++i ){
        indexs.push_back(std::vector<int>());
        for ( int j=0; j<overlapped[i].size(); ++j ) {
            if( overlapped[i][j] ){
                indexs[i].push_back(j);
            }
        }
        if( indexs[i].size() == 0 ){
            indexs[i].push_back(-100);
        }
    }
    return indexs;
}

bool getNext( const std::vector< std::vector<int>> & allCombinations, std::vector<int> & combination ){
    size_t i;
    if( combination[0] <0 ){
        for ( i=0; i<combination.size(); ++i ){
            combination[i] = 0;
        }
        return true;
    } else {
        i =0;
        while( i<allCombinations.size() ){
            if ( combination[i] < (allCombinations[i].size()-1) ){
                combination[i] = combination[i] + 1;
                return true;
            }else{
                combination[i]=0;
                ++i;
            }
        }
        return false;
    }
}

void getCnsForMultipleSpecies ( int8_t ** seqs, int8_t ** seq_revs, std::vector<uint32_t> & lengths,
                                uint32_t & windowsSize, uint32_t & mini_cns_seed_size, int32_t & mini_cns_score,
                                const int & matrix_boundary_distance, const int & _open_gap_penalty,
                                const int & _extend_gap_penalty, const int & matchingScore,
                                const int & mismatchingPenalty, const bool & onlySyntenic, const std::string & output,
                                std::map<std::string, std::string>& sequences, std::vector<std::string> & seqNames,
                                const uint32_t & step_size, std::vector<std::string> & seqs_string,
                                const int & minimumNumberOfSpecies, int8_t m[][5], const uint32_t & mini_cns_size){

    double outputWithLengthPercentage = 0.8;

    if( lengths.size()<3 ){
        std::cerr << "you should have at least 3 sequences in your input fasta file" << std::endl;
        return;
    }

    //pair-wise sequence alignment begin
    std::vector<std::vector<PairedSimilarFragment>> allPairedSimilarFragments(lengths.size()-1);
    // the first sequence is used as reference and do not perform pairwise sequence alignment without the reference
    int8_t * seq1 = seqs[0];
    int8_t * seq1_rev_com = seq_revs[0];
    uint32_t length1 = lengths[0];
    for ( int i=1; i< lengths.size(); ++i ){

        int8_t * seq2 = seqs[i];
        int8_t * seq2_rev_com = seq_revs[i];
        uint32_t length2 = lengths[i];

        std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                                                                            windowsSize, mini_cns_seed_size, mini_cns_score,
                                                                                                            matrix_boundary_distance, _open_gap_penalty, _extend_gap_penalty, matchingScore, mismatchingPenalty, m, step_size, seqs_string[0], seqs_string[i], 1.0);

        if ( onlySyntenic ){
            std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
            pairedSimilarFragments0 = pairedSimilarFragments;
        }
        allPairedSimilarFragments[i-1] = pairedSimilarFragments0;
    }
    //pair-wise sequence alignment end


    //sort according to the reference coordinate begin //remember those have been sorted, should be used for speeding up
    for ( int i=0; i< allPairedSimilarFragments.size(); ++i ){
        std::sort(allPairedSimilarFragments[i].begin(), allPairedSimilarFragments[i].end(), [](PairedSimilarFragment a, PairedSimilarFragment b) {
            return a.getStart1() < b.getStart1();
        });
    }
    //sort according to the reference coordinate end

    std::vector<std::vector<bool>> overlapped(allPairedSimilarFragments.size()-1);
    for ( int i=1; i< allPairedSimilarFragments.size(); ++i ) {
        overlapped[i-1].resize(allPairedSimilarFragments[i].size(), false);
    }
//    std::cout << " line 114" << std::endl;

    std::ofstream ofile;
    ofile.open(output);
    /**
     * If there are multiple records in the 3,4,...n species that overlapped with the record in the 2(second) one.
     * Firstly only care about the first overlapped one, then take the second third .... overlapping as reference to run everything again
     * **/
    for( int i=0; i<allPairedSimilarFragments[0].size(); ++i ){ // todo warning the alignment must present in the pair-wise alignment of first two sequences

        // set all the if overlap values begin
        for ( int j=1; j<allPairedSimilarFragments.size(); ++j ){
            for ( int k=0; k<allPairedSimilarFragments[j].size(); ++k ){
                if( pairedSimilarFragmentsOverlap(allPairedSimilarFragments[0][i], allPairedSimilarFragments[j][k] ) ) {
                    overlapped[j-1][k]=true;
                } else if ( allPairedSimilarFragments[j][k].getStart1() > allPairedSimilarFragments[0][i].getEnd1() ) { // take the advantage of the vector has been sorted to speed up
                    for(;k<allPairedSimilarFragments[j].size(); ++k){
                        overlapped[j-1][k] = false;
                    }
                } else {
                    overlapped[j-1][k] = false;
                }
            }
        }
        // set all the if overlap values end

        // design a recursion algorithm for this purpose
        std::vector< std::vector<int>> allCombinations = generateAllCombinations(overlapped);  // matrix size: ### X (n-2)



        //overlapped has the same size with allPairedSimilarFragments
        std::vector<int> combination (allCombinations.size(), -1);
        while( getNext( allCombinations, combination ) ){
            std::set<int> starts;
            std::set<int> ends;
            starts.insert( allPairedSimilarFragments[0][i].getStart1() );
            ends.insert( allPairedSimilarFragments[0][i].getEnd1() );
            for(int k=0; k<combination.size(); ++k ) {
                if (allCombinations[k][combination[k]] > -1) { // this is a real overlap
                    starts.insert( allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getStart1() );
                    ends.insert( allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getEnd1() );
                }
            }

            int refStart = 0;
            int refEnd = 0;
            int length = 0;
            int numberOfSpecies = 0;
            for ( int start : starts ) {
                for ( int end : ends ) {
                    int thisNumberOfSpecies = 0;
                    int newLength = end - start + 1;
                    if( length<=newLength && mini_cns_size<=newLength ){
                        if( allPairedSimilarFragments[0][i].getStart1()<=start && allPairedSimilarFragments[0][i].getEnd1() >= end ){
                            thisNumberOfSpecies++;
                        }
                        for( int k=0; k<combination.size(); ++k ) {
                            if (allCombinations[k][combination[k]] > -1) {
                                if( allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getStart1()<=start && allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getEnd1() >= end ){
                                    thisNumberOfSpecies++;
                                }
                            }
                        }
                        if ( thisNumberOfSpecies > minimumNumberOfSpecies /*&& thisNumberOfSpecies > numberOfSpecies*/ ){
                            numberOfSpecies = thisNumberOfSpecies;
                            refStart = start;
                            refEnd = end;
                            length = newLength;
                        }
                    }
                }
            }

            if( (refEnd > refStart) && ((refEnd - refStart + 1) >= mini_cns_size) ){
                int newRefStart = refEnd;
                int newRedEnd = refStart;

                std::set<std::string> species;
                std::set<int> badIndex;
                int thisStart = allPairedSimilarFragments[0][i].getStart1() > refStart ? allPairedSimilarFragments[0][i].getStart1() : refStart;
                int thisEnd = allPairedSimilarFragments[0][i].getEnd1() < refEnd ? allPairedSimilarFragments[0][i].getEnd1() : refEnd;
                double thisLength = thisEnd - thisStart + 1;
                if( thisLength < outputWithLengthPercentage*double(length) ){
                    badIndex.insert(-1);
                }else{
                    if ( thisStart< newRefStart){
                        newRefStart = thisStart;
                    }
                    if( thisEnd > newRedEnd ){
                        newRedEnd=thisEnd;
                    }
                    std::string spe = seqNames[1];   // todo reformat the meta line of fasta and check this line
                    spe = spe.substr(18,5);
                    species.insert(spe);
                }

                for( int k=0; k<combination.size(); ++k ){
                    if( allCombinations[k][combination[k]] > -1 ){ // this is a real overlap
                        thisStart = allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getStart1() > refStart ? allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getStart1() : refStart;
                        thisEnd = allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getEnd1() < refEnd ? allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getEnd1() : refEnd;
                        thisLength = thisEnd - thisStart + 1;
    //                    std::cout << " line 216 thisLength: " << thisLength << " length " << length << std::endl;
                        if( thisLength < outputWithLengthPercentage*double(length) ){
                            badIndex.insert(k);
                        }else{
                            if ( thisStart< newRefStart){
                                newRefStart = thisStart;
                            }
                            if( thisEnd > newRedEnd ){
                                newRedEnd=thisEnd;
                            }
                            std::string spe = seqNames[k+2];
                            spe = spe.substr(18,5);
                            species.insert(spe);
                        }
                    }
                }
                if( newRefStart > refStart ){
                    refStart = newRefStart;
                }
                if( newRedEnd < refEnd ){
                    refEnd = newRedEnd;
                }

                // std::cout << " line 226  " << i << " species.size " << species.size() << " refStart " << refStart << " refEnd " << refEnd << std::endl;
                //allPairedSimilarFragments[0] never output
                if ( /*allMatched*/ species.size()>=minimumNumberOfSpecies &&  (refEnd > refStart) && ((refEnd - refStart + 1) >= mini_cns_size) ) {
                    std::vector<std::string> alignmentNames;
                    std::vector<std::string> alignmentSeqs;

    //                std::cout << " line 232" << std::endl;
                    if( badIndex.find(-1) == badIndex.end() ){
    //                    std::cout << " line 234" << std::endl;

                        int start2 = allPairedSimilarFragments[0][i].getStart2();
                        int end2 = allPairedSimilarFragments[0][i].getEnd2();

                        int endReference = allPairedSimilarFragments[0][i].getEnd1();
                        int positionReference = allPairedSimilarFragments[0][i].getStart1();
                        --positionReference;

                        std::string alignSeq;
                        std::string refSeq;
                        int numberOfAltChar = 0;
                        bool alignRegionStart = false;
                        bool alignRegionEnd = true;
    //                    std::cout << " line 248" << std::endl;
                        for ( int w=0; w<allPairedSimilarFragments[0][i].getAlignment1().length(); ++w ){
                            if( allPairedSimilarFragments[0][i].getAlignment1()[w] != '-' ){
                                ++positionReference;
                            }
                            if( !alignRegionStart ){
    //                            std::cout << " line 254" << std::endl;
                                if( positionReference>refStart ){
                                    int h = positionReference - refStart;
                                    while( h>0 ){
                                        alignSeq += ".";
                                        refSeq += seqs_string[0][positionReference-h-1];
                                        --h;
                                    }
                                }
    //                            std::cout << " line 263" << std::endl;
                                if( positionReference>=refStart && (!alignRegionStart) ){
                                    alignRegionStart =true;
                                    start2 += numberOfAltChar;
                                    end2 = start2;
                                }
    //                            std::cout << " line 248" << std::endl;
                            }
                            if( alignRegionStart && alignRegionEnd ){
                                refSeq += allPairedSimilarFragments[0][i].getAlignment1()[w];
                                alignSeq += allPairedSimilarFragments[0][i].getAlignment2()[w];
                            }
                            if( allPairedSimilarFragments[0][i].getAlignment2()[w] != '-' ){
                                if( alignRegionStart && alignRegionEnd ){
                                    ++end2;
                                }
                                ++numberOfAltChar;
                            }
    //                        std::cout << " line 281" << std::endl;
                            if( alignRegionEnd ){
    //                            std::cout << " line 283" << std::endl;
                                if( positionReference == endReference && positionReference<refEnd){
                                    int h = refEnd - positionReference;
                                    for( int z=0; z<h; ++z ){
                                        alignSeq += ".";
                                        refSeq += seqs_string[0][positionReference+z];
                                    }
                                    alignRegionEnd = false;
                                }
    //                            std::cout << " line 292" << std::endl;
                                if( alignRegionEnd && positionReference>=refEnd ){
                                    alignRegionEnd = false;
                                }
    //                            std::cout << " line 296" << std::endl;
                            }
                        }
                        --end2;

                        alignmentNames.push_back(seqNames[0] + ":" + std::to_string(refStart) + "-" + std::to_string(refEnd)); //the reference one
                        alignmentNames.push_back(seqNames[1] + ":" + std::to_string(start2) + "-" + std::to_string(end2));

                        alignmentSeqs.push_back(refSeq);
                        alignmentSeqs.push_back(alignSeq);
                    }

    //                std::cout << " line 307" << std::endl;
                    for (int k = 0; k < combination.size(); ++k) {
                        if( allCombinations[k][combination[k]] > -1 && badIndex.find(k)==badIndex.end() ){
                            int start2 = allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getStart2();
                            int end2 = allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getEnd2();

                            int endReference = allPairedSimilarFragments[0][i].getEnd1();
                            int positionReference = allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getStart1();
                            --positionReference;

                            std::string alignSeq;
                            std::string refSeq;
                            int numberOfAltChar = 0;
                            bool alignRegionStart = false;
                            bool alignRegionEnd = true;
                            for ( int w=0; w<allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getAlignment1().length(); ++w ){
                                if( allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getAlignment1()[w] != '-' ){
                                    ++positionReference;
                                }
                                if( !alignRegionStart ){
                                    if( positionReference>refStart && (!alignRegionStart) ){
                                        int h = positionReference - refStart;
                                        while( h>0 ){
                                            alignSeq += ".";
                                            refSeq += seqs_string[0][positionReference-h-1];
                                            h--;
                                        }
                                    }
                                    if( positionReference>=refStart ){
                                        alignRegionStart =true;
                                        start2 += numberOfAltChar;
                                        end2 = start2;
                                    }
                                }
                                if( alignRegionStart && alignRegionEnd ){
                                    alignSeq += allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getAlignment2()[w];
                                    refSeq += allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getAlignment1()[w];
                                }
                                if( allPairedSimilarFragments[k+1][allCombinations[k][combination[k]]].getAlignment2()[w] != '-' ){
                                    if( alignRegionStart && alignRegionEnd ){
                                        ++end2;
                                    }
                                    ++numberOfAltChar;
                                }
                                if( alignRegionEnd ){
                                    if( positionReference == endReference &&  positionReference< refEnd){
                                        int h = refEnd - positionReference;
                                        for( int z=0; z<h; ++z ){
                                            alignSeq += ".";
                                            refSeq += seqs_string[0][positionReference+z];
                                        }
                                        alignRegionEnd = false;
                                    }
                                    if( positionReference>=refEnd ){
                                        alignRegionEnd = false;
                                    }
                                }
                            }
                            --end2;

    //                        std::cout << " line 366" << std::endl;

                            if ( alignmentSeqs.empty()  ){
                                alignmentNames.push_back(seqNames[0] + ":" + std::to_string(refStart) + "-" + std::to_string(refEnd)); //the reference one
                                alignmentSeqs.push_back(refSeq);
                            } else { //up date the alignment sequence by comparing the reference sequence
                                int lengthq = refSeq.length();
                                for ( int o=0; o<lengthq; ++o ){
                                    if( alignmentSeqs[0].length() > o ){ // this is always the reference one
                                        if( refSeq[o] == '-' && alignmentSeqs[0][o] !='-' ){
                                            for( int p=0; p<alignmentSeqs.size(); ++p ){
                                                alignmentSeqs[p].insert(o, "-");
                                            }
                                        }else if( refSeq[o] != '-' && alignmentSeqs[0][o] =='-' ){
                                            refSeq.insert(o, "-");
                                            alignSeq.insert(o, "-");
                                            ++lengthq;
                                        }
                                    }
                                }
                                while( refSeq.length() < alignmentSeqs[0].length() ){
                                    refSeq.insert(refSeq.length(), "-");
                                    alignSeq.insert(alignSeq.length(), "-");
                                }
                                while( refSeq.length() > alignmentSeqs[0].length() ){
                                    for( int p=0; p<alignmentSeqs.size(); ++p ){
                                        alignmentSeqs[p].insert(alignmentSeqs[p].length(), "-");
                                    }
                                }
                            }
    //                        std::cout << " line 395" << std::endl;
                            alignmentNames.push_back(seqNames[k+2] + ":" + std::to_string(start2) + "-" + std::to_string(end2));
                            alignmentSeqs.push_back(alignSeq);
                        }
                    }
                    for ( int o=0; o<alignmentNames.size(); ++o ){
                        ofile << ">" << alignmentNames[o]<< std::endl << alignmentSeqs[o] << std::endl;
                    }
                    ofile << std::endl << std::endl;
                }
            }
        }
    }
    ofile.close();
}
