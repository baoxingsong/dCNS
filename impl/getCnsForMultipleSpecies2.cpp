//
// Created by Baoxing song on 2019-01-09.
//

#include "getCnsForMultipleSpecies.h"

bool pairedSimilarFragmentsOverlap( PairedSimilarFragment pairedSimilarFragment1, PairedSimilarFragment pairedSimilarFragment2){
    if( pairedSimilarFragment1.getStart1()>=pairedSimilarFragment2.getStart1() && pairedSimilarFragment1.getStart1()<pairedSimilarFragment2.getEnd1() ){
        return true;
    }else if( pairedSimilarFragment2.getStart1()>=pairedSimilarFragment1.getStart1() && pairedSimilarFragment2.getStart1()<pairedSimilarFragment1.getEnd1() ){
        return true;
    }
    return false;
}

//the return is a matrix which lists all the combinations, like [[2,3,5], [4,5,7]]   2,3,5 is the index from allPairedSimilarFragments[1], allPairedSimilarFragments[2] and allPairedSimilarFragments[3]
// allPairedSimilarFragments[0] is the reference one, so do not use it here
std::vector< std::vector<int>> generateAllCombinations(std::vector<std::vector<bool>> & overlapped, std::vector<std::vector<PairedSimilarFragment>> & allPairedSimilarFragments, int index){
    if( (index+1) == (allPairedSimilarFragments.size()-1) ){ // this is the last one
        std::vector< std::vector<int>> indexs;
        for( int ii=0; ii<allPairedSimilarFragments[index+1].size(); ++ii ){
            if( overlapped[index][ii] ){
                indexs.push_back(std::vector<int>(0));
                indexs[indexs.size()-1].push_back(ii);
            }
        }
        if( indexs.size() == 0 ){
            indexs.push_back(std::vector<int>(0));
            indexs[indexs.size()-1].push_back(-100);
        }
        return indexs;
    } else {
        std::vector< std::vector<int>> indexs = generateAllCombinations( overlapped, allPairedSimilarFragments, index+1);
        std::vector< std::vector<int>> indexs2;
        int j;
        bool ifEverOverlap = false;
        for( int ii=0; ii<allPairedSimilarFragments[index+1].size(); ++ii ){
            if( overlapped[index][ii] ){
                ifEverOverlap = true;
                for( j=0; j<indexs.size(); ++j ){
                    indexs2.push_back(std::vector<int>(0));
                    indexs2[indexs2.size()-1].push_back(ii);
                    indexs2[indexs2.size()-1].insert(indexs2[indexs2.size()-1].end(), indexs[j].begin(), indexs[j].end());
                }
            }
        }
        if( !ifEverOverlap ){
            for( j=0; j<indexs.size(); ++j ){
                indexs2.push_back(std::vector<int>(0));
                indexs2[indexs2.size()-1].push_back(-100);
                indexs2[indexs2.size()-1].insert(indexs2[indexs2.size()-1].end(), indexs[j].begin(), indexs[j].end());
            }
        }
        return indexs2;
    }
}

void getCnsForMultipleSpecies(int8_t ** seqs, int8_t ** seq_revs, std::vector<uint32_t> & lengths,
                                uint32_t & windowsSize, uint32_t & mini_cns_size, int32_t & mini_cns_score,
                                const int & matrix_boundary_distance, const int & _open_gap_penalty,
                                const int & _extend_gap_penalty, const int & matchingScore,
                                const int & mismatchingPenalty, const std::string & output,
                                std::map<std::string, std::string>& sequences, std::vector<std::string> & seqNames,
                                const uint32_t & step_size, std::vector<std::string> & seqs_string,
                                const int & minimumNumberOfSpecies, int8_t m[][5]){
    if( lengths.size()<3 ){
        std::cerr << "you should have at least 3 sequences in your input fasta file" << std::endl;
        return;
    }
    int i, j, k;
    std::vector<std::vector<PairedSimilarFragment>> allPairedSimilarFragments(lengths.size()-1);
    // the first sequence is used as reference and do not perform pairwise sequence alignment without the reference
    int8_t * seq1 = seqs[0];
    int8_t * seq1_rev_com = seq_revs[0];
    uint32_t length1 = lengths[0];
    for ( i=1; i< lengths.size(); ++i ){
        int8_t * seq2 = seqs[i];
        int8_t * seq2_rev_com = seq_revs[i];
        uint32_t length2 = lengths[i];
        std::cout << " line 80 " << i << std::endl;

        allPairedSimilarFragments[i-1] = findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                                                windowsSize, mini_cns_size, mini_cns_score,
                                                                                matrix_boundary_distance, _open_gap_penalty, _extend_gap_penalty, matchingScore, mismatchingPenalty, m, step_size, seqs_string[0], seqs_string[i]);
    }

    for ( i=0; i< allPairedSimilarFragments.size(); ++i ){
        std::sort(allPairedSimilarFragments[i].begin(), allPairedSimilarFragments[i].end(), [](PairedSimilarFragment a, PairedSimilarFragment b) {
            return a.getStart1() < b.getStart1();
        });
    }

    std::vector<std::vector<bool>> overlapped(allPairedSimilarFragments.size()-1);
    for ( i=1; i< allPairedSimilarFragments.size(); ++i ) {
        for (j = 0; j < allPairedSimilarFragments[i].size(); ++j) {
            overlapped[i-1].push_back(false);
        }
    }
    std::cout << " line 99" << std::endl;
    std::ofstream ofile;
    ofile.open(output);
    /**
     * If there are multiple records in the 3,4,...n species that overlapped with the record in the 2(second) one.
     * Firstly only care about the first overlapped one, then take the second third .... overlapping as reference to run everything again
     * **/
    std::cout << " line 106   allPairedSimilarFragments[0].size() " << allPairedSimilarFragments[0].size() << std::endl;
    for( i=0; i<allPairedSimilarFragments[0].size(); ++i ){
        std::cout << " line 108  " << i << std::endl;
        //firstly set all the if overlap as false begin
        for ( j=1; j<allPairedSimilarFragments.size(); ++j ){
            for (k=0; k< allPairedSimilarFragments[j].size(); ++k){
                overlapped[j-1][k] = false;
            }
        } // set all the if overlap as false end

        for ( j=1; j<allPairedSimilarFragments.size(); ++j ){
            for ( k=0; k<allPairedSimilarFragments[j].size(); ++k ){
                if( pairedSimilarFragmentsOverlap(allPairedSimilarFragments[0][i], allPairedSimilarFragments[j][k] ) ) {
                    overlapped[j-1][k]=true;
                }
            }
        }
        // design a recursion algorithm for this purpose
        std::cout << " line 124  " << i << std::endl;
        std::vector< std::vector<int>> allCombinations = generateAllCombinations(overlapped, allPairedSimilarFragments, 0);  // ### X (n-2)
        std::cout << " line 126  " << i << std::endl;
        //overlapped has the same size with allPairedSimilarFragments
        for( j=0; j<allCombinations.size(); ++j ){
            int refStart = allPairedSimilarFragments[0][i].getStart1();
            int refEnd = allPairedSimilarFragments[0][i].getEnd1();
            bool allMatched = true; //every sequence ever match this similar fragments
            std::set<std::string> species;
            std::set<int> badIndex;
            for( k=0; k<allCombinations[j].size(); ++k ){
                if( allCombinations[j][k] > -1 ){ // there is real overlap
                    int oldRefStart = refStart;
                    int oldRefEnd = refEnd;
                    if( refStart < allPairedSimilarFragments[k+1][allCombinations[j][k]].getStart1() ){
                        refStart = allPairedSimilarFragments[k+1][allCombinations[j][k]].getStart1();
                    }
                    if( refEnd > allPairedSimilarFragments[k+1][allCombinations[j][k]].getEnd1()  ){
                        refEnd = allPairedSimilarFragments[k+1][allCombinations[j][k]].getEnd1();
                    }
                    if( (refEnd <= refStart) || (refEnd - refStart + 1) < mini_cns_size ){
                        refStart = oldRefStart;
                        refEnd = oldRefEnd;
                        badIndex.insert(k);
                    }else{
                        std::string spe = seqNames[k+2];
                        spe = spe.substr(18,5);
                        species.insert(spe);
                    }
                }else{
                    allMatched = false;
                }
            }
            //allPairedSimilarFragments[0] never outputed
            if ( /*allMatched*/ species.size()>=minimumNumberOfSpecies && (refEnd > refStart) && ((refEnd - refStart + 1) >= mini_cns_size) ) {
                std::vector<std::string> alignmentNames;
                std::vector<std::string> alignmentSeqs;

                alignmentNames.push_back(seqNames[0] + ":" + std::to_string(refStart) + "-" + std::to_string(refEnd)); //the reference one
                {
                    int start2 = allPairedSimilarFragments[0][i].getStart2();
                    int end2 = allPairedSimilarFragments[0][i].getEnd2();

                    int positionReference = allPairedSimilarFragments[0][i].getStart1();
                    --positionReference;

                    std::string alignSeq;
                    std::string refSeq;
                    int numberOfAltChar = 0;
                    bool alignRegionStart = false;
                    bool alignRegionEnd = true;
                    for ( int m=0; m<allPairedSimilarFragments[0][i].getAlignment1().length(); ++m ){
                        if( allPairedSimilarFragments[0][i].getAlignment1()[m] != '-' ){
                            ++positionReference;
                        }
                        if( !alignRegionStart ){
                            if( positionReference>=refStart ){
                                alignRegionStart =true;
                                start2 += numberOfAltChar;
                                end2 = start2;
                            }
                        }
                        if( alignRegionStart && alignRegionEnd ){
                            alignSeq += allPairedSimilarFragments[0][i].getAlignment2()[m];
                            refSeq += allPairedSimilarFragments[0][i].getAlignment1()[m];
                        }
                        if( allPairedSimilarFragments[0][i].getAlignment2()[m] != '-' ){
                            if( alignRegionStart && alignRegionEnd ){
                                ++end2;
                            }
                            ++numberOfAltChar;
                        }
                        if( alignRegionEnd ){
                            if( positionReference>=refEnd ){
                                alignRegionEnd = false;
                            }
                        }
                    }
                    --end2;

                    alignmentNames.push_back(seqNames[1] + ":" + std::to_string(start2) + "-" + std::to_string(end2));

                    alignmentSeqs.push_back(refSeq);
                    alignmentSeqs.push_back(alignSeq);
                }

                for (k = 0; k < allCombinations[j].size(); ++k) {
                    if( allCombinations[j][k] > -1 && badIndex.find(k)==badIndex.end() ){
                        int start2 = allPairedSimilarFragments[k+1][allCombinations[j][k]].getStart2();
                        int end2 = allPairedSimilarFragments[k+1][allCombinations[j][k]].getEnd2();

                        int positionReference = allPairedSimilarFragments[k+1][allCombinations[j][k]].getStart1();
                        --positionReference;

                        std::string alignSeq;
                        std::string refSeq;
                        int numberOfAltChar = 0;
                        bool alignRegionStart = false;
                        bool alignRegionEnd = true;
                        for ( int m=0; m<allPairedSimilarFragments[k+1][allCombinations[j][k]].getAlignment1().length(); ++m ){
                            if( allPairedSimilarFragments[k+1][allCombinations[j][k]].getAlignment1()[m] != '-' ){
                                ++positionReference;
                            }
                            if( !alignRegionStart ){
                                if( positionReference>=refStart ){
                                    alignRegionStart =true;
                                    start2 += numberOfAltChar;
                                    end2 = start2;
                                }
                            }
                            if( alignRegionStart && alignRegionEnd ){
                                alignSeq += allPairedSimilarFragments[k+1][allCombinations[j][k]].getAlignment2()[m];
                                refSeq += allPairedSimilarFragments[k+1][allCombinations[j][k]].getAlignment1()[m];
                            }
                            if( allPairedSimilarFragments[k+1][allCombinations[j][k]].getAlignment2()[m] != '-' ){
                                if( alignRegionStart && alignRegionEnd ){
                                    ++end2;
                                }
                                ++numberOfAltChar;
                            }
                            if( alignRegionEnd ){
                                if( positionReference>=refEnd ){
                                    alignRegionEnd = false;
                                }
                            }
                        }
                        --end2;

                        alignmentNames.push_back(seqNames[k+2] + ":" + std::to_string(start2) + "-" + std::to_string(end2));
                        {
                            int lengthq = refSeq.length();
                            for ( int o=0; o<lengthq; ++o ){
                                if( alignmentSeqs[0].length() > o ){
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
                        alignmentSeqs.push_back(alignSeq);
                    }
                }

                for ( int o=0; o<alignmentNames.size(); ++o ){
                    std::cout << ">" << alignmentNames[o]<< std::endl << alignmentSeqs[o] << std::endl;
                    ofile << ">" << alignmentNames[o]<< std::endl << alignmentSeqs[o] << std::endl;
                }
                ofile << std::endl << std::endl;
                std::cout << std::endl << std::endl;
            }
        }
    }
    ofile.close();
}
