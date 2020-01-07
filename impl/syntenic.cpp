//
// Created by Baoxing song on 2019-01-04.
//

#include "syntenic.h"


/**
 * This is a unweighted longest increasing subsequence that does not consider node long and node distance,
 * Try to implemented other method in the future
 *
 *
 * In the Mummer paper they said they are using a weighted longest Increasing Subsequence,
 * while I have not find more details to figure out how to implemented "weighted" longest Increasing Subsequence yet.
 * */

std::vector<PairedSimilarFragment> longestIncreasingSubsequence ( std::vector<PairedSimilarFragment> & pairedSimilarFragments){
    int maxLength = 1, bestEnd = 0;
    int DP[pairedSimilarFragments.size()];
    int prev[pairedSimilarFragments.size()];
    DP[0] = 1;
    prev[0] = -1;

    for (int i = 1; i < pairedSimilarFragments.size(); i++) {
        DP[i] = 1;
        prev[i] = -1;

        for (int j = i - 1; j >= 0; j--)
            if (DP[j] + 1 > DP[i] && pairedSimilarFragments[j].getStart2() < pairedSimilarFragments[i].getStart2()){
                DP[i] = DP[j] + 1;
                prev[i] = j;
            }
        if (DP[i] > maxLength){
            bestEnd = i;
            maxLength = DP[i];
        }
    }
    std::vector<PairedSimilarFragment> sorted_pairedSimilarFragments;
    int i=bestEnd;
    sorted_pairedSimilarFragments.push_back(pairedSimilarFragments[i]);
    int j = prev[i];
    while( j>=0 ){
        sorted_pairedSimilarFragments.push_back(pairedSimilarFragments[j]);
        j=prev[j];
    }

    // todo think about how to merge records
    std::vector<PairedSimilarFragment> filtered_sorted_pairedSimilarFragments;
    for( i=sorted_pairedSimilarFragments.size()-1; i>=0; --i ){
        filtered_sorted_pairedSimilarFragments.push_back(sorted_pairedSimilarFragments[i]);
    }
    return filtered_sorted_pairedSimilarFragments;
}

/**
 * This is the method described in the LAGAN/MLAGAN paper
 * But I did not care about the position overlap between neighboring records
 * */

std::vector<PairedSimilarFragment> longestIncreasingSubsequenceLAGAN_v0 ( std::vector<PairedSimilarFragment> & pairedSimilarFragments){
    // then for the seed-to-chain should check the overlap of pairedSimilarFragment
    int32_t maxSore = pairedSimilarFragments[0].getLength(), bestEnd = 0;
    int32_t DP[pairedSimilarFragments.size()];
    int prev[pairedSimilarFragments.size()];
    DP[0] = pairedSimilarFragments[0].getLength();
    prev[0] = -1;
    for (int i = 1; i < pairedSimilarFragments.size(); ++i) {
        DP[i] = pairedSimilarFragments[i].getLength();
        prev[i] = -1;
        for (int j = i - 1; j >= 0; --j){
            if (DP[j] + pairedSimilarFragments[i].getLength() > DP[i] && pairedSimilarFragments[j].getStart2() < pairedSimilarFragments[i].getStart2() && pairedSimilarFragments[j].getEnd2() < pairedSimilarFragments[i].getEnd2() && pairedSimilarFragments[j].getEnd1() < pairedSimilarFragments[i].getEnd1() ){
                DP[i] = DP[j] + pairedSimilarFragments[i].getLength();
                prev[i] = j;
            }
        }
        if (DP[i] > maxSore){
            bestEnd = i;
            maxSore = DP[i];
        }
    }

    std::vector<PairedSimilarFragment> sorted_pairedSimilarFragments;
    int i=bestEnd;
    sorted_pairedSimilarFragments.push_back(pairedSimilarFragments[i]);
    int j = prev[i];
    while( j>=0 ){
        sorted_pairedSimilarFragments.push_back(pairedSimilarFragments[j]);
        j=prev[j];
    }

    std::vector<PairedSimilarFragment> filtered_sorted_pairedSimilarFragments;
    for( i=sorted_pairedSimilarFragments.size()-1; i>=0; --i ){
        filtered_sorted_pairedSimilarFragments.push_back(sorted_pairedSimilarFragments[i]);
    }
    return filtered_sorted_pairedSimilarFragments;
}

std::vector<PairedSimilarFragment> longestIncreasingSubsequenceLAGAN ( std::vector<PairedSimilarFragment> & pairedSimilarFragments){
    // then for the seed-to-chain should check the overlap of pairedSimilarFragment
    int32_t maxSore = pow(pairedSimilarFragments[0].getScore(), 2), bestEnd = 0;
    int32_t DP[pairedSimilarFragments.size()];
    int prev[pairedSimilarFragments.size()];
    DP[0] = pow(pairedSimilarFragments[0].getScore(), 2);
    prev[0] = -1;
    for (int i = 1; i < pairedSimilarFragments.size(); ++i) {
        DP[i] = pow(pairedSimilarFragments[i].getScore(), 2);
        prev[i] = -1;
        for (int j = i - 1; j >= 0; --j){
            if (DP[j] + pow(pairedSimilarFragments[i].getScore(), 2) > DP[i] &&
                    pairedSimilarFragments[j].getStart2() < pairedSimilarFragments[i].getStart2() &&
                    (pairedSimilarFragments[j].getEnd2() < pairedSimilarFragments[i].getEnd2() ||
                    pairedSimilarFragments[j].getEnd1() < pairedSimilarFragments[i].getEnd1() )){
                DP[i] = DP[j] + pow(pairedSimilarFragments[i].getScore(), 2);
                prev[i] = j;
            }
        }
        if (DP[i] > maxSore){
            bestEnd = i;
            maxSore = DP[i];
        }
    }

    std::vector<PairedSimilarFragment> sorted_pairedSimilarFragments;
    int i=bestEnd;
    sorted_pairedSimilarFragments.push_back(pairedSimilarFragments[i]);
    int j = prev[i];
    while( j>=0 ){
        sorted_pairedSimilarFragments.push_back(pairedSimilarFragments[j]);
        j=prev[j];
    }

    std::vector<PairedSimilarFragment> filtered_sorted_pairedSimilarFragments;
    for( i=sorted_pairedSimilarFragments.size()-1; i>=0; --i ){
        filtered_sorted_pairedSimilarFragments.push_back(sorted_pairedSimilarFragments[i]);
    }
    return filtered_sorted_pairedSimilarFragments;
}

std::vector<PairedSimilarFragment> syntenic ( std::vector<PairedSimilarFragment> & pairedSimilarFragments){
    if( pairedSimilarFragments.size() == 0 ){
        return pairedSimilarFragments;
    }
    std::sort(pairedSimilarFragments.begin(), pairedSimilarFragments.end(), [](PairedSimilarFragment a, PairedSimilarFragment b) {
        return a.getStart1() < b.getStart1();
    });
    return longestIncreasingSubsequenceLAGAN(pairedSimilarFragments);
}
