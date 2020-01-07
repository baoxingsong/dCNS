//
// Created by Baoxing song on 2019-01-02.
//

/*
 * This one could only be used to get the maximum score and maximum score position
 * It could not be used to get the sequence alignment, since score matrix is not enough to trace back
 * **/


#include "sequenceAlignment.h"

/**
 * this is a standarded smith waterman implementation
 * should be used as benchmark to test other implementations
 * reverseAlignment is used to align the reverse sequence for CNS detection
 * if it is true, the algorithm try hardly to start the alignment from the first base pair for both query and target sequence
 * */

std::vector<uint32_t> SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                   bool reverseAlignment, bool returnCigar){

    int32_t ** M = new int32_t *[length1 + 1];
    int32_t ** E = new int32_t *[length1 + 1];
    int32_t ** F = new int32_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i) {
        M[i] = new int32_t [length2 + 1];
        E[i] = new int32_t [length2 + 1];
        F[i] = new int32_t [length2 + 1];
        std::fill_n(M[i], length2+1, 0);
        std::fill_n(E[i], length2+1, 0);
        std::fill_n(F[i], length2+1, 0);
    }

    endPosition1=0;
    endPosition2=0;
    maxScore = 0;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            E[i][j] = (_open_gap_penalty + M[i][j-1]) > (_extend_gap_penalty + E[i][j-1]) ? (_open_gap_penalty + M[i][j-1]) : (_extend_gap_penalty + E[i][j-1]);
//            E[i][j] = E[i][j] > (_open_gap_penalty + F[i][j-1]) ? E[i][j] : (_open_gap_penalty + F[i][j-1]);

            F[i][j] = (_open_gap_penalty + M[i-1][j]) > (_extend_gap_penalty + F[i-1][j]) ?( _open_gap_penalty + M[i-1][j]) : (_extend_gap_penalty + F[i-1][j]);
//            F[i][j] = F[i][j] > (_open_gap_penalty + E[i-1][j]) ? F[i][j] : (_open_gap_penalty + E[i-1][j]);

            M[i][j] = (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) > E[i][j] ? (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) : E[i][j];
            M[i][j] = M[i][j] > F[i][j] ? M[i][j] : F[i][j];
            M[i][j] = M[i][j] > 0 ? M[i][j] : 0;
            F[i][j] = F[i][j] > 0 ? F[i][j] : 0;
            E[i][j] = E[i][j] > 0 ? E[i][j] : 0;

            if( M[i][j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                // >= will omit the first similar fragments
                maxScore = M[i][j];
                endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2=j;
            }
        }
    }
    //std::cout << "maxScore " << maxScore << std::endl;
    std::vector<uint32_t> cigar;
    if(returnCigar ){ // cigar is only needed for reverse alignment
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int ii = endPosition1;
        int jj = endPosition2;
        while (ii>0 && jj>0) {
            if (ii > 0 && jj > 0 &&
                M[ii][jj] == M[ii - 1][jj - 1] + m.getScore(seq1[ii - 1], seq2[jj - 1]) ) {
                --ii;
                --jj;
                op = 0;
            } else if (jj > 0 && M[ii][jj] == E[ii][jj]) {
                --jj;
                op = 1;
            } else if (ii > 0 && M[ii][jj] == F[ii][jj]) {
                --ii;
                op = 2;
            } else if(M[ii][jj]==0){
                break;
            } else {
                std::cout << "M[ii][jj] " << M[ii][jj] << " E[ii][jj] " << E[ii][jj] << " F[ii][jj] " << F[ii][jj] << std::endl;
                std::cout << "there is something wrong with smith-waterman algorithm in line 114" << std::endl;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        if (!reverseAlignment){
            std::reverse(cigar.begin(),cigar.end());
        }
        // trace back end
    }

    for (i = 0; i <= length1; ++i) {
        delete[] M[i];
        delete[] E[i];
        delete[] F[i];
    }
    delete[] M;
    delete[] E;
    delete[] F;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}


// if you only want the maximum score and the position then use this one
void SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                    const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                                    int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & mi, bool positions){
    int8_t **m = mi.getScore();
    int32_t * M1 = new int32_t [length2 + 1];
    int32_t * M2 = new int32_t [length2 + 1];
    int32_t * F = new int32_t [length2 + 1];
    std::fill_n(M1, length2+1, 0);
    std::fill_n(F, length2+1, 0);
    std::fill_n(M2, length2+1, 0);

    int32_t e;
    int32_t * t;

    int32_t i, j;
    maxScore = 0;
    if (positions){
        endPosition1=0;
        endPosition2=0;
        for ( i=1; i<=length1; ++i ){
            e=0;
            for ( j=1; j<=length2; ++j ){
                e = (_open_gap_penalty + M2[j-1]) > (_extend_gap_penalty + e) ? (_open_gap_penalty + M2[j-1]) : (_extend_gap_penalty + e);

                F[j] = (_open_gap_penalty + M1[j]) >  (_extend_gap_penalty + F[j]) ?(_open_gap_penalty + M1[j]) :  (_extend_gap_penalty + F[j]);

                M2[j] = (m[seq1[i-1]][seq2[j-1]] + M1[j-1]) > e ? (m[seq1[i-1]][seq2[j-1]] + M1[j-1]) : e;
                M2[j] = M2[j] > F[j] ? M2[j] : F[j];

                M2[j] = M2[j] > 0 ? M2[j] : 0;
                F[j] = F[j] > 0 ? F[j] : 0;
                e = e > 0 ? e : 0;

                if( M2[j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2=j;
                }
            }
            t = M1;
            M1 = M2;
            M2 = t;
        }
        --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
        --endPosition2;
    }else{
        for ( i=1; i<=length1; ++i ){
            e=0;
            for ( j=1; j<=length2; ++j ){
                e = _open_gap_penalty + M2[j-1] > _extend_gap_penalty + e ? _open_gap_penalty + M2[j-1] : _extend_gap_penalty + e;
                F[j] = (_open_gap_penalty + M1[j]) >  (_extend_gap_penalty + F[j]) ?(_open_gap_penalty + M1[j]) :  (_extend_gap_penalty + F[j]);

                M2[j] = (m[seq1[i-1]][seq2[j-1]] + M1[j-1]) > e ? (m[seq1[i-1]][seq2[j-1]] + M1[j-1]) : e;
                M2[j] = M2[j] > F[j] ? M2[j] : F[j];

                M2[j] = M2[j] > 0 ? M2[j] : 0;
                F[j] = F[j] > 0 ? F[j] : 0;
                e = e > 0 ? e : 0;
                maxScore = maxScore > M2[j] ? maxScore : M2[j];
            }
            t = M1;
            M1 = M2;
            M2 = t;
        }
    }
    delete[] M1;
    delete[] F;
    delete[] M2;
}

// this function does not work correctly
void ssw_int8(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                                const int32_t &length2, const int8_t &_open_gap_penalty, const int8_t &_extend_gap_penalty,
                                                int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m, bool positions){

    // Put the largest number of the 16 numbers in vm into m.
#define max16(m, vm) (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 8)); \
					  (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 4)); \
					  (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 2)); \
					  (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 1)); \
					  (m) = _mm_extract_epi16((vm), 0)


    uint8_t bias = 0;
    size_t n = 5;  //ATCGN
    size_t segLen = ( length2 + 15)/16;

    // vProfile begain
    int32_t i, j, ai, k;

    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int8_t* t = (int8_t*)vProfile;

    for( ai=0; ai<n; ++ai ){
        for( i=0; i<segLen; ++i ){
            j = i;
            for( k=0; k<16; ++k ){
                *t++ = j>= length2 ? bias : m.getScore(ai, seq2[j-1]);
                j += segLen;
            }
        }
    }

    maxScore = 0;
    int32_t end_read = length2 - 1;
    int32_t end_ref = -1; /* 0_based best alignment ending point; Initialized as isn't aligned -1. */

    /* array to record the largest score of each reference position */
    uint8_t* maxColumn = (uint8_t*) calloc(length1, 1);

    /* array to record the alignment read ending position of the largest score of each reference position */
    int32_t* end_read_column = (int32_t*) calloc(length1, sizeof(int32_t));

    __m128i vZero = _mm_set1_epi32(0);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));

//    __m128i* vP;
    __m128i vGapO = _mm_set1_epi8(-_open_gap_penalty);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi8(-_extend_gap_penalty);

    __m128i vBias = _mm_set1_epi8(bias);

    __m128i vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxMark = vZero; /* Trace the highest score till the previous column. */
    __m128i vTemp;

    for (i = 0; i<length1; ++i) {//important

        int32_t cmp;
        __m128i e, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0.
							   Any errors to vH values will be corrected in the Lazy_F loop.
							 */

        __m128i vH = pvHStore[segLen - 1];
        vH = _mm_slli_si128 (vH, 1); /* Shift the 128-bit value in vH left by 1 byte. */
        const __m128i* vP = vProfile + seq1[i]*(segLen); /* Right part of the vProfile */

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */

        for (j = 0; j < segLen; ++j) {
            vH = _mm_adds_epu8(vH, _mm_load_si128(vP + j));
            vH = _mm_subs_epu8(vH, vBias); /* vH will be always > 0 */

            /* Get max from vH, vE and vF. */
            e = _mm_load_si128(pvE + j);
            vH = _mm_max_epu8(vH, e);
            vH = _mm_max_epu8(vH, vF);
            vMaxColumn = _mm_max_epu8(vMaxColumn, vH);

            /* Save vH values. */
            _mm_store_si128(pvHStore + j, vH);

            /* Update vE value. */
            vH = _mm_subs_epu8(vH, vGapO); /* saturation arithmetic, result >= 0 */
            e = _mm_subs_epu8(e, vGapE);
            e = _mm_max_epu8(e, vH);
            _mm_store_si128(pvE + j, e);

            /* Update vF value. */
            vF = _mm_subs_epu8(vF, vGapE);
            vF = _mm_max_epu8(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + j);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
        for (k = 0; (k < 16); ++k) {
            vF = _mm_slli_si128 (vF, 1);
            for (j = 0; (j < segLen); ++j) {
                vH = _mm_load_si128(pvHStore + j);
                vH = _mm_max_epu8(vH, vF);
                vMaxColumn = _mm_max_epu8(vMaxColumn, vH);	// newly added line
                _mm_store_si128(pvHStore + j, vH);
                vH = _mm_subs_epu8(vH, vGapO);
                vF = _mm_subs_epu8(vF, vGapE);
                if (( _mm_movemask_epi8(_mm_cmpgt_epi8(vF, vH)))) goto end; //todo think about here
            }
        }
        end:

        vMaxScore = _mm_max_epu8(vMaxScore, vMaxColumn);
        vTemp = _mm_cmpeq_epi8(vMaxMark, vMaxScore);
        cmp = _mm_movemask_epi8(vTemp);
        if (cmp != 0xffff) {
            uint8_t temp;
            vMaxMark = vMaxScore;
            max16(temp, vMaxScore);
            vMaxScore = vMaxMark;

            if ((temp > maxScore)) {
//                std::cout << maxScore << std::endl;
                maxScore = temp;
//                std::cout << maxScore << std::endl;
                if (maxScore + bias >= 255) break;	//overflow
                end_ref = i;

                /* Store the column with the highest alignment score in order to trace the alignment ending position on read. */
                for (j = 0; (j < segLen); ++j) pvHmax[j] = pvHStore[j];
            }
        }

        /* Record the max score of current column. */
        max16(maxColumn[i], vMaxColumn);
        //if (maxColumn[i] == 0) break; todo
    }
    if( positions ){


        /* Trace the alignment ending position on read. */
        uint8_t *t = (uint8_t*)pvHmax;
        int32_t column_len = segLen * 16;
        for (i = 0; (i < column_len); ++i, ++t) {
            int32_t temp;
            if (*t == maxScore) {
                temp = i / 16 + i % 16 * segLen;
                if (temp < end_read) end_read = temp;
            }
        }

        endPosition1 = end_ref;
        endPosition2 = end_read;

    }
    free(pvHmax);
    free(pvE);
    free(pvHLoad);
    free(pvHStore);

    free(maxColumn);
    free(end_read_column);

    free(vProfile);
}





// this is a 2-piece affine gap cost smithwaterman method which return cigar
std::vector<uint32_t> SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                    const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                                    const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                                    int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                                    bool reverseAlignment, bool returnCigar){

    int32_t ** M = new int32_t *[length1 + 1];
    int32_t ** E1 = new int32_t *[length1 + 1];
    int32_t ** E2 = new int32_t *[length1 + 1];
    int32_t ** F1 = new int32_t *[length1 + 1];
    int32_t ** F2 = new int32_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i) {
        M[i] = new int32_t [length2 + 1];
        E1[i] = new int32_t [length2 + 1];
        E2[i] = new int32_t [length2 + 1];
        F1[i] = new int32_t [length2 + 1];
        F2[i] = new int32_t [length2 + 1];
        std::fill_n(M[i], length2+1, 0);
        std::fill_n(E1[i], length2+1, 0);
        std::fill_n(E2[i], length2+1, 0);
        std::fill_n(F1[i], length2+1, 0);
        std::fill_n(F2[i], length2+1, 0);
    }

    endPosition1=0;
    endPosition2=0;
    maxScore = 0;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            E1[i][j] = (_open_gap_penalty1 + M[i][j-1]) > (_extend_gap_penalty1 + E1[i][j-1]) ? (_open_gap_penalty1 + M[i][j-1]) : (_extend_gap_penalty1 + E1[i][j-1]);
            E2[i][j] = (_open_gap_penalty2 + M[i][j-1]) > (_extend_gap_penalty2 + E2[i][j-1]) ? (_open_gap_penalty2 + M[i][j-1]) : (_extend_gap_penalty2 + E2[i][j-1]);

            F1[i][j] = (_open_gap_penalty1 + M[i-1][j]) > (_extend_gap_penalty1 + F1[i-1][j]) ? (_open_gap_penalty1 + M[i-1][j]) : (_extend_gap_penalty1 + F1[i-1][j]);
            F2[i][j] = (_open_gap_penalty2 + M[i-1][j]) >  (_extend_gap_penalty2 + F2[i-1][j]) ?(_open_gap_penalty2 + M[i-1][j]) : (_extend_gap_penalty2 + F2[i-1][j]);

            M[i][j] = (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) > E1[i][j] ? (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) : E1[i][j];
            M[i][j] = M[i][j] > F1[i][j] ? M[i][j] : F1[i][j];
            M[i][j] = M[i][j] > F2[i][j] ? M[i][j] : F2[i][j];
            M[i][j] = M[i][j] > E2[i][j] ? M[i][j] : E2[i][j];

            M[i][j] = M[i][j] > 0 ? M[i][j] : 0;
            F1[i][j] = F1[i][j] > 0 ? F1[i][j] : 0;
            E1[i][j] = E1[i][j] > 0 ? E1[i][j] : 0;
            F2[i][j] = F2[i][j] > 0 ? F2[i][j] : 0;
            E2[i][j] = E2[i][j] > 0 ? E2[i][j] : 0;

            if( M[i][j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                maxScore = M[i][j];
                endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2=j;
            }
        }
    }
    //std::cout << "maxScore " << maxScore << std::endl;
    std::vector<uint32_t> cigar;

    if( returnCigar ){ // cigar is only needed for reverse alignment
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int ii = endPosition1;
        int jj = endPosition2;
        while (ii>0 && jj>0) {
            if (ii > 0 && jj > 0 &&
                M[ii][jj] == M[ii - 1][jj - 1] + m.getScore(seq1[ii - 1], seq2[jj - 1]) ) {
                --ii;
                --jj;
                op = 0;
            } else if (jj > 0 && M[ii][jj] == E1[ii][jj]) {
                --jj;
                op = 1;
            } else if (ii > 0 && M[ii][jj] == F1[ii][jj]) {
                --ii;
                op = 2;
            } else if (ii > 0 && M[ii][jj] == F2[ii][jj]) {
                --ii;
                op = 2;
            } else if (ii > 0 && M[ii][jj] == E2[ii][jj]) {
                --jj;
                op = 1;
            } else if(M[ii][jj]==0){
                //std::cerr << "there is something wrong with smith-waterman algorithm in line 110" << std::endl;
                // here should never run, there is some problem with the code
                break;
            } else {
//                std::cout << "M[ii][jj] " << M[ii][jj] << " E[ii][jj] " << E[ii][jj] << " F[ii][jj] " << F[ii][jj] << std::endl;
                std::cout << "there is something wrong with smith-waterman algorithm in line 114" << std::endl;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        if (!reverseAlignment){
            std::reverse(cigar.begin(),cigar.end());
        }
    }

    for (i = 0; i <= length1; ++i) {
        delete[] M[i];
        delete[] E1[i];
        delete[] F1[i];
        delete[] E2[i];
        delete[] F2[i];
    }
    delete[] M;
    delete[] E1;
    delete[] F1;
    delete[] E2;
    delete[] F2;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}



void SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                   const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m){

    int32_t * M1 = new int32_t [length2 + 1];
    int32_t * F11 = new int32_t [length2 + 1];
    int32_t * F12 = new int32_t [length2 + 1];

    int32_t * M2 = new int32_t [length2 + 1];
    int32_t * F21 = new int32_t [length2 + 1];
    int32_t * F22 = new int32_t [length2 + 1];

    int32_t * t;

    std::fill_n(M1, length2+1, 0);
    std::fill_n(F11, length2+1, 0);
    std::fill_n(F12, length2+1, 0);

    std::fill_n(M2, length2+1, 0);
    std::fill_n(F21, length2+1, 0);
    std::fill_n(F22, length2+1, 0);

    endPosition1=0;
    endPosition2=0;
    int32_t i, j, e1, e2;
    maxScore = 0;

    for ( i=1; i<=length1; ++i ){
        e1=0, e2=0;
        for ( j=1; j<=length2; ++j ){
            e1 = (_open_gap_penalty1 + M2[j-1]) > (_extend_gap_penalty1 + e1) ? (_open_gap_penalty1 + M2[j-1]) : (_extend_gap_penalty1 + e1);
            F21[j] = (_open_gap_penalty1 + M1[j]) >  (_extend_gap_penalty1 + F11[j]) ?(_open_gap_penalty1 + M1[j]) :  (_extend_gap_penalty1 + F11[j]);

            e2 = (_open_gap_penalty2 + M2[j-1]) > (_extend_gap_penalty2 + e2) ? (_open_gap_penalty2 + M2[j-1]) : (_extend_gap_penalty2 + e2);
            F22[j] = (_open_gap_penalty2 + M1[j]) >  (_extend_gap_penalty2 + F12[j]) ?(_open_gap_penalty2 + M1[j]) :  (_extend_gap_penalty2 + F12[j]);

            M2[j] = (m.getScore(seq1[i-1], seq2[j-1]) + M1[j-1]) > e1 ? (m.getScore(seq1[i-1], seq2[j-1]) + M1[j-1]) : e1;
            M2[j] = M2[j] > F21[j] ? M2[j] : F21[j];
            M2[j] = M2[j] > F22[j] ? M2[j] : F22[j];
            M2[j] = M2[j] > e2 ? M2[j] :e2;

            M2[j] = M2[j] > 0 ? M2[j] : 0;
            F21[j] = F21[j] > 0 ? F21[j] : 0;
            e1 = e1 > 0 ? e1 : 0;
            F22[j] = F22[j] > 0 ? F22[j] : 0;
            e2 = e2 > 0 ? e2 : 0;

            if( M2[j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                maxScore = M2[j];
                endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2=j;
            }
        }
        t = M1;
        M1 = M2;
        M2 = t;

        t = F11;
        F11 = F21;
        F21 = t;

        t = F12;
        F12 = F22;
        F22 = t;
    }

    delete[] M1;
    delete[] F11;
    delete[] F12;

    delete[] M2;
    delete[] F21;
    delete[] F22;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
}



// not well read yet, and not tested
std::vector<uint32_t> NeedlemanWunsch(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                    const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                                    const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                                    const Scorei & m,
                                    bool reverseAlignment, bool returnCigar){

    int32_t ** M = new int32_t *[length1 + 1];
    int32_t ** E1 = new int32_t *[length1 + 1];
    int32_t ** E2 = new int32_t *[length1 + 1];
    int32_t ** F1 = new int32_t *[length1 + 1];
    int32_t ** F2 = new int32_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i) {
        M[i] = new int32_t [length2 + 1];
        E1[i] = new int32_t [length2 + 1];
        E2[i] = new int32_t [length2 + 1];
        F1[i] = new int32_t [length2 + 1];
        F2[i] = new int32_t [length2 + 1];
        std::fill_n(M[i], length2+1, 0);
        std::fill_n(E1[i], length2+1, 0);
        std::fill_n(E2[i], length2+1, 0);
        std::fill_n(F1[i], length2+1, 0);
        std::fill_n(F2[i], length2+1, 0);
    }
    i=0;
    for(j=1; j<(length2+1); ++j){
        F1[i][j] = j * _extend_gap_penalty1;
        F2[i][j] = j * _extend_gap_penalty2;
    }
    j=0;
    for(i=1; i<(length1+1); ++i){
        E1[i][j] = i * _extend_gap_penalty1;
        E2[i][j] = i * _extend_gap_penalty2;
    }
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            E1[i][j] = (_open_gap_penalty1 + M[i][j-1]) > (_extend_gap_penalty1 + E1[i][j-1]) ? (_open_gap_penalty1 + M[i][j-1]) : (_extend_gap_penalty1 + E1[i][j-1]);

            E2[i][j] = (_open_gap_penalty2 + M[i][j-1]) > (_extend_gap_penalty2 + E2[i][j-1]) ? (_open_gap_penalty2 + M[i][j-1]) : (_extend_gap_penalty2 + E2[i][j-1]);

            F1[i][j] = (_open_gap_penalty1 + M[i-1][j]) > (_extend_gap_penalty1 + F1[i-1][j]) ? (_open_gap_penalty1 + M[i-1][j]) : (_extend_gap_penalty1 + F1[i-1][j]);

            F2[i][j] = (_open_gap_penalty2 + M[i-1][j]) >  (_extend_gap_penalty2 + F2[i-1][j]) ?(_open_gap_penalty2 + M[i-1][j]) : (_extend_gap_penalty2 + F2[i-1][j]);

            M[i][j] = (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) > E1[i][j] ? (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) : E1[i][j];
            M[i][j] = M[i][j] > F1[i][j] ? M[i][j] : F1[i][j];
            M[i][j] = M[i][j] > F2[i][j] ? M[i][j] : F2[i][j];
            M[i][j] = M[i][j] > E2[i][j] ? M[i][j] : E2[i][j];

        }
    }
    //std::cout << "maxScore " << maxScore << std::endl;
    std::vector<uint32_t> cigar;
    if( reverseAlignment && returnCigar ){ // cigar is only needed for reverse alignment
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int ii = length1;
        int jj = length2;
        while (ii>0 && jj>0) {
            if (ii > 0 && jj > 0 &&
                M[ii][jj] == M[ii - 1][jj - 1] + m.getScore(seq1[ii - 1], seq2[jj - 1])) {
                --ii;
                --jj;
                op = 0;
            } else if (jj > 0 && M[ii][jj] == E1[ii][jj]) {
                --jj;
                op = 1;
            } else if (ii > 0 && M[ii][jj] == F1[ii][jj]) {
                --ii;
                op = 2;
            } else if (ii > 0 && M[ii][jj] == F2[ii][jj]) {
                --ii;
                op = 2;
            } else if (ii > 0 && M[ii][jj] == E2[ii][jj]) {
                --jj;
                op = 1;
            } else if(M[ii][jj]==0){
                //std::cerr << "there is something wrong with smith-waterman algorithm in line 110" << std::endl;
                // here should never run, there is some problem with the code
                break;
            } else {
//                std::cout << "M[ii][jj] " << M[ii][jj] << " E[ii][jj] " << E[ii][jj] << " F[ii][jj] " << F[ii][jj] << std::endl;
                std::cout << "there is something wrong with smith-waterman algorithm in line 716" << std::endl;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
    }

    for (i = 0; i <= length1; ++i) {
        delete[] M[i];
        delete[] E1[i];
        delete[] F1[i];
        delete[] E2[i];
        delete[] F2[i];
    }
    delete[] M;
    delete[] E1;
    delete[] F1;
    delete[] E2;
    delete[] F2;
    return cigar;
}


std::vector<uint32_t> SemiGlobal(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                                 const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & mi,
                                 bool returnCigar, const int32_t & zdrop, const int32_t & w, uint8_t **T){
    int8_t ** m = mi.getScore();
    int32_t i=0, j=0, mscore=0;
//
//    for ( i=0; i<60; ++i ){
//        std::cout << static_cast<int16_t>(*(seq1+i));
//    }
//    std::cout << std::endl;
//    for ( i=0; i<60; ++i ){
//        std::cout << static_cast<int16_t>(*(seq2+i));
//    }
//    std::cout << std::endl;


    uint8_t d=0;
    if( returnCigar ) {
        for (i = 0; i < (length1 + 1); ++i) {
            std::fill_n(T[i], length2 + 1, 0);
        }
    }
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F1 = new int32_t [length2 + 1]; //F1 and F2 is for gapscore1 and gapscore2
    int32_t * F2 = new int32_t [length2 + 1];

    std::fill_n(M1, length2+1, -536870912);
    std::fill_n(M2, length2+1, -536870912);
    std::fill_n(F1, length2+1, -536870912);
    std::fill_n(F2, length2+1, -536870912);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = 0;
    endPosition2 = 0;

    maxScore = 0;
    int32_t e1, e2;
    for ( i=1; i<=length1; ++i ){
        int32_t thisMax = 0;
        int thisMaxj=0;

        int32_t start = i - w;
        if ( start < 1 ){
            start = 1;
        }
        int32_t end = i + w;
        if( end > length2 ){
            end = length2;
        }
/*
        e1 = -536870912;
        e2 = -536870912;

        F1[start] = -536870912;
        F1[end] = -536870912;
        F2[start] = -536870912;
        F2[end] = -536870912;
*/

        e1 = -536870912;
        e2 = -536870912;

        if( returnCigar ) {
            for (j = 1; j <= length2; ++j) {
                mscore = m[(*(seq1+(i - 1)))][(*(seq2+(j - 1)))] + M1[j - 1];

                d = mscore > F1[j] ? 0 : 1;
                M2[j] = mscore > F1[j] ? mscore : F1[j];

                d = M2[j] > e1 ? d : 2;
                M2[j] = M2[j] > e1 ? M2[j] : e1;

                d = M2[j] > F2[j] ? d : 3;
                M2[j] = M2[j] > F2[j] ? M2[j] : F2[j];

                d = M2[j] > e2 ? d : 4;
                M2[j] = M2[j] > e2 ? M2[j] : e2;



                if ( M2[j] > maxScore) { // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2 = j;
                }

                if (j>=endPosition2 && M2[j] > thisMax){
                    thisMax = M2[j];
                    thisMaxj = j;
                }

                int32_t h = M2[j] + _open_gap_penalty1;
                int32_t f = F1[j] + _extend_gap_penalty1;
                d |= f >= h ? 1 << 3 : 0;
                f = f >= h ? f : h;
                F1[j] = f;

                e1 += _extend_gap_penalty1;
                d |= e1 >= h ? 1 << 4 : 0;
                e1 = e1 >= h ? e1 : h;

                int32_t h2 = M2[j] + _open_gap_penalty2;
                int32_t f2 = F2[j] + _extend_gap_penalty2;
                d |= f2 >= h2 ? 1 << 5 : 0;
                f2 = f2 >= h2 ? f2 : h2;
                F2[j] = f2;

                e2 += _extend_gap_penalty2;
                d |= e2 >= h2 ? 1 << 6 : 0;
                e2 = e2 >= h2 ? e2 : h2;
                T[i][j] = d;
            }
        }else{
            for ( j=1; j<=length2; ++j ){
                mscore = m[(*(seq1+(i - 1)))][(*(seq2+(j - 1)))] + M1[j - 1];

                M2[j] = mscore > F1[j] ? mscore : F1[j];

                M2[j] = M2[j] > e1 ? M2[j] : e1;

                M2[j] = M2[j] > F2[j] ? M2[j] : F2[j];

                M2[j] = M2[j] > e2 ? M2[j] : e2;


                if( M2[j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2=j;
                }

                if (j>=endPosition2 && M2[j] > thisMax){
                    thisMax = M2[j];
                    thisMaxj = j;
                }

                int32_t h = M2[j] + _open_gap_penalty1;
                int32_t f = F1[j] + _extend_gap_penalty1;
                f  = f >= h? f    : h;
                F1[j] = f;

                e1 += _extend_gap_penalty1;
                e1  = e1 >= h? e1  : h;

                int32_t h2 = M2[j] + _open_gap_penalty2;
                int32_t f2 = F2[j] + _extend_gap_penalty2;
                f2 = f2 >= h2? f2 : h2;
                F2[j] = f2;

                e2 += _extend_gap_penalty2;
                e2 = e2 >= h2? e2 : h2;
            }
        }

        if( thisMaxj > endPosition2 && i > endPosition1 ){
            int32_t l = (i - endPosition1) - (thisMaxj-endPosition2);
            l = l>0 ? l : -l;
            if( maxScore-thisMax > zdrop + l * _extend_gap_penalty2) {
                break;
            }
        }
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
    }
//    std::cout << "line 861 endPosition1 " << endPosition1 << " endPosition2: " << endPosition2 << " maxScore: " << maxScore << std::endl;
    std::vector<uint32_t> cigar;
    if( returnCigar ){
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int ii = endPosition1;
        int jj = endPosition2;
        uint32_t tmp;
        int state=0;
        while (ii>0 && jj>0) {
            tmp = T[ii][jj];
            if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
            else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
            if (state == 0) state = tmp & 7;
            if (state == 0){
                op=0;
                --ii;
                --jj;
            }else if (state == 1 || state == 3){
                op =2;
                --ii;
            }else {
                op =1;
                --jj;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        std::reverse(cigar.begin(),cigar.end());
        // for debug purpose todo remove it
        if ( 0==ii && 0==jj ){

        }else{
            std::cout << "ii:" << ii << " jj:" << jj << std::endl;
        }
        // for debug end
    }

//    std::cout << " line 511 maxScore " << maxScore << std::endl;
    delete[] M1;
    delete[] M2;
    delete[] F2;
    delete[] F1;

    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
//    std::cout << "line 889 endPosition1 " << endPosition1 << " endPosition2: " << endPosition2 << " maxScore: " << maxScore << std::endl;

    return cigar;
}



void SemiGlobal_xextend(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                                 const int32_t & xdrop, const int32_t & w){
    int32_t i, j, mscore;
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F = new int32_t [length2 + 1]; //F1 and F2 is for gapscore1 and gapscore2

    std::fill_n(M1, length2+1, -536870912);
    std::fill_n(M2, length2+1, -536870912);
    std::fill_n(F, length2+1, -536870912);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = 0;
    endPosition2 = 0;

    maxScore = 0;
    int32_t e;
    for ( i=1; i<=length1; ++i ){
        int32_t thisMax = 0;
        int thisMaxj=0;
        int32_t start = i - w;
        if ( start < 1 ){
            start = 1;
        }
        int32_t end = i + w;
        if( end > length2 ){
            end = length2;
        }
        e = -536870912;
//        F[start] = -536870912;
//        F[end] = -536870912;
        for ( j=start; j<=end; ++j ){
            mscore = m.getScore(seq1[i-1], seq2[j-1]) + M1[j-1];
            M2[j] = mscore > F[j] ? mscore : F[j];
            M2[j] = M2[j] > e ? M2[j] : e;
            if( M2[j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                // >= will omit the first similar fragments
                maxScore = M2[j];
                endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2=j;
            }
            if (M2[j] > thisMax){
                thisMax = M2[j];
                thisMaxj = j;
            }
            int32_t h = M2[j] + _open_gap_penalty;
            int32_t f = F[j] + _extend_gap_penalty;
            f  = f >= h? f    : h;
            F[j] = f;
            e += _extend_gap_penalty;
            e  = e >= h? e  : h;
        }
        if( thisMaxj > endPosition2 && i > endPosition1 ){
            if( maxScore-thisMax > xdrop ) {
                break;
            }
        }
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
    }
    delete[] M1;
    delete[] M2;
    delete[] F;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
}




/*************************************************************************

a weighted dynamic programming sequence alignment method ZDP (Zebric dynamic programming)

The standard needleman-wunsch algorithm set up a score matrix and generate alignment by tracing back
Since we have different score weight for each base-pair, the tracing back step would query the weighted score again, which is time consuming.
Here when set up the score matrix, we setup a trace matrix at the same time, and we could trace back using the trace matrix
        (this is a kind of approximation, might could not get the identical alignment result comparing the the standard traceback approach,
while the standard traceback might not work for the weighted sequence alignment approach, since different weight could generate indistinct trace path)

************************************************************************/



// seqAchar is the reference sequence
// seqBChar is the query sequence

// lengthA the the length of the reference sequence (if the sequence is a java/cpp string, this variable is not necessary)
// lengthA the the length of the query sequence

// weight is a vector of weight of each reference sequence bp
// Score is a custom class, together with the 'weight' variable,  we could get difference match/mis-match/insertion/deletion scores for each base-pair

// A is the alignment result using CIGAR letters [MDI]


// this is the semiglobal weighted sequence alignment approach
std::vector<uint32_t> SemiGlobal(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int16_t * weights, Score & score,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2,
                                 bool returnCigar, const int32_t & z, const int32_t & w, uint8_t **T){
    int32_t i, j, mscore;
    uint8_t d;
    if( returnCigar ) {
        for (i = 0; i < (length1 + 1); ++i) {
            std::fill_n(T[i], length2 + 1, 0);
        }
    }
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F1 = new int32_t [length2 + 1]; //F1 and F2 is for gapscore1 and gapscore2
    int32_t * F2 = new int32_t [length2 + 1];

    std::fill_n(M1, length2+1, -536870912);
    std::fill_n(M2, length2+1, -536870912);
    std::fill_n(F1, length2+1, -536870912);
    std::fill_n(F2, length2+1, -536870912);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = 0;
    endPosition2 = 0;

    maxScore = 0;
    int32_t e1, e2, zdrop, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2;
    int32_t ** m;
    int32_t thisMax;
    int thisMaxj;
    int32_t start;
    int32_t end;
    for ( i=1; i<=length1; ++i ){
        thisMax = 0;
        thisMaxj=0;
        start = i - w > 1 ? i - w : 1;
        end = i + w < length2 ?  i + w : length2;
        e1 = -536870912;
        e2 = -536870912;
//        F1[start] = -536870912;
//        F1[end] = -536870912;
//        F2[start] = -536870912;
//        F2[end] = -536870912;

        _open_gap_penalty1 = score.getOpenGapPenalty1(weights[i-1]);
        _extend_gap_penalty1 = score.getExtendGapPenalty1(weights[i-1]);
        _open_gap_penalty2 = score.getOpenGapPenalty2(weights[i-1]);
        _extend_gap_penalty2 = score.getExtendGapPenalty2(weights[i-1]);
        m = score.getM(weights[i-1]);

        //zdrop = z * (weights[i-1]/2); // 1/2===1
        zdrop = score.getZdrop(weights[i-1]);

        if( returnCigar ) {
            for (j = start; j <= end; ++j) {
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];

                d = mscore > F1[j] ? 0 : 1;
                M2[j] = mscore > F1[j] ? mscore : F1[j];

                d = M2[j] > e1 ? d : 2;
                M2[j] = M2[j] > e1 ? M2[j] : e1;

                d = M2[j] > F2[j] ? d : 3;
                M2[j] = M2[j] > F2[j] ? M2[j] : F2[j];

                d = M2[j] > e2 ? d : 4;
                M2[j] = M2[j] > e2 ? M2[j] : e2;


                if ( M2[j] > maxScore) { // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2 = j;
                }

                if (j>=endPosition2 && M2[j] > thisMax){
                    thisMax = M2[j];
                    thisMaxj = j;
                }

                int32_t h = M2[j] + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F1[j] + score.getExtendGapPenalty1(weights[i]);
                d |= f >= h ? 1 << 3 : 0;
                f = f >= h ? f : h;
                F1[j] = f;

                h = M2[j] + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                d |= e1 >= h ? 1 << 4 : 0;
                e1 = e1 >= h ? e1 : h;

                int32_t h2 = M2[j] + score.getOpenGapPenalty2(weights[i]);
                int32_t f2 = F2[j] + score.getExtendGapPenalty2(weights[i]);
                d |= f2 >= h2 ? 1 << 5 : 0;
                f2 = f2 >= h2 ? f2 : h2;
                F2[j] = f2;

                h2 = M2[j] + _open_gap_penalty2;
                e2 += _extend_gap_penalty2;
                d |= e2 >= h2 ? 1 << 6 : 0;
                e2 = e2 >= h2 ? e2 : h2;
                T[i][j] = d; // z[i,j]
            }
        }else{
            for ( j=start; j<=end; ++j ){
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];

                M2[j] = mscore > F1[j] ? mscore : F1[j];

                M2[j] = M2[j] > e1 ? M2[j] : e1;

                M2[j] = M2[j] > F2[j] ? M2[j] : F2[j];

                M2[j] = M2[j] > e2 ? M2[j] : e2;

                if( M2[j] >0 && M2[j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2=j;
                }

                if (j>=endPosition2 && M2[j] > thisMax){
                    thisMax = M2[j];
                    thisMaxj = j;
                }

                int32_t h = M2[j] + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F1[j] + score.getExtendGapPenalty1(weights[i]);
                f  = f >= h? f    : h;
                F1[j] = f;

                h = M2[j] + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                e1  = e1 >= h? e1  : h;

                int32_t h2 = M2[j] + score.getOpenGapPenalty2(weights[i]);
                int32_t f2 = F2[j] + score.getExtendGapPenalty2(weights[i]);
                f2 = f2 >= h2? f2 : h2;
                F2[j] = f2;

                h2 = M2[j] + _open_gap_penalty2;
                e2 += _extend_gap_penalty2;
                e2 = e2 >= h2? e2 : h2;
            }
        }

        if( thisMaxj > endPosition2 && i > endPosition1 ){
            int32_t l = (i - endPosition1) - (thisMaxj-endPosition2);
            l = l>0 ? l : -l;
            if( maxScore-thisMax > zdrop + l * _extend_gap_penalty2) {
                break;
            }
        }
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
    }

    std::vector<uint32_t> cigar;

    if( returnCigar ){
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int ii = endPosition1;
        int jj = endPosition2;
        int tmp, state=0;

        while (ii>0 && jj>0) {
            tmp = T[ii][jj];
            if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
            else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
            if (state == 0) state = tmp & 7;
            if (state == 0){
                op=0;
                --ii;
                --jj;
            }else if (state == 1 || state == 3){
                op =2;
                --ii;
            }else {
                op =1;
                --jj;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        std::reverse(cigar.begin(),cigar.end());
    }

    delete[] M1;
    delete[] M2;
    delete[] F2;
    delete[] F1;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}



// this is the semiglobal weighted sequence alignment approach
std::vector<uint32_t> SemiGlobal_single_gap_penalty1(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int16_t * weights, Score & score,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2,
                                 bool returnCigar, const int32_t & z, const int32_t & w, uint8_t **T){
    int32_t i, j, mscore;
    uint8_t d;
    if( returnCigar ) {
        for (i = 0; i < (length1 + 1); ++i) {
            std::fill_n(T[i], length2 + 1, 0);
        }
    }
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F1 = new int32_t [length2 + 1]; //F1 and F2 is for gapscore1 and gapscore2
    int32_t * F2 = new int32_t [length2 + 1];

    std::fill_n(M1, length2+1, -536870912);
    std::fill_n(M2, length2+1, -536870912);
    std::fill_n(F1, length2+1, -536870912);
    std::fill_n(F2, length2+1, -536870912);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = 0;
    endPosition2 = 0;

    maxScore = 0;
    int32_t e1, e2, zdrop, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2;
    int32_t ** m;
    int32_t thisMax;
    int thisMaxj;
    int32_t start;
    int32_t end;
    for ( i=1; i<=length1; ++i ){
        thisMax = 0;
        thisMaxj = 0;
        start = i - w > 1 ? i - w : 1;
        end = i + w < length2 ?  i + w : length2;
        e1 = -536870912;
        e2 = -536870912;
//        F1[start] = -536870912;
//        F1[end] = -536870912;
//        F2[start] = -536870912;
//        F2[end] = -536870912;

        _open_gap_penalty1 = score.getOpenGapPenalty1(weights[i-1]);
        _extend_gap_penalty1 = score.getExtendGapPenalty1(weights[i-1]);
        _open_gap_penalty2 = score.getOpenGapPenalty1(weights[i-1]);
        _extend_gap_penalty2 = score.getExtendGapPenalty1(weights[i-1]);
        m = score.getM(weights[i-1]);

        //zdrop = z * (weights[i-1]/2); // 1/2===1
        zdrop = score.getZdrop(weights[i-1]);

        if( returnCigar ) {
            for (j = start; j <= end; ++j) {
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];

                d = mscore > F1[j] ? 0 : 1;
                M2[j] = mscore > F1[j] ? mscore : F1[j];

                d = M2[j] > e1 ? d : 2;
                M2[j] = M2[j] > e1 ? M2[j] : e1;

                d = M2[j] > F2[j] ? d : 3;
                M2[j] = M2[j] > F2[j] ? M2[j] : F2[j];

                d = M2[j] > e2 ? d : 4;
                M2[j] = M2[j] > e2 ? M2[j] : e2;


                if ( M2[j] > maxScore) { // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2 = j;
                }

                if (j>=endPosition2 && M2[j] > thisMax){
                    thisMax = M2[j];
                    thisMaxj = j;
                }

                int32_t h = M2[j] + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F1[j] + score.getExtendGapPenalty1(weights[i]);
                d |= f >= h ? 1 << 3 : 0;
                f = f >= h ? f : h;
                F1[j] = f;

                h = M2[j] + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                d |= e1 >= h ? 1 << 4 : 0;
                e1 = e1 >= h ? e1 : h;

                int32_t h2 = M2[j] + score.getOpenGapPenalty2(weights[i]);
                int32_t f2 = F2[j] + score.getExtendGapPenalty2(weights[i]);
                d |= f2 >= h2 ? 1 << 5 : 0;
                f2 = f2 >= h2 ? f2 : h2;
                F2[j] = f2;

                h2 = M2[j] + _open_gap_penalty2;
                e2 += _extend_gap_penalty2;
                d |= e2 >= h2 ? 1 << 6 : 0;
                e2 = e2 >= h2 ? e2 : h2;
                T[i][j] = d; // z[i,j]
            }
        }else{
            for ( j=start; j<=end; ++j ){
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j-1];

                M2[j] = mscore > F1[j] ? mscore : F1[j];

                M2[j] = M2[j] > e1 ? M2[j] : e1;

                M2[j] = M2[j] > F2[j] ? M2[j] : F2[j];

                M2[j] = M2[j] > e2 ? M2[j] : e2;

                if( M2[j] >0 && M2[j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2=j;
                }

                if (j>=endPosition2 && M2[j] > thisMax){
                    thisMax = M2[j];
                    thisMaxj = j;
                }

                int32_t h = M2[j] + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F1[j] + score.getExtendGapPenalty1(weights[i]);
                f  = f >= h? f    : h;
                F1[j] = f;

                h = M2[j] + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                e1  = e1 >= h? e1  : h;

                int32_t h2 = M2[j] + score.getOpenGapPenalty2(weights[i]);
                int32_t f2 = F2[j] + score.getExtendGapPenalty2(weights[i]);
                f2 = f2 >= h2? f2 : h2;
                F2[j] = f2;

                h2 = M2[j] + _open_gap_penalty2;
                e2 += _extend_gap_penalty2;
                e2 = e2 >= h2? e2 : h2;
            }
        }

        if( thisMaxj > endPosition2 && i > endPosition1 ){
            int32_t l = (i - endPosition1) - (thisMaxj-endPosition2);
            l = l>0 ? l : -l;
            if( maxScore-thisMax > zdrop + l * _extend_gap_penalty2) {
                break;
            }
        }
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
    }

    std::vector<uint32_t> cigar;

    if( returnCigar ){
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int ii = endPosition1;
        int jj = endPosition2;
        int tmp, state=0;

        while (ii>0 && jj>0) {
            tmp = T[ii][jj];
            if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
            else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
            if (state == 0) state = tmp & 7;
            if (state == 0){
                op=0;
                --ii;
                --jj;
            }else if (state == 1 || state == 3){
                op =2;
                --ii;
            }else {
                op =1;
                --jj;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        std::reverse(cigar.begin(),cigar.end());
    }

    delete[] M1;
    delete[] M2;
    delete[] F2;
    delete[] F1;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}



// this is the semiglobal weighted sequence alignment approach
std::vector<uint32_t> SemiGlobal_single_gap_penalty(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int16_t * weights, Score & score,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2,
                                 bool returnCigar, const int32_t & z, const int32_t & w, uint8_t **T){
    int32_t i, j, mscore;
    uint8_t d;
    if( returnCigar ) {
        for (i = 0; i < (length1 + 1); ++i) {
            std::fill_n(T[i], length2 + 1, 0);
        }
    }
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F = new int32_t [length2 + 1]; //F1 and F2 is for gapscore1 and gapscore2

    std::fill_n(M1, length2+1, -536870912);
    std::fill_n(M2, length2+1, -536870912);
    std::fill_n(F, length2+1, -536870912);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = 0;
    endPosition2 = 0;

    maxScore = 0;
    int32_t e1, zdrop, _open_gap_penalty1, _extend_gap_penalty1;
    int32_t ** m;

    for ( i=1; i<=length1; ++i ){
        int32_t thisMax = 0;
        int thisMaxj=0;
        int32_t start = i - w;
        if ( start < 1 ){
            start = 1;
        }
        int32_t end = i + w;
        if( end > length2 ){
            end = length2;
        }

        e1 = -536870912;

//        F[start] = -536870912;
//        F[end] = -536870912;

        _open_gap_penalty1 = score.getOpenGapPenalty1(weights[i-1]);
        _extend_gap_penalty1 = score.getExtendGapPenalty1(weights[i-1]);
        m = score.getM(weights[i-1]);

        //zdrop = z * (weights[i-1]/2); // 1/2===1
        zdrop = score.getZdrop(weights[i-1]);

        if( returnCigar ) {
            for (j = start; j <= end; ++j) {
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];

                d = mscore > F[j] ? 0 : 1;
                M2[j] = mscore > F[j] ? mscore : F[j];

                d = M2[j] > e1 ? d : 2;
                M2[j] = M2[j] > e1 ? M2[j] : e1;


                if ( M2[j] > maxScore) { // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2 = j;
                }

                if (j>=endPosition2 && M2[j] > thisMax){
                    thisMax = M2[j];
                    thisMaxj = j;
                }

                int32_t h = M2[j] + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F[j] + score.getExtendGapPenalty1(weights[i]);
                d |= f >= h ? 1 << 3 : 0;
                f = f >= h ? f : h;
                F[j] = f;

                h = M2[j] + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                d |= e1 >= h ? 1 << 4 : 0;
                e1 = e1 >= h ? e1 : h;

                T[i][j] = d; // z[i,j]
            }
        }else{
            for ( j=start; j<=end; ++j ){
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];
                M2[j] = mscore > F[j] ? mscore : F[j];
                M2[j] = M2[j] > e1 ? M2[j] : e1;

                if( M2[j] >0 && M2[j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                    // >= will omit the first similar fragments
                    maxScore = M2[j];
                    endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2=j;
                }

                if (j>=endPosition2 && M2[j] > thisMax){
                    thisMax = M2[j];
                    thisMaxj = j;
                }

                int32_t h = M2[j] + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F[j] + score.getExtendGapPenalty1(weights[i]);
                f  = f >= h? f    : h;
                F[j] = f;

                h = M2[j] + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                e1  = e1 >= h? e1  : h;
            }
        }
        if( thisMaxj > endPosition2 && i > endPosition1 ){
            if( maxScore-thisMax > zdrop ) {
                break;
            }
        }
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
    }

    std::vector<uint32_t> cigar;

    if( returnCigar ){
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int ii = endPosition1;
        int jj = endPosition2;
        int tmp, state=0;

        while (ii>0 && jj>0) {
            tmp = T[ii][jj];
            if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
            else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
            if (state == 0) state = tmp & 7;
            if (state == 0){
                op=0;
                --ii;
                --jj;
            }else if (state == 1 || state == 3){
                op =2;
                --ii;
            }else {
                op =1;
                --jj;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        std::reverse(cigar.begin(),cigar.end());
    }

    delete[] M1;
    delete[] M2;
    delete[] F;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}
