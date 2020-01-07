//
// Created by Baoxing song on 2019-01-02.
//


#include "findSimilarFragmentsForPairedSequence.h"


// this function assume start1, end1, start2 and end2 are 1 based coordinate
std::vector<uint32_t> x_extend(int32_t & start1, int32_t & end1, int32_t & start2, int32_t & end2,
        int32_t & maxScore, seqan::CharString seq1, seqan::CharString seq2, std::string & alignment1,
        std::string & alignment2, const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                               const int & mismatchingPenalty, const double & length1d, const double & length2d,
                               const int & lDiag, const int & uDiag, const int & xDrop,
                               const double & kValue, const double & lambda, const double & pvalues, double & eValue, double & pvalue ){
    //std::cout << "before " << start1 << " " << end1 << " " << start2 << " " << end2 << " " << maxScore << std::endl;
    --start1;
    --start2;
    seqan::Score<int> sc(matchingScore, mismatchingPenalty, _extend_gap_penalty, _open_gap_penalty);
    seqan::Align<seqan::Infix<seqan::CharString const>::Type> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), infix(seq1, start1, end1));
    assignSource(row(align, 1), infix(seq2, start2, end2));
    int score = globalAlignment(align, sc);
    //std::cout << align;

    seqan::Tuple<unsigned, 4> positions = { {reinterpret_cast<uint32_t &>(start1), reinterpret_cast<uint32_t &>(start2), reinterpret_cast<uint32_t &>(end1), reinterpret_cast<uint32_t &>(end2)} };

//    We keep track of the best alignment score, denoted T , detected for a grid point
//lying on an antidiagonal before the current one. If we discover that S(i, j) , T ¡ X , where X ¶ 0 is
//user-speci ed, then we set S(i, j) to ¡ 1 so as to guarantee that S(i, j) will not play a role in subsequent
//evaluations of S.


    score = extendAlignment(align, score, seq1, seq2, positions, seqan::EXTEND_BOTH, lDiag, uDiag, xDrop, sc);
    // this is the reference document, I am using for those euqations
    // https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
    //std::cout << align;
    eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(score));
    pvalue = 1.0 - exp(-1.0*eValue);

    start1 = clippedBeginPosition(row(align, 0)) + 1;
    end1 = endPosition(row(align, 0));
    start2 = clippedBeginPosition(row(align, 1)) + 1;
    end2 = endPosition(row(align, 1));
    maxScore = score;

    std::vector<uint32_t> cigar;
    if( pvalue < pvalues ){
        for ( size_t i=0; i<length(row(align, 0)); ++i) {
            uint32_t op = 0;
            uint32_t length = 1;
            if( row(align, 0)[i] == '-' ){
                op = 1;
            }else if( row(align, 1)[i] == '-' ){
                op = 2;
            }
            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
            alignment1 += row(align, 0)[i];
            alignment2 += row(align, 1)[i];
        }
    }
    return cigar;
}



// do not use this function anymore, since it could not deal with copy number very well

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & mini_cns_size, const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                                                                           const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                                                                           const int & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int & lDiag, const int & uDiag, const int & xDrop){

    if( windowsSize > 2* mini_cns_size ){
        std::cerr << "the windows size is larger than 2*mini_cns_size, this might missing some similar fragments" << std::endl;
    }
    if( step_size*2 > mini_cns_score ){
        std::cerr << "step_size* is larger than mini_cns_score, the sequence alignment method would not work properly" << std::endl;
    }



    seqan::CharString charString_seq1 = seq1_string;
    seqan::CharString charString_seq2 = seq2_string;

    std::vector<PairedSimilarFragment> pairedSimilarFragments;

    int32_t startPosition2 = 0;
    int32_t this_windowsSize = windowsSize > length2 ? length2 : windowsSize;
    int32_t maxScore;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;


    while( ((startPosition2 + this_windowsSize) < length2) ) { // could not put == here, since the first one is startPosition2 + 0
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, endPosition1, endPosition2, m, true);
        while ((this_windowsSize - endPosition2) <= matrix_boundary_distance &&
               ((startPosition2 + this_windowsSize) < length2) && maxScore < mini_cns_score) { // once it reach mini_cns_score, it suggest this region could be aligned, and begin to extend it
            this_windowsSize += step_size;
            this_windowsSize = this_windowsSize < (length2-startPosition2) ? this_windowsSize : (length2-startPosition2);
            // std::cout << "line 212 startPosition2: " << startPosition2 << " this_windowsSize: " << this_windowsSize << std::endl;
            SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                          maxScore, endPosition1, endPosition2, m, true);
        }

        if ( mini_cns_score <= maxScore ) { //
            int32_t oldMaxScore = maxScore;
            //std::cout << "line 221 maxScore " << maxScore << std::endl << std::endl;
            SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + ( length2 - (startPosition2 + endPosition2 ) - 1 ),
                    1+endPosition1, endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1_rc, endPosition2_rc, m, true);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            // std::cout << "line 225 maxScore " << maxScore << std::endl << std::endl;
            if( (endPosition1_rc+1)>mini_cns_size && (endPosition2_rc+1)>mini_cns_size && endPosition1>=endPosition1_rc
                && endPosition2>=endPosition2_rc){ // why endPosition1>=endPosition1_rc and endPosition2>=endPosition2_rc
                int32_t start1 = endPosition1+1-endPosition1_rc+1;    // this 1 based position
                int32_t end1 = endPosition1+1;
                int32_t start2 = startPosition2+endPosition2-endPosition2_rc+1+1; // this 1 based position
                int32_t end2 = startPosition2+endPosition2+1;
                std::string alignment1;
                std::string alignment2;

                std::vector<uint32_t > cigar = x_extend(start1, end1, start2, end2, maxScore, charString_seq1, charString_seq2,
                        alignment1, alignment2, _open_gap_penalty, _extend_gap_penalty, matchingScore, mismatchingPenalty,
                        length1d, length2d, lDiag, uDiag, xDrop, kValue, lambda, pvalues, eValue, pvalue);

                if( pvalue < pvalues ){
                    PairedSimilarFragment pairedSimilarFragment(start1, end1, start2, end2, maxScore, cigar, alignment1, alignment2, pvalue, eValue);
                    // here the ordinate put into pairedSimilarFragment is 1 based
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                    // std::cout << "line 230 maxScore " << maxScore << std::endl << std::endl;
                }
                startPosition2 = end2-1; // here is a 2 base pair overlap
            }else{
                startPosition2 += endPosition2;
            }
            this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;
        }else{
            startPosition2 += windowsSize/2;
            this_windowsSize=windowsSize;
        }
    }
    // std::cout << "line 238" << std::endl;
    this_windowsSize = length2-startPosition2;
    if( this_windowsSize > mini_cns_size ){
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, endPosition1, endPosition2, m, true);
        // std::cout << "line 243" << std::endl;
        if ( mini_cns_score <= maxScore ) { //
            int32_t oldMaxScore = maxScore;
            SmithWaterman(seq1_rev_com+(length1-1-endPosition1), seq2_rev_com+(length2-(startPosition2+endPosition2)-1), 1+endPosition1,
                          endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore,
                          endPosition1_rc, endPosition2_rc, m, true);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            if( (endPosition1_rc+1)>mini_cns_size && (endPosition2_rc+1)>mini_cns_size ){
                int32_t start1 = endPosition1+1-endPosition1_rc + 1; // this 1 based position
                int32_t end1 = endPosition1+1;
                int32_t start2 = startPosition2+endPosition2-endPosition2_rc+1 + 1;// this 1 based position
                int32_t end2 = startPosition2+endPosition2+1;
                std::string alignment1;
                std::string alignment2;

                std::vector<uint32_t > cigar = x_extend(start1, end1, start2, end2, maxScore, charString_seq1, charString_seq2,
                                               alignment1, alignment2, _open_gap_penalty, _extend_gap_penalty,
                                               matchingScore, mismatchingPenalty, length1d, length2d, lDiag, uDiag,
                                               xDrop, kValue, lambda, pvalues, eValue, pvalue);
                if( pvalue < pvalues ) {
                    PairedSimilarFragment pairedSimilarFragment(start1, end1, start2, end2, maxScore, cigar, alignment1,
                                                                alignment2, pvalue, eValue);
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                }
            }
        }
    }
    return pairedSimilarFragments;
}



class StartAlignment{
    private:
        int32_t refPosition;
        int32_t queryPosition;
    public:
        StartAlignment(int32_t & _refPosition, int32_t & _queryPosition){
            refPosition = _refPosition;
            queryPosition = _queryPosition;
        }
        int32_t getRefPosition() const {
            return refPosition;
        }
        void setRefPosition(int32_t & _refPosition) {
            refPosition = _refPosition;
        }
        int32_t getQueryPosition() const {
            return queryPosition;
        }
        void setQueryPosition(int32_t & _queryPosition) {
            queryPosition = _queryPosition;
        }
        bool operator< (const StartAlignment & c) const {
            if ( this->refPosition < c.getRefPosition() ){
                return true;
            }else if( this->refPosition == c.getRefPosition() ){
                return this->queryPosition < c.getQueryPosition();
            }
            return false;
        }
};


class Seed{
    private:
        int32_t start1;
        int32_t end1;
        int32_t start2;
        int32_t end2;
    public:
        Seed(int32_t & _start1, int32_t & _end1, int32_t & _start2, int32_t & _end2){
            start1 = _start1;
            end1 = _end1;
            start2 = _start2;
            end2 = _end2;
        }
        int32_t getStart1() const {
            return start1;
        }
        void setStart1(int32_t & start1) {
            Seed::start1 = start1;
        }
        int32_t getEnd1() const {
            return end1;
        }
        void setEnd1(int32_t & end1) {
            Seed::end1 = end1;
        }
        int32_t getStart2() const {
            return start2;
        }
        void setStart2(int32_t & start2) {
            Seed::start2 = start2;
        }
        int32_t getEnd2() const {
            return end2;
        }
        void setEnd2(int32_t & end2) {
            Seed::end2 = end2;
        }
        bool operator< (const Seed & c) const {
            if( this->start1 == c.start1 && this->start2 == c.start2 ){
                return false;
            }
            if( this->end1 == c.end1 && this->end2 == c.end2 ){
                return false;
            }
            if ( this->start1 < c.start1 ){
                return true;
            }else if( this->start1 == c.start1 ){
                return this->start2 < c.start2;
            }
            return false;
        }
};


std::set<Seed> getSeeds( int8_t * seq1, int8_t * seq2, int8_t * seq1_rev_com, int8_t * seq2_rev_com,
                           const int32_t & length1, const int32_t & length2, const int & _open_gap_penalty1,
                           const int & _extend_gap_penalty1,
                           const int32_t & startPosition2, const int32_t & mini_cns_score,
                         const int32_t & windowsSize, const Scorei & m, const int32_t & step_size ){

    std::set<Seed> seeds;
    int32_t maxScore;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;
    int32_t this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;

    SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1, _extend_gap_penalty1,
                  maxScore, endPosition1, endPosition2, m, false);
//    while ((this_windowsSize - endPosition2) <= matrix_boundary_distance &&
//           ((startPosition2 + this_windowsSize) < length2) && maxScore < mini_cns_score) {
//        // once it reach mini_cns_score, it suggest this region could be aligned, and begin to extend it
//        this_windowsSize += step_size;
//        this_windowsSize = this_windowsSize < (length2 - startPosition2) ? this_windowsSize : (length2 - startPosition2);
//        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1,
//                      _extend_gap_penalty1, maxScore, endPosition1, endPosition2, m, false);
//    }

    if (mini_cns_score <= maxScore) {
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1,
                      _extend_gap_penalty1, maxScore, endPosition1, endPosition2, m, true);
        //int32_t oldMaxScore = maxScore;
        SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1),
                      seq2_rev_com + (length2 - (startPosition2 + endPosition2) - 1),
                      1 + endPosition1, endPosition2 + 1, _open_gap_penalty1, _extend_gap_penalty1, maxScore,
                      endPosition1_rc, endPosition2_rc, m, true);// this is the reverse alignment
        //assert(oldMaxScore == maxScore);
        //if(seed_max_size >= endPosition1_rc && seed_max_size>=endPosition2_rc){
            int32_t end1 = endPosition1; // 0based
            int32_t end2 = startPosition2+endPosition2; //0-based
            int32_t start1 = endPosition1-endPosition1_rc; //0-based coordinate
            int32_t start2 = startPosition2+endPosition2-endPosition2_rc; //0-based coordinate
            Seed seed0(start1, end1, start2, end2);
            seeds.insert(seed0);
            std::set<Seed> seeds1 = getSeeds( seq1, seq2, seq1_rev_com+(length1-start1), seq2_rev_com, start1,
                    length2, _open_gap_penalty1, _extend_gap_penalty1,
                    startPosition2, mini_cns_score, windowsSize, m, step_size );
            if ( !seeds1.empty() ){
                for( Seed seed1 : seeds1 ){
                    seeds.insert(seed1);
                }
            }
            std::set<Seed> seeds2 = getSeeds( seq1+end1+1, seq2, seq1_rev_com, seq2_rev_com, length1-end1-1,length2,
                    _open_gap_penalty1, _extend_gap_penalty1, startPosition2, mini_cns_score, windowsSize, m, step_size );
            if ( !seeds2.empty() ){
                for( Seed seed : seeds2 ){
                    int32_t start11 = seed.getStart1()+end1+1;
                    int32_t start21 = seed.getStart2();
                    int32_t end11 = seed.getEnd1()+end1+1;
                    int32_t end21 = seed.getEnd2();
                    Seed seed2(start11, end11, start21, end21);
                    seeds.insert(seed2);
                }
            }
        //}
    }
    return seeds;
}


void getSeeds( int8_t * seq1, int8_t * seq2, int8_t * seq1_rev_com, int8_t * seq2_rev_com, const int32_t & length1,
        const int32_t & length2, const int & _open_gap_penalty1, const int & _extend_gap_penalty1,
        const int & matrix_boundary_distance, const int32_t & startPosition2, const int32_t & mini_cns_score,
        const int32_t & windowsSize, int32_t & this_windowsSize_return, const Scorei & m,
        const int32_t & step_size, std::set<Seed> & seeds){

    int32_t maxScore;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;
    int32_t this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;

    SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1, _extend_gap_penalty1,
                  maxScore, endPosition1, endPosition2, m, false);
    while ((this_windowsSize - endPosition2) <= matrix_boundary_distance &&
           ((startPosition2 + this_windowsSize) < length2) && maxScore < mini_cns_score) { // once it reach mini_cns_score, it suggest this region could be aligned, and begin to extend it
        this_windowsSize += step_size;
        this_windowsSize = this_windowsSize < (length2 - startPosition2) ? this_windowsSize : (length2 - startPosition2);
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1,
                      _extend_gap_penalty1, maxScore, endPosition1, endPosition2, m, false);
    }

    if (mini_cns_score <= maxScore) {
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1,
                      _extend_gap_penalty1, maxScore, endPosition1, endPosition2, m, true);
        //int32_t oldMaxScore = maxScore;
        SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1),
                      seq2_rev_com + (length2 - (startPosition2 + endPosition2) - 1),
                      1 + endPosition1, endPosition2 + 1, _open_gap_penalty1, _extend_gap_penalty1, maxScore,
                      endPosition1_rc, endPosition2_rc, m, true);// this is the reverse alignment
        //assert(oldMaxScore == maxScore);
        //if ((endPosition1_rc ) < mini_cns_size && (endPosition2_rc ) < mini_cns_size) {
        int32_t end1 = endPosition1; // 0based
        int32_t end2 = startPosition2+endPosition2; //0-based
        int32_t start1 = endPosition1-endPosition1_rc; //0-based coordinate
        int32_t start2 = startPosition2+endPosition2-endPosition2_rc; //0-based coordinate

        Seed seed0(start1, end1, start2 , end2); //to
        if( seeds.find(seed0) == seeds.end()){
            seeds.insert(seed0);
            std::set<Seed> seeds1 = getSeeds( seq1, seq2, seq1_rev_com+(length1-start1), seq2_rev_com, start1,
                    length2, _open_gap_penalty1, _extend_gap_penalty1,
                    startPosition2, mini_cns_score, windowsSize, m, step_size );
            if ( !seeds1.empty() ){
                for( Seed seed1 : seeds1 ){
                    seeds.insert(seed1);
                }
            }
            std::set<Seed> seeds2 = getSeeds( seq1+end1+1, seq2, seq1_rev_com, seq2_rev_com, length1-end1-1,
                    length2, _open_gap_penalty1, _extend_gap_penalty1,
                    startPosition2, mini_cns_score, windowsSize, m, step_size );
            if ( !seeds2.empty() ){
                for( Seed seed : seeds2 ){
                    int32_t start11 = seed.getStart1()+end1+1;
                    int32_t start21 = seed.getStart2();
                    int32_t end11 = seed.getEnd1()+end1+1;
                    int32_t end21 = seed.getEnd2();
                    Seed seed2(start11, end11, start21, end21);
                    seeds.insert(seed2);
                }
            }
        }
        //}
    }
    this_windowsSize_return = this_windowsSize;
}


void x_extend_seed(uint32_t & start1, uint32_t & end1, uint32_t & start2, uint32_t & end2,
                               int32_t & maxScore, seqan::CharString & seq1, seqan::CharString & seq2,
                               seqan::Score<int> & sc, seqan::Align<seqan::Infix<seqan::CharString const>::Type> & align,
                               const double & length1d, const double & length2d,
                               const int & lDiag, const int & uDiag, const int & xDrop,
                               const double & kValue, const double & lambda, double & eValue, double & pvalue ){

    assignSource(row(align, 0), infix(seq1, start1, end1+1));
    assignSource(row(align, 1), infix(seq2, start2, end2+1));
    maxScore = globalAlignment(align, sc);
    seqan::Tuple<unsigned, 4> positions = { {start1, start2, end1+1, end2+1} };
    maxScore = extendAlignment(align, maxScore, seq1, seq2, positions, seqan::EXTEND_BOTH, lDiag, uDiag, xDrop, sc);
    // this is the reference document, I am using for those euqations
    // https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
    eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
    pvalue = 1.0 - exp(-1.0*eValue);

    start1 = clippedBeginPosition(row(align, 0)); //0 based
    end1 = endPosition(row(align, 0))-1; //0 based
    start2 = clippedBeginPosition(row(align, 1)); //0 based
    end2 = endPosition(row(align, 1))-1; //0 based
}


// not sure this better or above one better
void x_extend_seed(int32_t & start1, int32_t & end1, int32_t & start2, int32_t & end2,
                   int32_t & maxScore, int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com,
                   const double & length1d, const double & length2d,const int32_t & length1,
                   const int32_t & length2,const int & _open_gap_penalty1, const int & _extend_gap_penalty1,
                   const int32_t & xdrop, const int32_t & w,
                   const double & kValue, const double & lambda, double & eValue, double & pvalue, const Scorei & m ){

    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;
    int32_t s1 = start1;
    int32_t s2 = start2;

    SemiGlobal_xextend(seq1+start1, seq2+start2, length1-start1, length2-start2, _open_gap_penalty1, _extend_gap_penalty1,
    maxScore, endPosition1, endPosition2, m, xdrop, w);

    SemiGlobal_xextend(seq1_rev_com + (length1 - 1 - (start1+endPosition1)),
                  seq2_rev_com + (length2 - 1 - (start2 + endPosition2)),
                  start1+1 + endPosition1, start2+1+endPosition2, _open_gap_penalty1, _extend_gap_penalty1, maxScore,
                  endPosition1_rc, endPosition2_rc, m,  xdrop, w);// this is the reverse alignment

    end1 = s1 + endPosition1; // 0based
    end2 = s2 + endPosition2; //0-based
    start1 = s1+endPosition1-endPosition1_rc; //0-based coordinate
    start2 = s2+endPosition2-endPosition2_rc; //0-based coordinate
    eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
    pvalue = 1.0 - exp(-1.0*eValue);
}




// merge seeds into larger seeds
void link_seeds( std::set<Seed> & seeds ){
    std::vector<Seed> old_seeds;
    for( Seed s : seeds ){
        old_seeds.push_back(s);
    }
    std::sort(old_seeds.begin(), old_seeds.end(), [](Seed a, Seed b) {
        return a.getStart1() < b.getStart1();
    });

    bool everChange = true;
    int i, j;
    while(everChange){
        std::vector<Seed> new_seeds;
        everChange = false;
        std::set<int> usedList;
        for( i=0; i<old_seeds.size(); ++i ){
            if( usedList.find(i) == usedList.end() ){
                Seed s1 = old_seeds[i];
                bool ever_s1_used = false;
                for( j=i+1; j<old_seeds.size(); ++j ){
                    if( usedList.find(j) == usedList.end() ) {
                        Seed s2 = old_seeds[j];
                        if( s1.getStart2() < s2.getStart2() && s1.getEnd1()>s2.getEnd1() && s1.getEnd2()>s2.getEnd2()){ // s2 is a part of s1
                            usedList.insert(j);
                        }else if ( s1.getStart2() < s2.getStart2() && s1.getEnd1()<s2.getEnd1() && s1.getEnd2()<s2.getEnd2()  ) {
                            int a = s2.getStart1() - s1.getStart1();
                            if ( (a==s2.getStart2()-s1.getStart2()) && ( a==s2.getEnd1()-s1.getEnd1() ) && (a==s2.getEnd2()-s1.getEnd2()) ){ // merger s1 and s2
                                int32_t start1 = s1.getStart1();
                                int32_t start2 = s1.getStart2();
                                int32_t end1 = s2.getEnd1();
                                int32_t end2 = s2.getEnd2();
                                Seed s(start1, end1, start2, end2);
                                new_seeds.push_back(s);
                                usedList.insert(i);
                                usedList.insert(j);
                                ever_s1_used = true;
                                j = old_seeds.size(); //stop here
                            }
                        }else if ( s2.getStart1() > s1.getEnd1() || s2.getStart2() > s1.getEnd2() ) {
                            break;
                        }
                    }
                }
                if(  ever_s1_used ){
                    everChange = true;
                }else{
                    new_seeds.push_back(s1);
                }
            }
        }
        old_seeds = new_seeds;
    }
    seeds.clear();
    for( Seed s : old_seeds ){
        seeds.insert(s);
    }
}

void getAllExtendSeed( int8_t * seq1, int8_t * seq1_rev_com,
                       int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                       const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                       const int & _open_gap_penalty1, const int & _extend_gap_penalty1, const int & matchingScore,
                       const int & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                       std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                       const double & lambda, const double & kValue, const int32_t & zDrop, const int32_t & bandwidth,
                       const int & lDiag, const int & uDiag, const int & xDrop,  std::vector<Seed> & x_extend_seeds ){


    seqan::Align<seqan::Infix<seqan::CharString const>::Type> align;
    resize(rows(align), 2);

    int32_t startPosition2 = 0;
    int32_t this_windowsSize_return;
    int32_t maxScore;


    int32_t start1;
    int32_t end1;
    int32_t start2 ;
    int32_t end2;

    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;

    std::set<Seed> seeds;
    std::set<Seed> linked_seeds;

    seqan::CharString charString_seq1 = seq1_string;
    seqan::CharString charString_seq2 = seq2_string;
    seqan::Score<int> sc(matchingScore, mismatchingPenalty, _extend_gap_penalty1, _open_gap_penalty1);

    bool notEnd = true;
    while( notEnd ) { // could not put == here, since the first one is startPosition2 + 0
        getSeeds(seq1, seq2, seq1_rev_com, seq2_rev_com, length1, length2,
                 _open_gap_penalty1, _extend_gap_penalty1, matrix_boundary_distance,
                 startPosition2, mini_cns_score, windowsSize,
                 this_windowsSize_return, m, step_size, seeds);
        startPosition2 += (this_windowsSize_return - windowsSize  + step_size);
//        startPosition2 += step_size;
        if( (startPosition2 + windowsSize) > length2 ){
            getSeeds(seq1, seq2, seq1_rev_com, seq2_rev_com, length1, length2,
                     _open_gap_penalty1, _extend_gap_penalty1, matrix_boundary_distance,
                     startPosition2, mini_cns_score, windowsSize,
                     this_windowsSize_return, m, step_size, seeds);
            notEnd = false;
            break;
        }
    }
    link_seeds( seeds );  // this function is very fast, by using less seeds, the program is significantly faster

    for (Seed seed : seeds) {
        start1 = seed.getStart1();
        end1 = seed.getEnd1();
        start2 = seed.getStart2();
        end2 = seed.getEnd2();
//        std::cout << "line 550:" <<  start1 << " " << end1 << " " << start2 << " " << end2 << std::endl;
//        x_extend_seed(start1, end1, start2, end2, maxScore, charString_seq1, charString_seq2, sc, align,
//                            length1d, length2d, lDiag, uDiag, xDrop, kValue, lambda, eValue, pvalue );


        x_extend_seed(start1, end1, start2, end2, maxScore, seq1, seq1_rev_com, seq2, seq2_rev_com,
                      length1d, length2d, length1, length2, _open_gap_penalty1, _extend_gap_penalty1, xDrop, uDiag,
                      kValue, lambda, eValue, pvalue, m );

//        std::cout << "line 553:" <<  start1 << " " << end1 << " " << start2 << " " << end2 << std::endl;
        if (pvalue < pvalues ) {
            bool never_used = true;
            for ( int i=0; i<x_extend_seeds.size(); ++i ){
                if( x_extend_seeds[i].getStart1() == start1 && x_extend_seeds[i].getStart2() == start2 ){
                    never_used = false;
                    if( end1 > x_extend_seeds[i].getEnd1() && end2 > x_extend_seeds[i].getEnd2() ){
                        x_extend_seeds[i].setEnd1(end1);
                        x_extend_seeds[i].setEnd2(end2);
                    }
                }else if( x_extend_seeds[i].getEnd1() == end1 && x_extend_seeds[i].getEnd2() == end2 ){
                    never_used = false;
                }
            }
            if( never_used ){
                Seed x_seed(start1, end1, start2, end2);
                x_extend_seeds.push_back(x_seed);
            }
        }
    }
}




std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                   const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                   const int & _open_gap_penalty1, const int & _extend_gap_penalty1,const int & matchingScore,
                   const int & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                   std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                   const double & lambda, const double & kValue, const int & lDiag, const int & uDiag, const int & xDrop){

    std::vector<PairedSimilarFragment> pairedSimilarFragments;

    std::set<StartAlignment> startAlignments;
    uint8_t **T = new uint8_t *[length1+1];
    for (int32_t i = 0; i < (length1 + 1); ++i) {
        T[i] = new uint8_t[length2 + 1];
    }

    seqan::Align<seqan::Infix<seqan::CharString const>::Type> align;
    resize(rows(align), 2);

    int32_t startPosition2 = 0;
    int32_t this_windowsSize_return;
    int32_t maxScore;

    int32_t start1;
    int32_t end1;
    int32_t start2 ;
    int32_t end2;

    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;

    std::set<Seed> seeds;
    std::set<Seed> linked_seeds;

    seqan::CharString charString_seq1 = seq1_string;
    seqan::CharString charString_seq2 = seq2_string;
//    std::cout << "begin to generate seeds" << std::endl;
    bool notEnd = true;
    while( notEnd ) { // could not put == here, since the first one is startPosition2 + 0
        getSeeds(seq1, seq2, seq1_rev_com, seq2_rev_com, length1, length2,
                 _open_gap_penalty1, _extend_gap_penalty1, matrix_boundary_distance,
                 startPosition2, mini_cns_score, windowsSize,
                 this_windowsSize_return, m, step_size, seeds);
        startPosition2 += (this_windowsSize_return - windowsSize  + step_size);
        if( (startPosition2 + windowsSize) > length2 ){
            getSeeds(seq1, seq2, seq1_rev_com, seq2_rev_com, length1, length2,
                     _open_gap_penalty1, _extend_gap_penalty1, matrix_boundary_distance,
                     startPosition2, mini_cns_score, windowsSize,
                     this_windowsSize_return, m, step_size, seeds);
            notEnd = false;
            break;
        }
    }
    link_seeds( seeds );  // this function is very fast, by using less seeds, the program is significantly faster
//    std::cout << "generating seeds done" << std::endl;
    for (Seed seed : seeds) {
        start1 = seed.getStart1()+1;
        end1 = seed.getEnd1()+1;
        start2 = seed.getStart2()+1;
        end2 = seed.getEnd2()+1;

        std::string alignment1;
        std::string alignment2;

        std::vector<uint32_t > cigar = x_extend(start1, end1, start2, end2, maxScore, charString_seq1, charString_seq2,
                                                alignment1, alignment2, _open_gap_penalty1, _extend_gap_penalty1, matchingScore, mismatchingPenalty,
                                                length1d, length2d, lDiag, uDiag, xDrop, kValue, lambda, pvalues, eValue, pvalue);


        if (pvalue < pvalues ) {
            bool never_used = true;
            PairedSimilarFragment pairedSimilarFragment(start1, end1, start2, end2, maxScore, cigar, alignment1, alignment2, pvalue, eValue);
            for ( int i=0; i<pairedSimilarFragments.size(); ++i ){
                if( pairedSimilarFragments[i].getStart1() == start1 && pairedSimilarFragments[i].getStart2() == start2 ){
                    never_used = false;
                    if( end1 > pairedSimilarFragments[i].getEnd1() && end2 > pairedSimilarFragments[i].getEnd2() ){
                        pairedSimilarFragments[i] = pairedSimilarFragment;
                    }
                    break;
                }else if( pairedSimilarFragments[i].getEnd1() == end1 && pairedSimilarFragments[i].getEnd2() == end2 ){
                    never_used = false;
                    break;
                }
            }
            if( never_used ){
                pairedSimilarFragments.push_back(pairedSimilarFragment);
            }
        }
    }

    for (int32_t i = 0; i < (length1 + 1); ++i) {
        delete[] T[i];
    }
    delete[] T;

    return pairedSimilarFragments;
}


std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                const int & _open_gap_penalty1, const int & _extend_gap_penalty1,
                const int & _open_gap_penalty2, const int & _extend_gap_penalty2, const int & matchingScore,
                const int & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                const double & lambda, const double & kValue, const int32_t & zDrop, const int32_t & bandwidth,
                const int & _lDiag, const int & _uDiag, const int & xDrop){

    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    int32_t maxScore;

    int32_t start1;
    int32_t end1;
    int32_t start2;
    int32_t end2;
    std::set<StartAlignment> startAlignments;
    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;
    uint8_t **T = new uint8_t *[length1+1];
    for (int32_t i = 0; i < (length1 + 1); ++i) {
        T[i] = new uint8_t[length2 + 1];
    }

    std::vector<Seed> x_extend_seeds;

    getAllExtendSeed(seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2, windowsSize,
                     mini_cns_score, matrix_boundary_distance, _open_gap_penalty1, _extend_gap_penalty1,
                     matchingScore, mismatchingPenalty, m, step_size,
                     seq1_string, seq2_string, pvalues, lambda, kValue, zDrop, bandwidth, _lDiag, _uDiag,
                     xDrop, x_extend_seeds );

//    std::cout << "line 769 " << std::endl;
    for( Seed x_seed :  x_extend_seeds){

        endPosition1 = x_seed.getEnd1();
        endPosition2 = x_seed.getEnd2();

        SemiGlobal(seq1_rev_com + (length1 - 1 - endPosition1),
                   seq2_rev_com + (length2 - endPosition2 - 1),
                   1 + endPosition1, endPosition2 + 1, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2,
                   _extend_gap_penalty2, maxScore, endPosition1_rc, endPosition2_rc, m, false, zDrop, bandwidth, T);// this is the reverse alignment
//        std::cout << "line 779 " << maxScore << std::endl;
        start1 = endPosition1 - endPosition1_rc; //0 based
        start2 = endPosition2 - endPosition2_rc ; // 0 based

        StartAlignment startAlignment(start1, start2);
        if (startAlignments.find(startAlignment) == startAlignments.end()) {
            startAlignments.insert(startAlignment);
//            std::cout <<"line 786" << std::endl;
            std::vector<uint32_t> cigar = SemiGlobal(seq1 + start1, seq2 + start2, length1 - start1,
                                                     length2 - start2, _open_gap_penalty1, _extend_gap_penalty1,
                                                     _open_gap_penalty2, _extend_gap_penalty2,
                                                     maxScore, end1, end2, m, true, zDrop, bandwidth, T);
            /*std::cout << "line 791 maxScore:" << maxScore << " cigar size:" << cigar.size() << "\t";
            for( int j=0; j<cigar.size(); ++j ){
                uint32_t cigarLength = cigar[j]>>4;
                uint32_t cigarType = cigar[j]&0xf;
                std::cout << cigarLength << "MIDDSHI=XB"[cigarType];
            }*/
//            std::cout << std::endl << "line 797 end1: " << end1 << " end2: " << end2 << std::endl;
            end1 += start1; // 0based
            end2 += start2; // 0based

            int toRemove = -1;
            bool found = false;
            for ( int i=0; i< pairedSimilarFragments.size(); ++i) {
                if ( pairedSimilarFragments[i].getEnd1() == end1+1 && pairedSimilarFragments[i].getEnd2() == end2+1 ){
                    found = true;
                    if( start1+1 < pairedSimilarFragments[i].getStart1() && start2+1 < pairedSimilarFragments[i].getStart2() ){
                        toRemove = i;
                    }
                }
            }
            if( !found or toRemove>=0 ){

                eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                pvalue = 1.0 - exp(-1.0*eValue);

                PairedSimilarFragment pairedSimilarFragment(start1+1, end1 + 1, start2+1, end2 + 1, maxScore, cigar,
                        pvalue, eValue); // start1， start2, end1 and end2 are 1 based
                // here the ordinate put into pairedSimilarFragment is 1 based
                if( toRemove>=0  ){
                    pairedSimilarFragments[toRemove] = pairedSimilarFragment;
                }else{
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                }
            }
        }
    }

    for (int32_t i = 0; i < (length1 + 1); ++i) {
        delete[] T[i];
    }
    delete[] T;

    return pairedSimilarFragments;
}





std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_wighted ( int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                   const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                   const int & _open_gap_penalty1, const int & _extend_gap_penalty1, const int & _open_gap_penalty2,
                   const int & _extend_gap_penalty2, const int & matchingScore, const int & mismatchingPenalty,
                   const Scorei & m, const int32_t & step_size, std::string & seq1_string, std::string & seq2_string,
                   const double & pvalues, const double & lambda, const double & kValue, const int32_t & zDrop,
                   const int32_t & bandwidth, int & _lDiag, int & _uDiag, const int & xDrop, Score & score,
                   int16_t * weight, int16_t * weight_rev){

    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    int32_t maxScore;

    int32_t start1;
    int32_t start2;
    std::set<StartAlignment> startAlignments;
    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;
    uint8_t **T = new uint8_t *[length1+1]; // this matrix if for trace backing, so int8_t is enough
    for (int32_t i = 0; i < (length1 + 1); ++i) {
        T[i] = new uint8_t[length2 + 1];
    }

    std::vector<Seed> x_extend_seeds;

    getAllExtendSeed(seq1, seq1_rev_com, seq2, seq2_rev_com, reinterpret_cast<int32_t &>(length1), length2, windowsSize,
                     mini_cns_score, matrix_boundary_distance, _open_gap_penalty1, _extend_gap_penalty1,
                     matchingScore, mismatchingPenalty, m, step_size,
                     seq1_string, seq2_string, pvalues, lambda, kValue, zDrop, bandwidth, _lDiag, _uDiag,
                     xDrop, x_extend_seeds );
//    std::cout << "line 868:" << x_extend_seeds.size() << std::endl;
    for( Seed x_seed :  x_extend_seeds){

        endPosition1 = x_seed.getEnd1();
        endPosition2 = x_seed.getEnd2();
        SemiGlobal(seq1_rev_com + (length1 - 1 - endPosition1),
                                                 seq2_rev_com + (length2 - endPosition2 - 1),
                                                 1 + endPosition1, endPosition2 + 1,
                                                 weight_rev + (length1 - 1 - endPosition1), score, maxScore,
                                                 endPosition1_rc, endPosition2_rc, false, zDrop, bandwidth, T);
//        std::cout << "line 878: maxScore:" << maxScore << " endPosition1_rc:" << endPosition1_rc << " endPosition2_rc:" << endPosition2_rc << std::endl;
        start1 = endPosition1 - endPosition1_rc; //0 based
        start2 = endPosition2 - endPosition2_rc ; // 0 based
//        std::cout << "line 881 start1:" << start1 << " start2:" << start2 << std::endl;
        StartAlignment startAlignment(start1, start2);
        if (startAlignments.find(startAlignment) == startAlignments.end()) {

            startAlignments.insert(startAlignment);

            std::vector<uint32_t> cigar = SemiGlobal(seq1 + start1, seq2 + start2, length1 - start1, length2 - start2,
                weight + start1, score, maxScore, endPosition1, endPosition2, true, zDrop, bandwidth, T);
//            std::cout << "line 889: endPosition1:" << endPosition1 << " endPosition2:" << endPosition2  << std::endl;
            endPosition1 += start1; // 0based
            endPosition2 += start2; // 0based
//            std::cout << "line 892: endPosition1:" << endPosition1 << " endPosition2:" << endPosition2  << std::endl;
            int toRemove = -1;
            bool found = false;
            for ( int i=0; i< pairedSimilarFragments.size(); ++i) {
                if ( pairedSimilarFragments[i].getEnd1() == endPosition1+1 && pairedSimilarFragments[i].getEnd2() == endPosition2+1 ){
                    found = true;
                    if( start1+1 < pairedSimilarFragments[i].getStart1() && start2+1 < pairedSimilarFragments[i].getStart2() ){
                        toRemove = i;
                    }
                }
            }

            if( !found or toRemove>=0 ){
                eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                pvalue = 1.0 - exp(-1.0*eValue);

                PairedSimilarFragment pairedSimilarFragment(start1+1, endPosition1 + 1, start2+1, endPosition2 + 1, maxScore, cigar,
                                                            pvalue, eValue); // start1， start2, end1 and end2 are 1 based
                // here the ordinate put into pairedSimilarFragment is 1 based
                if( toRemove>=0  ){
                    pairedSimilarFragments[toRemove] = pairedSimilarFragment;
                }else{
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                }
            }
        }
    }

    for (int32_t i = 0; i < (length1 + 1); ++i) {
        delete[] T[i];
    }
    delete[] T;

    return pairedSimilarFragments;
}





std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_wighted_1gap ( int8_t * seq1, int8_t * seq1_rev_com,
                       int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                       const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                       const int & _open_gap_penalty1, const int & _extend_gap_penalty1, const int & _open_gap_penalty2,
                       const int & _extend_gap_penalty2, const int & matchingScore, const int & mismatchingPenalty,
                       const Scorei & m, const int32_t & step_size, std::string & seq1_string, std::string & seq2_string,
                       const double & pvalues, const double & lambda, const double & kValue, const int32_t & zDrop,
                       const int32_t & bandwidth, int & _lDiag, int & _uDiag, const int & xDrop, Score & score,
                       int16_t * weight, int16_t * weight_rev){

    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    int32_t maxScore;

    int32_t start1;
    int32_t start2;
    std::set<StartAlignment> startAlignments;
    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;
    uint8_t **T = new uint8_t *[length1+1]; // this matrix if for trace backing, so int8_t is enough
    for (int32_t i = 0; i < (length1 + 1); ++i) {
        T[i] = new uint8_t[length2 + 1];
    }

    std::vector<Seed> x_extend_seeds;

    getAllExtendSeed( seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2, windowsSize,
                      mini_cns_score, matrix_boundary_distance, _open_gap_penalty1, _extend_gap_penalty1,
                      matchingScore, mismatchingPenalty, m, step_size,
                      seq1_string, seq2_string, pvalues, lambda, kValue, zDrop, bandwidth, _lDiag, _uDiag,
                      xDrop, x_extend_seeds );

//    std::cout << "line 970:" << x_extend_seeds.size() << std::endl;
    for( Seed x_seed :  x_extend_seeds){

        endPosition1 = x_seed.getEnd1();
        endPosition2 = x_seed.getEnd2();
//        std::cout << "line 975: endPosition1:" << endPosition1 << " endPosition2:" << endPosition2  << std::endl;
        SemiGlobal_single_gap_penalty(seq1_rev_com + (length1 - 1 - endPosition1),
                   seq2_rev_com + (length2 - endPosition2 - 1),
                   1 + endPosition1, endPosition2 + 1,
                   weight_rev + (length1 - 1 - endPosition1), score, maxScore,
                   endPosition1_rc, endPosition2_rc, false, zDrop, bandwidth, T);
//        std::cout << "line 981: maxScore:" << maxScore << " endPosition1_rc:" << endPosition1_rc << " endPosition2_rc:" << endPosition2_rc << std::endl;

        start1 = endPosition1 - endPosition1_rc; //0 based
        start2 = endPosition2 - endPosition2_rc ; // 0 based
//        std::cout << "line 985 start1:" << start1 << " start2:" << start2 << std::endl;
        StartAlignment startAlignment(start1, start2);
        if (startAlignments.find(startAlignment) == startAlignments.end()) {

            startAlignments.insert(startAlignment);

            std::vector<uint32_t> cigar = SemiGlobal_single_gap_penalty(seq1 + start1, seq2 + start2, length1 - start1, length2 - start2,
                                                     weight + start1, score, maxScore, endPosition1, endPosition2, true, zDrop, bandwidth, T);
//            std::cout << "line 993 maxScore:" << maxScore << " endPosition1:" << endPosition1 << " endPosition2:" << endPosition2 << std::endl;

            endPosition1 += start1; // 0based
            endPosition2 += start2; // 0based

            int toRemove = -1;
            bool found = false;
            for ( int i=0; i< pairedSimilarFragments.size(); ++i) {
                if ( pairedSimilarFragments[i].getEnd1() == endPosition1+1 && pairedSimilarFragments[i].getEnd2() == endPosition2+1 ){
                    found = true;
                    if( start1+1 < pairedSimilarFragments[i].getStart1() && start2+1 < pairedSimilarFragments[i].getStart2() ){
                        toRemove = i;
                    }
                }
            }
//            std::cout << "line 808" << std::endl;
            if( !found or toRemove>=0 ){
                eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                pvalue = 1.0 - exp(-1.0*eValue);

                PairedSimilarFragment pairedSimilarFragment(start1+1, endPosition1 + 1, start2+1, endPosition2 + 1, maxScore, cigar,
                                                            pvalue, eValue); // start1， start2, end1 and end2 are 1 based
                // here the ordinate put into pairedSimilarFragment is 1 based
                if( toRemove>=0  ){
                    pairedSimilarFragments[toRemove] = pairedSimilarFragment;
                }else{
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                }
            }
//            std::cout << "line 822" << std::endl;
        }
    }

    for (int32_t i = 0; i < (length1 + 1); ++i) {
        delete[] T[i];
    }
    delete[] T;

    return pairedSimilarFragments;
}







/*

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2, uint32_t & windowsSize,
                                                                           const uint32_t & mini_cns_size, const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                                                                           const int & _open_gap_penalty1, const int & _extend_gap_penalty1,
                                                                           const int & _open_gap_penalty2, const int & _extend_gap_penalty2, const int & matchingScore,
                                                                           const int & mismatchingPenalty, const Scorei & m, const uint32_t & step_size,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const uint32_t & zDrop, const int32_t & bandwidth){

    if( windowsSize > 2* mini_cns_size ){
        std::cerr << "the windows size is larger than 2*mini_cns_size, this might missing some similar fragments" << std::endl;
    }
    if( step_size*2 > mini_cns_score ){
        std::cerr << "step_size* is larger than mini_cns_score, the sequence alignment method would not work properly" << std::endl;
    }

    std::vector<PairedSimilarFragment> pairedSimilarFragments;

    uint32_t startPosition2 = 0;
    uint32_t this_windowsSize = windowsSize > length2 ? length2 : windowsSize;
    int32_t maxScore;
    uint32_t endPosition1;
    uint32_t endPosition2;
    uint32_t endPosition1_rc;
    uint32_t endPosition2_rc;

    uint32_t start1;
    uint32_t end1;
    uint32_t start2 ;
    uint32_t end2;
    std::string alignment1;
    std::string alignment2;


    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;

    int8_t **T = new int8_t *[length1+1];
    for (int32_t i = 0; i < (length1 + 1); ++i) {
        T[i] = new int8_t[length2 + 1];
    }
    std::set<StartAlignment> startAlignments;

    while( ((startPosition2 + this_windowsSize) < length2) ) { // could not put == here, since the first one is startPosition2 + 0
        // here we only try the largest score of the seed, maybe not sensitive enough
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1, _extend_gap_penalty1,
                      maxScore, endPosition1, endPosition2, m, true);
        while ((this_windowsSize - endPosition2) <= matrix_boundary_distance &&
               ((startPosition2 + this_windowsSize) < length2) && maxScore < mini_cns_score) { // once it reach mini_cns_score, it suggest this region could be aligned, and begin to extend it
            this_windowsSize += step_size;
            this_windowsSize = this_windowsSize < (length2-startPosition2) ? this_windowsSize : (length2-startPosition2);
            SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1, _extend_gap_penalty1,
                          maxScore, endPosition1, endPosition2, m, true); //todo end position1 and endposition2 is not correct here
        } //todo if there are multiple copies in the seq1, maybe could not be reported here

        if ( mini_cns_score <= maxScore ) {

                SemiGlobal(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + ( length2 - ( endPosition2 ) - 1 ),
                              1+endPosition1, endPosition2+1, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2,
                              _extend_gap_penalty2, maxScore, endPosition1_rc, endPosition2_rc, m, false, zDrop, bandwidth, T);// this is the reverse alignment

                start1 = endPosition1+1-endPosition1_rc; //1 based
                start2 = startPosition2+endPosition2-endPosition2_rc+1; // since startPosition2 is 0 based so add 1



                StartAlignment startAlignment(start1, start2);
                if( startAlignments.find(startAlignment) != startAlignments.end() ){

                }else {
                    startAlignments.insert(startAlignment);
                }

                SemiGlobal(seq1 + start1-1, seq2 + start2-1, length1-start1+1, length2-start2+1, _open_gap_penalty1, _extend_gap_penalty1,
                                                            _open_gap_penalty2, _extend_gap_penalty2,maxScore, end1, end2, m, false, zDrop, bandwidth, T);

                eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                pvalue = 1.0 - exp(-1.0*eValue);

                if( pvalue < pvalues ){
                    std::vector<uint32_t> cigar = SemiGlobal(seq1 + start1-1, seq2 + start2-1, length1-start1+1, length2-start2+1, _open_gap_penalty1, _extend_gap_penalty1,
                                                             _open_gap_penalty2, _extend_gap_penalty2,
                                                             maxScore, end1, end2, m, true, zDrop, bandwidth, T);
                    end1 += start1-1;
                    end2 += start2-1;
                    PairedSimilarFragment pairedSimilarFragment(start1, end1+1, start2, end2+1, maxScore, cigar, alignment1, alignment2, pvalue, eValue); // start1， start2, end1 and end2 are 1 based
                    // here the ordinate put into pairedSimilarFragment is 1 based
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                    // std::cout << "line 230 maxScore " << maxScore << std::endl << std::endl;
                }else{
                    end1 += start1-1;
                    end2 += start2-1;
                }
            startPosition2 += endPosition2;
            this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;
        }else{
            startPosition2 += windowsSize/2;
            this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;
        }
    }
    for (int32_t i = 0; i < (length1 + 1); ++i) {
        delete[] T[i];
    }
    delete[] T;

    return pairedSimilarFragments;
}
*/
