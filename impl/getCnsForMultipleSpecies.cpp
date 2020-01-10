//
// Created by Baoxing song on 2019-01-09.
//

#include "getCnsForMultipleSpecies.h"

/** TODO: should be updated
 * The problem with this one is that the MSA output must have the reference sequence there
 * The alignment without reference would not be used
 */


class Combination{
    private:
        std::vector<int> combine;
    public:
        Combination(std::vector<int> & c){
            combine = c;
        }
        const std::vector<int> & getCombine() const{
            return combine;
        }
        bool operator< (const Combination & c) const {
            assert(c.getCombine().size() == combine.size());
            for( int32_t i=0; i<combine.size(); ++i ){
                if (combine[i] < c.getCombine()[i]) {
                    return true;
                }
            }
            return false;
        }
    friend std::ostream& operator<<(std::ostream& os, const Combination& combine);
};

std::ostream& operator<<(std::ostream& os, const Combination& combine) {
    for( int32_t i=0; i<combine.getCombine().size(); ++i ){
        os<< "::" << combine.getCombine()[i];
    }
    os<<std::endl;
    return os;
}

/*
 * take to PairedSimilarFragment and check if they overlap with each other in the terms of reference range
 * return a bool value
 * **/
bool pairedSimilarFragmentsOverlap( const PairedSimilarFragment & pairedSimilarFragment1, const PairedSimilarFragment & pairedSimilarFragment2){
    if( pairedSimilarFragment2.getStart1()<=pairedSimilarFragment1.getStart1() && pairedSimilarFragment1.getStart1()<pairedSimilarFragment2.getEnd1() ){
        return true;
    }else if( pairedSimilarFragment1.getStart1()<=pairedSimilarFragment2.getStart1() && pairedSimilarFragment2.getStart1()<pairedSimilarFragment1.getEnd1() ){
        return true;
    }
    return false;
}


bool pairedSimilarFragmentsOverlap( const PairedSimilarFragment2 & pairedSimilarFragment1, const PairedSimilarFragment2 & pairedSimilarFragment2){
    if( pairedSimilarFragment2.getStart1()<=pairedSimilarFragment1.getStart1() && pairedSimilarFragment1.getStart1()<pairedSimilarFragment2.getEnd1() ){
        return true;
    }else if( pairedSimilarFragment1.getStart1()<=pairedSimilarFragment2.getStart1() && pairedSimilarFragment2.getStart1()<pairedSimilarFragment1.getEnd1() ){
        return true;
    }
    return false;
}

//the return is a matrix which lists all the combinations, like [[2,3,5], [4,5,7]]   2,3,5 is the index from allPairedSimilarFragments[1], allPairedSimilarFragments[2] and allPairedSimilarFragments[3]
// allPairedSimilarFragments[0] is the reference one, so do not use it here
std::vector< std::vector<int>> generateAllCombinations(std::vector<std::vector<bool>> & overlapped){
    std::vector< std::vector<int>> indexs(overlapped.size());
    for ( int32_t i=0; i < overlapped.size(); ++i ){
        for ( int32_t j=0; j<overlapped[i].size(); ++j ) {
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

void getCnsForMultipleSpecies ( int8_t ** seqs, int8_t ** seq_revs, std::vector<int32_t> & lengths,
                                int32_t & windowsSize, int32_t & mini_cns_seed_size, int32_t & mini_cns_score,
                                const int32_t & matrix_boundary_distance, const int32_t & _open_gap_penalty,
                                const int32_t & _extend_gap_penalty, const int32_t & matchingScore,
                                const int32_t & mismatchingPenalty, const bool & onlySyntenic, const std::string & output,
                                std::map<std::string, std::string>& sequences, std::vector<std::string> & seqNames,
                                const int32_t & step_size, std::vector<std::string> & seqs_string,
                                const int32_t & minimumNumberOfSpecies, const Scorei & m, const int32_t & mini_cns_size,
                                const double & outputWithMinimumLengthPercentage, const double & lambda,
                                const double & kValue, const int32_t & w, const int32_t & xDrop){

    if( lengths.size()<3 ){
        std::cerr << "you should have at least 3 sequences in your input fasta file" << std::endl;
        return;
    }

    //pair-wise sequence alignment begin
    std::vector<std::vector<PairedSimilarFragment>> allPairedSimilarFragments(lengths.size()-1);
    // the first sequence is used as reference and do not perform pairwise sequence alignment without the reference
    int8_t * seq1 = seqs[0];
    int8_t * seq1_rev_com = seq_revs[0];
    int32_t length1 = lengths[0];
    for ( int32_t ii=1; ii< lengths.size(); ++ii ){

        int8_t * seq2 = seqs[ii];
        int8_t * seq2_rev_com = seq_revs[ii];
        int32_t length2 = lengths[ii];

        std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence
                                                (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                windowsSize, mini_cns_score,
                                                matrix_boundary_distance, _open_gap_penalty, _extend_gap_penalty,
                                                matchingScore, mismatchingPenalty, m, step_size, seqs_string[0], seqs_string[ii], 1.0,
                                                lambda, kValue, w, xDrop);

        if ( onlySyntenic ){
            std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
            pairedSimilarFragments0 = pairedSimilarFragments;
        }
        for( int32_t iii=0; iii<pairedSimilarFragments0.size(); ++iii ){
            std::string alignment1 = "";
            std::string alignment2 = "";
            int32_t s1 = pairedSimilarFragments0[iii].getStart1()-1;
            int32_t s2 = pairedSimilarFragments0[iii].getStart2()-1;

            //transform those cigar to string begin
            for( int32_t j=0; j<pairedSimilarFragments0[iii].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[iii].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[iii].getCigar()[j]&0xf;
                if( cigarType == 0 ){
                    alignment1 += sequences[seqNames[0]].substr(s1, cigarLength);
                    alignment2 += sequences[seqNames[ii]].substr(s2, cigarLength);;
                    s1 += cigarLength;
                    s2 += cigarLength;
                }else if( cigarType == 1 ){
                    alignment1 += std::string(cigarLength, '-');
                    alignment2 += sequences[seqNames[ii]].substr(s2, cigarLength);;
                    s2 += cigarLength;
                }else if( cigarType == 2 ){
                    alignment1 += sequences[seqNames[0]].substr(s1, cigarLength);
                    alignment2 += std::string(cigarLength, '-');
                    s1 += cigarLength;
                }
            }
            //transform those cigar to string end
            pairedSimilarFragments0[iii].setAlignment1(alignment1);
            pairedSimilarFragments0[iii].setAlignment2(alignment2);
        }
        allPairedSimilarFragments[ii-1] = pairedSimilarFragments0;
    }
    //pair-wise sequence alignment end

    //sort according to the reference coordinate begin //remember those have been sorted, should be used for speeding up
    for ( int32_t i=0; i< allPairedSimilarFragments.size(); ++i ){
        std::sort(allPairedSimilarFragments[i].begin(), allPairedSimilarFragments[i].end(), [](PairedSimilarFragment a, PairedSimilarFragment b) {
            return a.getStart1() < b.getStart1();
        });
    }
    //sort according to the reference coordinate end

    std::set<Combination> allCombinationsSetUsed;
    std::ofstream ofile;
    ofile.open(output);

    std::vector<std::vector<bool>> overlapped(allPairedSimilarFragments.size());
    for (int32_t i = 0; i < allPairedSimilarFragments.size(); ++i) {
        overlapped[i].resize(allPairedSimilarFragments[i].size(), false);
    }
    for ( int32_t ai=0; ai< allPairedSimilarFragments.size(); ++ai ) {
        /**
         * If there are multiple records in the 3,4,...n species that overlapped with the record in the 2(second) one.
         * Firstly only care about the first overlapped one, then take the second third .... overlapping as reference to run everything again
         * **/

        for (int32_t i=0; i <allPairedSimilarFragments[ai].size(); ++i) {
            // the alignment must present in the reference, maybe not good
            // BUT ED suggest it make no sense to detect CNS no present in the reference sequence

            // set all the if overlap values begin
            for (int32_t j = 0; j < allPairedSimilarFragments.size(); ++j) {
                if( j==ai ){
                    overlapped[j].resize(allPairedSimilarFragments[j].size(), false);
                    overlapped[j][i] = true;
                }else{
                    for (int32_t k = 0; k < allPairedSimilarFragments[j].size(); ++k) {
                        if (pairedSimilarFragmentsOverlap(allPairedSimilarFragments[ai][i], allPairedSimilarFragments[j][k])) {
                            overlapped[j][k] = true;
                        } else if (allPairedSimilarFragments[j][k].getStart1() > allPairedSimilarFragments[ai][i].getEnd1()) {
                            // take the advantage of the vector has been sorted to speed up
                            for (; k < allPairedSimilarFragments[j].size(); ++k) {
                                overlapped[j][k] = false;
                            }
                        } else {
                            overlapped[j][k] = false;
                        }
                    }
                }
            }
            // set all the if overlap values end

            std::vector<std::vector<int>> allCombinations = generateAllCombinations(overlapped);  // matrix size: ### X (n-2)


            //overlapped has the same size with allPairedSimilarFragments
            std::vector<int> combination(allCombinations.size(), -1);
            while (getNext(allCombinations, combination)) {
                // if this combination has been checked, do not repeat it begin
                Combination c(combination);
                if (allCombinationsSetUsed.find(c) != allCombinationsSetUsed.end()) {
//                    std::cout << "line 192" << std::endl;
                    continue;
                }
                allCombinationsSetUsed.insert(c);
                // if this combination has been checked, do not repeat it end

                // find a good refStart and refEnd combination begin
                std::set<int> starts;
                std::set<int> ends;
                for (int32_t k = 0; k < combination.size(); ++k) {
                    if (allCombinations[k][combination[k]] > -1) { // this is a real overlap
                        starts.insert(allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1());
                        ends.insert(allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1());
                    }
                }

                int32_t refStart = 0;
                int32_t refEnd = 0;
                int32_t length = 0;
                for (int32_t start : starts) {
                    for (int32_t end : ends) {
                        std::set<std::string> species;
                        int32_t thisNumberOfSequences = 0;
                        int32_t newLength = end - start + 1;
                        if (length < newLength && mini_cns_size <= newLength) {
                            for (int32_t k = 0; k < combination.size(); ++k) {
                                if (allCombinations[k][combination[k]] > -1) {
                                    if (allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1() <=
                                        start && allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1() >= end) {
                                        ++thisNumberOfSequences;
                                        std::string spe = seqNames[k + 1]; // seqNames is 1 larger than allPairedSimilarFragments
                                        spe=spe.substr(0, spe.find("@", 0)-1);
                                        species.insert(spe);
                                    }
                                }
                            }
                            if (species.size() >=minimumNumberOfSpecies && thisNumberOfSequences >= minimumNumberOfSpecies /*&& thisNumberOfSpecies > numberOfSpecies*/ ) {
                                refStart = start;
                                refEnd = end;
                                length = newLength;
                            }
                        }
                    }
                }
                // find a good refStart and refEnd combination end

                if ((refEnd > refStart) && ((refEnd - refStart + 1) >= mini_cns_size)) {
                    int32_t newRefStart = refEnd; // since the one used for refEnd or refStart maybe not incldued in the final output, so recheck the refStart and refEnd here
                    int32_t newRedEnd = refStart;
                    std::set<std::string> species;
                    std::set<int> badIndex;

                    for (int32_t k = 0; k < combination.size(); ++k) {
                        if (allCombinations[k][combination[k]] > -1) { // this is a real overlap
                            int32_t thisStart = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1() > refStart ? allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1()
                                             : refStart;
                            int32_t thisEnd = allPairedSimilarFragments[k ][allCombinations[k][combination[k]]].getEnd1() <
                                      refEnd ? allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1()
                                             : refEnd;
                            double thisLength = thisEnd - thisStart + 1;
                            //                    std::cout << " line 216 thisLength: " << thisLength << " length " << length << std::endl;
                            if (thisLength < outputWithMinimumLengthPercentage * double(length)) {
                                badIndex.insert(k);
                            } else {
                                if (thisStart < newRefStart) {
                                    newRefStart = thisStart;
                                }
                                if (thisEnd > newRedEnd) {
                                    newRedEnd = thisEnd;
                                }
                                std::string spe = seqNames[k + 1];
                                spe = spe.substr(18, 5);
                                species.insert(spe);
                            }
                        }
                    }

                    if (newRefStart > refStart) {
                        refStart = newRefStart;
                    }
                    if (newRedEnd < refEnd) {
                        refEnd = newRedEnd;
                    }
                    // std::cout << " line 226  " << i << " species.size " << species.size() << " refStart " << refStart << " refEnd " << refEnd << std::endl;
                    //allPairedSimilarFragments[0] never output
                    if ( species.size() >= minimumNumberOfSpecies && (refEnd > refStart) && ((refEnd - refStart + 1) >= mini_cns_size)) {
                        std::vector<std::string> alignmentNames;
                        std::vector<std::string> alignmentSeqs;

                        //                std::cout << " line 307" << std::endl;
                        for (int32_t k = 0; k < combination.size(); ++k) {
                            if (allCombinations[k][combination[k]] > -1 && badIndex.find(k) == badIndex.end()) {
                                int32_t start2 = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart2();
                                int32_t end2 = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd2();

                                //int32_t endReference = allPairedSimilarFragments[ai][i].getEnd1();
                                int32_t endReference = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1();
                                int32_t positionReference = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1();
                                --positionReference;

                                std::string alignSeq;
                                std::string refSeq;
                                int32_t numberOfAltChar = 0;
                                bool alignRegionStart = false;
                                bool alignRegionEnd = true;
                                for (int32_t w = 0; w < allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment1().length(); ++w) {
                                    if (allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment1()[w] != '-') {
                                        ++positionReference;
                                    }
                                    if (!alignRegionStart) {
                                        if (positionReference > refStart && (!alignRegionStart)) {
                                            int32_t h = positionReference - refStart;
                                            while (h > 0) {
                                                alignSeq += ".";
                                                refSeq += seqs_string[0][positionReference - h - 1];
                                                --h;
                                            }
                                        }
                                        if (positionReference >= refStart) {
                                            alignRegionStart = true;
                                            start2 += numberOfAltChar;
                                            end2 = start2;
                                        }
                                    }
                                    if (alignRegionStart && alignRegionEnd) {
                                        alignSeq += allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment2()[w];
                                        refSeq += allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment1()[w];
                                    }
                                    if (allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment2()[w] != '-') {
                                        if (alignRegionStart && alignRegionEnd) {
                                            ++end2;
                                        }
                                        ++numberOfAltChar;
                                    }
                                    if (alignRegionEnd) {
                                        if (positionReference == endReference && positionReference < refEnd) {
                                            int32_t h = refEnd - positionReference;
                                            for (int32_t z = 0; z < h; ++z) {
                                                alignSeq += ".";
                                                refSeq += seqs_string[0][positionReference + z];
                                            }
                                            alignRegionEnd = false;
                                        }
                                        if (positionReference >= refEnd) {
                                            alignRegionEnd = false;
                                        }
                                    }
                                }
                                --end2;

                                if (alignmentSeqs.empty()) {
                                    alignmentNames.push_back(seqNames[0] + ":" + std::to_string(refStart) + "-" + std::to_string(refEnd)); //the reference one
                                    alignmentSeqs.push_back(refSeq);
                                } else { //update the alignment sequence by comparing the reference sequence
                                    int32_t lengthq = refSeq.length();
                                    for (int32_t o = 0; o < lengthq; ++o) {
                                        if (alignmentSeqs[0].length() > o) { // this is always the reference one
                                            if (refSeq[o] == '-' && alignmentSeqs[0][o] != '-') {
                                                for (int32_t p = 0; p < alignmentSeqs.size(); ++p) {
                                                    alignmentSeqs[p].insert(o, "-");
                                                }
                                            } else if (refSeq[o] != '-' && alignmentSeqs[0][o] == '-') {
                                                refSeq.insert(o, "-");
                                                alignSeq.insert(o, "-");
                                                ++lengthq;
                                            }
                                        }
                                    }
                                    while (refSeq.length() < alignmentSeqs[0].length()) {
                                        refSeq.insert(refSeq.length(), "-");
                                        alignSeq.insert(alignSeq.length(), "-");
                                    }
                                    while (refSeq.length() > alignmentSeqs[0].length()) {
                                        for (int32_t p = 0; p < alignmentSeqs.size(); ++p) {
                                            alignmentSeqs[p].insert(alignmentSeqs[p].length(), "-");
                                        }
                                    }
                                }
                                alignmentNames.push_back(seqNames[k + 1] + ":" + std::to_string(start2) + "-" + std::to_string(end2));
                                alignmentSeqs.push_back(alignSeq);
                            }
                        }
//                        ofile << c ;
                        for (int32_t o = 0; o < alignmentNames.size(); ++o) {
                            ofile << ">" << alignmentNames[o] << std::endl << alignmentSeqs[o] << std::endl;
                        }
                        ofile << std::endl << std::endl;
                    }
                }
            }
        }
    }
    ofile.close();
}



//todo set a parameter for each species, the maximum number of align records should be used for MSA
// should set up a score or something like that
//take sam files as input
void getCnsForMultipleSpecies ( const bool & onlySyntenic, const std::string & output,
                                //std::map<std::string, std::string> & sequences, /*species, fastaFile*/
                                std::map<std::string, std::string> & samFiles,/*species, sameFile*/
                                std::string & referenceGenomeFile,
                                const int32_t & minimumNumberOfSpecies, const int32_t & mini_cns_size,
                                const double & outputWithMinimumLengthPercentage){

    std::map<std::string, std::string> referenceGenome;
    readFastaFile( referenceGenomeFile, referenceGenome);
    std::cout << "reference genome reading done" << std::endl;

    std::map <std::string /*reference chr*/,  std::map<std::string /*species*/, std::vector<PairedSimilarFragment2>>> pairedSimilarFragments;

    for ( std::map<std::string, std::string>::iterator ii = referenceGenome.begin(); ii!=referenceGenome.end(); ++ii ){
        pairedSimilarFragments[ii->first] = std::map<std::string, std::vector<PairedSimilarFragment2>>();
        for ( std::map<std::string, std::string>::iterator i = samFiles.begin(); i!=samFiles.end(); ++i ){
            pairedSimilarFragments[ii->first][i->first] = std::vector<PairedSimilarFragment2>();
        }
    }

    for ( std::map<std::string, std::string>::iterator i = samFiles.begin(); i!=samFiles.end(); ++i ){
        std::string species = i->first;
        std::string samFile = i->second;

        std::ifstream infile(samFile);
        if( ! infile.good()){
            std::cerr << "error in opening sam file " << samFile << std::endl;
            exit (1);
        }

//        std::regex samRegex("^(.*?)\\s+(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\w+)\\s+(\\S+)\\s+(\\w+)\\s+(\\w+)\\s+(\\w+)\\s+");
//        std::smatch samMatch;
        std::regex cigarRegex ("(\\d+)([MIDNSH])");
        std::smatch cigarMatch;

        std::regex cigarRegex1 ("^(\\d+)([SH])");
        std::smatch cigarMatch1;

        int32_t refStart;
        int32_t refEnd;
        int32_t queryStart;
        int32_t queryEnd;
        int32_t s1;
        int32_t s2;
        int32_t cigarLength;
        std::string line;
        std::string cigarString;
        std::string cigarChar;
        std::string chromosomeName;

        std::cout << "begin to read " << samFile << std::endl;
        while (std::getline(infile, line)){
            if(line.size()<9 || line[0]=='@'){
            }else {
                std::vector<std::string> elems;
                char delim= '\t';
                split(line, delim, elems);

//                std::cout << "begin line 487: " << line << std::endl;
//                regex_search(line, samMatch, samRegex);
//                if (!samMatch.empty()) {
                if( elems.size() >10 ){
//                    cigarString = samMatch[6];
                    cigarString = elems[5];
//                    std::cout << "cigarString: " << cigarString << std::endl;
                    regex_search(cigarString, cigarMatch1, cigarRegex1);
                    if ( !cigarMatch1.empty() ) {
                        //std::string querySeq = samMatch[10];
                        std::string querySeq = elems[9];
//                        std::cout << "querySeq:" << querySeq << std::endl;

//                        std::string referenceChromosomeName = samMatch[3];
//                        std::string queryChromosomeName = samMatch[1];
//                        refStart = stoi(samMatch[4]);
                        std::string referenceChromosomeName = elems[2];
                        std::string queryChromosomeName = elems[0];
                        refStart = stoi(elems[3]);

                        queryStart = stoi(cigarMatch1[1]);
                        refEnd = refStart;
                        queryEnd = queryStart;
                        s1 = refStart-1;
                        s2 = 0;

                        std::string alignment1 = "";
                        std::string alignment2 = "";
                        while (std::regex_search(cigarString, cigarMatch, cigarRegex)) {
                            cigarChar = cigarMatch[2];
                            cigarLength = stoi(cigarMatch[1]);
                            if (cigarChar[0] == 'M'){
                                alignment1 += referenceGenome[referenceChromosomeName].substr(s1, cigarLength);
                                alignment2 += querySeq.substr(s2, cigarLength);

                                refEnd += cigarLength;
                                queryEnd += cigarLength;
                                s1 += cigarLength;
                                s2 += cigarLength;
                            }else if( cigarChar[0] == 'D' || cigarChar[0] == 'N') {
                                alignment1 += referenceGenome[referenceChromosomeName].substr(s1, cigarLength);
                                alignment2 += std::string(cigarLength, '-');
                                refEnd += cigarLength;
                                s1 += cigarLength;
                            }else if (cigarChar[0] == 'I') {
                                alignment1 += std::string(cigarLength, '-');
                                alignment2 += querySeq.substr(s2, cigarLength);
                                queryEnd += cigarLength;
                                s2 += cigarLength;
                            }
                            cigarString = cigarMatch.suffix().str();
                        }
                        int32_t cigarFlag = std::stoi(elems[1]);

//                        refStart = refStart - 1;
                        refEnd = refEnd - 1;
  //                      queryStart = queryStart - 1;
                        queryEnd = queryEnd - 1;

                        PairedSimilarFragment2 pairedSimilarFragment(species, queryChromosomeName, refStart, refEnd,
                                queryStart, queryEnd, alignment1, alignment2,1);

                        if ( 0 != cigarFlag % 32 ){ //negative strand
                            pairedSimilarFragment.setStrand(0);
                        }
                        pairedSimilarFragments[referenceChromosomeName][species].push_back(pairedSimilarFragment);

                    }
                }
            }
        }
        std::cout << samFile << " reading done" << std::endl;
    }
    // sam files reading done, get a data structure pairedSimilarFragments


    std::map <std::string /*reference chr*/,  std::vector<std::vector<PairedSimilarFragment2>>> pairedSimilarFragments1;
    for ( std::map <std::string,  std::map<std::string, std::vector<PairedSimilarFragment2>>>::iterator ii = pairedSimilarFragments.begin(); ii!=pairedSimilarFragments.end(); ++ii ){
        pairedSimilarFragments1[ii->first] = std::vector<std::vector<PairedSimilarFragment2>>();
        for ( std::map<std::string, std::vector<PairedSimilarFragment2>>::iterator i = ii->second.begin(); i!=ii->second.end(); ++i ){
            pairedSimilarFragments1[ii->first].push_back(i->second);
        }
    }
    // re-organized the data structure into pairedSimilarFragments1

    std::ofstream ofile;
    ofile.open(output);
    for( std::map <std::string, std::vector<std::vector<PairedSimilarFragment2>>>::iterator it=pairedSimilarFragments1.begin(); it!=pairedSimilarFragments1.end(); ++it ) {
        std::string referenceChr = it->first;
        std::cout << referenceChr << std::endl;
        std::vector<std::vector<PairedSimilarFragment2>> allPairedSimilarFragments = it->second; // each vector is the CNS in each species

        //sort according to the reference coordinate begin //remember those have been sorted, should be used for speeding up
        for ( int32_t i=0; i< allPairedSimilarFragments.size(); ++i ){
            std::sort(allPairedSimilarFragments[i].begin(), allPairedSimilarFragments[i].end(),
                      [](PairedSimilarFragment2 a, PairedSimilarFragment2 b) {
                          return a.getStart1() < b.getStart1();
                      });
        }
        //sort according to the reference coordinate end

        std::set<Combination> allCombinationsSetUsed;

        std::vector<std::vector<bool>> overlapped(allPairedSimilarFragments.size());
        for (int32_t i = 0; i < allPairedSimilarFragments.size(); ++i) {
            overlapped[i].resize(allPairedSimilarFragments[i].size(), false);
        }
        for ( int32_t ai=0; ai< allPairedSimilarFragments.size(); ++ai ) {
            /**
             * If there are multiple records in the 3,4,...n species that overlapped with the record in the 2(second) one.
             * Firstly only care about the first overlapped one, then take the second third .... overlapping as reference to run everything again
             * **/

            for (int32_t i = 0; i < allPairedSimilarFragments[ai].size(); ++i) {
                // the alignment must present in the reference, maybe not good
                // BUT ED suggest it make no sense to detect CNS no present in the reference sequence

                // here use allPairedSimilarFragments[ai][i] as reference, and check which CNS in different species overlapped with this one


                // set all the if overlap values begin
                for (int32_t j = 0; j < allPairedSimilarFragments.size(); ++j) {
                    if (j == ai) {
                        overlapped[j].resize(allPairedSimilarFragments[j].size(), false); // only use 1 in this accession
                        overlapped[j][i] = true;
                    } else {
                        for (int32_t k = 0; k < allPairedSimilarFragments[j].size(); ++k) {
                            if (pairedSimilarFragmentsOverlap(allPairedSimilarFragments[ai][i],
                                                              allPairedSimilarFragments[j][k])) {
                                overlapped[j][k] = true;
                            } else if (allPairedSimilarFragments[j][k].getStart1() >
                                       allPairedSimilarFragments[ai][i].getEnd1()) {
                                // take the advantage of the vector has been sorted to speed up
                                for (; k < allPairedSimilarFragments[j].size(); ++k) {
                                    overlapped[j][k] = false;
                                }
                            } else {
                                overlapped[j][k] = false;
                            }
                        }
                    }
                }
                // set all the if overlap values end


                // only keep one CNS for each species, else it would take too much RAM and CPU time
                // if you have 10 species, and for 8 species, there are 1000 CNS overlapped with the one taking as reference
                // the size of the index combination matrix would be 1000^8, take huge RAM. And the computing for all of those combinations would take a lot of time
                // for pratical purpose, for each species we only select the one which give longest overlap with the reference CNS
                for (int32_t j = 0; j < allPairedSimilarFragments.size(); ++j) {
                    if (j != ai) {
                        int32_t longest_length = 0;
                        int32_t longest_index = -1;
                        for (int32_t k = 0; k < allPairedSimilarFragments[j].size(); ++k) {
                            if( overlapped[j][k] ){
                                int32_t thisRefStart = allPairedSimilarFragments[ai][i].getStart1() > allPairedSimilarFragments[j][k].getStart1() ? allPairedSimilarFragments[ai][i].getStart1() : allPairedSimilarFragments[j][k].getStart1();
                                int32_t thisRefEnd = allPairedSimilarFragments[ai][i].getEnd1() < allPairedSimilarFragments[j][k].getEnd1() ? allPairedSimilarFragments[ai][i].getEnd1() : allPairedSimilarFragments[j][k].getEnd1();
                                int32_t thisRefLength = thisRefEnd - thisRefStart + 1;
                                if( thisRefLength > longest_length ){
                                    longest_length = thisRefLength;
                                    longest_index = k;
                                }
                            }
                        }
                        overlapped[j].resize(allPairedSimilarFragments[j].size(), false);
                        overlapped[j][longest_index] = true;
                    }
                }

                std::vector<std::vector<int32_t>> allCombinations = generateAllCombinations(overlapped);  // matrix size: ### X (n-2)
                // the size should be one now

                //overlapped has the same size with allPairedSimilarFragments
                std::vector<int32_t> combination(allCombinations.size(), -1);
                while (getNext(allCombinations, combination)) {
                    // if this combination has been checked, do not repeat it begin
                    Combination c(combination);
                    if (allCombinationsSetUsed.find(c) != allCombinationsSetUsed.end()) {
                        //                    std::cout << "line 192" << std::endl;
                        continue;
                    }
                    allCombinationsSetUsed.insert(c);
                    // if this combination has been checked, do not repeat it end

                    // find a good refStart and refEnd combination begin
                    std::set<int> starts;
                    std::set<int> ends;
                    for (int32_t k = 0; k < combination.size(); ++k) {
                        if (allCombinations[k][combination[k]] > -1) { // this is a real overlap
                            starts.insert(allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1());
                            ends.insert(allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1());
                        }
                    }

                    int32_t refStart = 0;
                    int32_t refEnd = 0;
                    int32_t length = 0; // get the longest CNS that fit the criteria
                    for (int32_t start : starts) {
                        for (int32_t end : ends) {
                            std::set<std::string> species;
                            int32_t thisNumberOfSequences = 0;
                            int32_t newLength = end - start + 1;
                            if (length < newLength && mini_cns_size <= newLength) {
                                for (int32_t k = 0; k < combination.size(); ++k) {
                                    if (allCombinations[k][combination[k]] > -1) {
                                        if (allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1() <=
                                            start &&
                                            allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1() >=
                                            end) {
                                            ++thisNumberOfSequences;
                                            species.insert(allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getSpecies());
                                        }
                                    }
                                }
                                if (species.size() >= minimumNumberOfSpecies && thisNumberOfSequences >=
                                                                                minimumNumberOfSpecies /*&& thisNumberOfSpecies > numberOfSpecies*/ ) {
                                    refStart = start;
                                    refEnd = end;
                                    length = newLength;
                                }
                            }
                        }
                    }
                    // find a good refStart and refEnd combination end

                    if ((refEnd > refStart) && ((refEnd - refStart + 1) >= mini_cns_size)) {
                        int32_t newRefStart = refEnd; // since the one used for refEnd or refStart maybe not incldued in the final output, so recheck the refStart and refEnd here
                        int32_t newRedEnd = refStart;
                        std::set<std::string> species;
                        std::set<int> badIndex;

                        for (int32_t k = 0; k < combination.size(); ++k) {
                            if (allCombinations[k][combination[k]] > -1) { // this is a real overlap
                                int32_t thisStart =
                                        allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1() >
                                        refStart
                                        ? allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1()
                                        : refStart;
                                int32_t thisEnd =
                                        allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1() <
                                        refEnd
                                        ? allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1()
                                        : refEnd;
                                double thisLength = thisEnd - thisStart + 1;
                                //                    std::cout << " line 216 thisLength: " << thisLength << " length " << length << std::endl;
                                if (thisLength < outputWithMinimumLengthPercentage * double(length)) {
                                    badIndex.insert(k);
                                } else {
                                    if (thisStart < newRefStart) {
                                        newRefStart = thisStart;
                                    }
                                    if (thisEnd > newRedEnd) {
                                        newRedEnd = thisEnd;
                                    }
                                    species.insert(allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getSpecies());
                                }
                            }
                        }

                        if (newRefStart > refStart) {
                            refStart = newRefStart;
                        }
                        if (newRedEnd < refEnd) {
                            refEnd = newRedEnd;
                        }
                        // std::cout << " line 226  " << i << " species.size " << species.size() << " refStart " << refStart << " refEnd " << refEnd << std::endl;
                        //allPairedSimilarFragments[0] never output
                        if (species.size() >= minimumNumberOfSpecies && (refEnd > refStart) &&
                            ((refEnd - refStart + 1) >= mini_cns_size)) {
                            std::vector<std::string> alignmentNames;
                            std::vector<std::string> alignmentSeqs;

                            //                std::cout << " line 307" << std::endl;
                            for (int32_t k = 0; k < combination.size(); ++k) {
                                if (allCombinations[k][combination[k]] > -1 && badIndex.find(k) == badIndex.end()) {
                                    int32_t start2 = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart2();
                                    int32_t end2 = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd2();

                                    //int32_t endReference = allPairedSimilarFragments[ai][i].getEnd1();
                                    int32_t endReference = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd1();
                                    int32_t positionReference = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart1();
                                    --positionReference;

                                    std::string alignSeq;
                                    std::string refSeq;
                                    int32_t numberOfAltChar = 0;
                                    bool alignRegionStart = false;
                                    bool alignRegionEnd = true;
                                    for (int32_t w = 0; w <
                                                    allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment1().length(); ++w) {
                                        if (allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment1()[w] !=
                                            '-') {
                                            ++positionReference;
                                        }
                                        if (!alignRegionStart) {
                                            if (positionReference > refStart && (!alignRegionStart)) {
                                                int32_t h = positionReference - refStart;
                                                while (h > 0) {
                                                    alignSeq += ".";
                                                    refSeq += referenceGenome[referenceChr][positionReference - h - 1];
                                                    --h;
                                                }
                                            }
                                            if (positionReference >= refStart) {
                                                alignRegionStart = true;
                                                start2 += numberOfAltChar;
                                                end2 = start2;
                                            }
                                        }
                                        if (alignRegionStart && alignRegionEnd) {
                                            alignSeq += allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment2()[w];
                                            refSeq += allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment1()[w];
                                        }
                                        if (allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getAlignment2()[w] !=
                                            '-') {
                                            if (alignRegionStart && alignRegionEnd) {
                                                ++end2;
                                            }
                                            ++numberOfAltChar;
                                        }
                                        if (alignRegionEnd) {
                                            if (positionReference == endReference && positionReference < refEnd) {
                                                int32_t h = refEnd - positionReference;
                                                for (int32_t z = 0; z < h; ++z) {
                                                    alignSeq += ".";
                                                    refSeq += referenceGenome[referenceChr][positionReference + z];
                                                }
                                                alignRegionEnd = false;
                                            }
                                            if (positionReference >= refEnd) {
                                                alignRegionEnd = false;
                                            }
                                        }
                                    }
                                    --end2;

                                    if (alignmentSeqs.empty()) {
                                        alignmentNames.push_back("reference:\t" + referenceChr + ":" + std::to_string(refStart) + "-" +
                                                                 std::to_string(refEnd)); //the reference one
                                        alignmentSeqs.push_back(refSeq);
                                    } else { //update the alignment sequence by comparing the reference sequence
                                        int32_t lengthq = refSeq.length();
                                        for (int32_t o = 0; o < lengthq; ++o) {
                                            if (alignmentSeqs[0].length() > o) { // this is always the reference one
                                                if (refSeq[o] == '-' && alignmentSeqs[0][o] != '-') {
                                                    for (int32_t p = 0; p < alignmentSeqs.size(); ++p) {
                                                        alignmentSeqs[p].insert(o, "-");
                                                    }
                                                } else if (refSeq[o] != '-' && alignmentSeqs[0][o] == '-') {
                                                    refSeq.insert(o, "-");
                                                    alignSeq.insert(o, "-");
                                                    ++lengthq;
                                                }
                                            }
                                        }
                                        while (refSeq.length() < alignmentSeqs[0].length()) {
                                            refSeq.insert(refSeq.length(), "-");
                                            alignSeq.insert(alignSeq.length(), "-");
                                        }
                                        while (refSeq.length() > alignmentSeqs[0].length()) {
                                            for (int32_t p = 0; p < alignmentSeqs.size(); ++p) {
                                                alignmentSeqs[p].insert(alignmentSeqs[p].length(), "-");
                                            }
                                        }
                                    }
                                    if( allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStrand() == 1 ){
                                        alignmentNames.push_back(allPairedSimilarFragments[k][0].getSpecies() + ":\t" + allPairedSimilarFragments[k][0].getQueryChr() + ":" + std::to_string(start2) + "-" +
                                                                 std::to_string(end2) +"\t+");

                                    }else{
                                        int32_t eee = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd2() - (start2 - allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart2());
                                        int32_t sss = allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getStart2() + (allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getEnd2() - end2);
                                        alignmentNames.push_back(allPairedSimilarFragments[k][0].getSpecies() + ":\t" + allPairedSimilarFragments[k][allCombinations[k][combination[k]]].getQueryChr() + ":" + std::to_string(sss) + "-" +
                                                                 std::to_string(eee) +"\t-");
                                    }
                                    alignmentSeqs.push_back(alignSeq);
                                }
                            }
                            //                        ofile << c ;
                            for (int32_t o = 0; o < alignmentNames.size(); ++o) {
                                ofile << ">" << alignmentNames[o] << std::endl << alignmentSeqs[o] << std::endl;
                            }
                            ofile << std::endl << std::endl;
                        }
                    }
                }
            }
        }
    }
    ofile.close();
}
