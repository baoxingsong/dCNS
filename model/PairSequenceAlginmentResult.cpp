//
// Created by bs674 on 5/22/19.
//

#include "PairSequenceAlginmentResult.h"

PairSequenceAlginmentResult::PairSequenceAlginmentResult(const std::string &refName, const std::string &altName,
                                                         const std::string &refSeq, const std::string &altSeq)
        : refName(refName), altName(altName), refSeq(refSeq), altSeq(altSeq) {}

const std::string &PairSequenceAlginmentResult::getRefName() const {
    return refName;
}

void PairSequenceAlginmentResult::setRefName(const std::string &refName) {
    PairSequenceAlginmentResult::refName = refName;
}

const std::string &PairSequenceAlginmentResult::getAltName() const {
    return altName;
}

void PairSequenceAlginmentResult::setAltName(const std::string &altName) {
    PairSequenceAlginmentResult::altName = altName;
}

const std::string &PairSequenceAlginmentResult::getRefSeq() const {
    return refSeq;
}

void PairSequenceAlginmentResult::setRefSeq(const std::string &refSeq) {
    PairSequenceAlginmentResult::refSeq = refSeq;
}

const std::string &PairSequenceAlginmentResult::getAltSeq() const {
    return altSeq;
}

void PairSequenceAlginmentResult::setAltSeq(const std::string &altSeq) {
    PairSequenceAlginmentResult::altSeq = altSeq;
}

void PairSequenceAlginmentResult::addDel(const int & position){
    this->refSeq.insert(position, "-");
    this->altSeq.insert(position, "-");
}
