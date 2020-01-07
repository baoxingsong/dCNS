//
// Created by bs674 on 5/22/19.
//

#ifndef SONG_CNS_PAIRSEQUENCEALGINMENTRESULT_H
#define SONG_CNS_PAIRSEQUENCEALGINMENTRESULT_H


#include <string>

class PairSequenceAlginmentResult {
private:
    std::string refName;
    std::string altName;
    std::string refSeq;
    std::string altSeq;
public:
    PairSequenceAlginmentResult(const std::string &refName, const std::string &altName, const std::string &refSeq,
                                const std::string &altSeq);

    const std::string &getRefName() const;

    void setRefName(const std::string &refName);

    const std::string &getAltName() const;

    void setAltName(const std::string &altName);

    const std::string &getRefSeq() const;

    void setRefSeq(const std::string &refSeq);

    const std::string &getAltSeq() const;

    void setAltSeq(const std::string &altSeq);

    void addDel(const int & position);
};


#endif //SONG_CNS_PAIRSEQUENCEALGINMENTRESULT_H
