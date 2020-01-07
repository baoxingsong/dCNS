//
// Created by bs674 on 8/23/19.
//

#ifndef AND_CNS_FASTAINDEXENTRY_H
#define AND_CNS_FASTAINDEXENTRY_H

#include <string>

class FastaIndexEntry {
private:
    std::string name; // sequence name
    int32_t length; // length of sequence
    int64_t offset; // bytes offset of sequence from start of file
    int32_t line_blen; // line length in bytes, sequence characters
    int32_t line_len;
public:
    FastaIndexEntry(const std::string &name, int32_t length, int64_t offset, int32_t lineBlen, int32_t lineLen);
    FastaIndexEntry();

    const std::string &getName() const;

    void setName(const std::string &name);

    int32_t getLength() const;

    void setLength(int32_t length);

    int64_t getOffset() const;

    void setOffset(int64_t offset);

    int32_t getLineBlen() const;

    void setLineBlen(int32_t lineBlen);

    int32_t getLineLen() const;

    void setLineLen(int32_t lineLen);

};


#endif //AND_CNS_FASTAINDEXENTRY_H
