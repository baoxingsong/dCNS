//
// Created by bs674 on 8/23/19.
//

#include "FastaIndexEntry.h"

FastaIndexEntry::FastaIndexEntry(const std::string &name, int32_t length, int64_t offset, int32_t lineBlen,
                                 int32_t lineLen) : name(name), length(length), offset(offset), line_blen(lineBlen),
                                                    line_len(lineLen) {}
FastaIndexEntry::FastaIndexEntry(){}

const std::string & FastaIndexEntry::getName() const {
    return name;
}

void FastaIndexEntry::setName(const std::string &name) {
    FastaIndexEntry::name = name;
}

int32_t FastaIndexEntry::getLength() const {
    return length;
}

void FastaIndexEntry::setLength(int32_t length) {
    FastaIndexEntry::length = length;
}

int64_t FastaIndexEntry::getOffset() const {
    return offset;
}

void FastaIndexEntry::setOffset(int64_t offset) {
    FastaIndexEntry::offset = offset;
}

int32_t FastaIndexEntry::getLineBlen() const {
    return line_blen;
}

void FastaIndexEntry::setLineBlen(int32_t lineBlen) {
    line_blen = lineBlen;
}

int32_t FastaIndexEntry::getLineLen() const {
    return line_len;
}

void FastaIndexEntry::setLineLen(int32_t lineLen) {
    line_len = lineLen;
}
