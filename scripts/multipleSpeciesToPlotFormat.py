#!python
import re

# bs674@cornell.edu

class Range:
    def __init__(self, refSeqName, refStart, refEnd, start, end, seqName, id, length):
        self.refSeqName=refSeqName
        self.refStart=refStart
        self.refEnd=refEnd
        self.start = start
        self.end = end
        self.seqName = seqName
        self.id = id
        self.length = length

def readFile(file):
    ranges = []
    referStart = 0
    referEnd = 0
    id = 0
    refSeqName = ""
    length = 0
    with open(file) as f:
        for line in f:
            m = re.search('>SORBI_3009G024600_maize_V4_\\+_8_135910122_136010121:(\d+)\-(\d+)', line)
            if (m != None):
                referStart = int(m.group(1))
                referEnd = int(m.group(2))
                refSeqName = "maize"
                length = referEnd - referStart + 1
            else:
                m = re.search('>SORBI_3009G024600_(.*?):(\d+)\-(\d+)', line)
                if (m != None):
                    id = id + 1
                    seqName = m.group(1)
                    start = int(m.group(2))
                    end = int(m.group(3))
                    range = Range(refSeqName, referStart, referEnd, start, end, seqName, id, length)
                    ranges.append(range)
                    referStart=start
                    referEnd=end
                    refSeqName=seqName
    return ranges


ranges = readFile("/home/bs674/Dropbox/andropogoneae-conservation/SSW/cmake-build-debug/test")
f = open("/home/bs674/maize_sorghum5p_MSA_plot",'w')
for range in ranges:
    f.write(str(range.refStart) + " " + range.refSeqName + "\t" + str(range.id) + "\t" + str(range.length) + "\n")
    f.write(str(range.refEnd) + " " + range.refSeqName + "\t" + str(range.id) + "\t" + str(range.length) + "\n")

    f.write(str(range.start) + " " + range.seqName + "\t" + str(range.id) + "\t" + str(range.length) + "\n")
    f.write(str(range.end) + " " + range.seqName + "\t" + str(range.id) + "\t" + str(range.length) + "\n")

f.close()
