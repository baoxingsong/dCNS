#!python
import re
import sys
import subprocess
from argparse import ArgumentParser
# bs674@cornell.edu


def runSamtools(bedfile, bamFile, outputFile):
    output = open(outputFile, 'w')
    outputFileDepth = outputFile + "Depth"
    outputFileMpileup = outputFile + "Mpileup"
    outputDepth = open(outputFileDepth + "", 'w')
    outputMpileup = open(outputFileMpileup + "Mpileup", 'w')
    for i in range(1, 11):
        chr = str(i)
        with open(bedfile) as f:
            for line in f:
                m = re.search('^'+ chr +'\s+(\d+)\s+(\d+)\s+(\S+)\s', line)
                if (m != None):
                    start = int(m.group(1))+1
                    end = int(m.group(2))+1
                    length = end - start + 1
                    region = chr + ':' + str(start) + "-" + str(end)
                    numberOfMpileupCoveredBases = 0
                    numberOfDepthCoveredBases = 0
                    for line2 in subprocess.run(['samtools', 'mpileup', '-r', region, bamFile], stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = 'utf8').stdout.split("\n"):
                        m2 = re.search('^' + chr + '\s+(\d+)\s+\S+\s+(\d+)', line2)
                        if (m2 != None):
                            if int(m2.group(2)) > 0:
                                numberOfMpileupCoveredBases = numberOfMpileupCoveredBases + 1
                                outputMpileup.write(chr + "\t" + m2.group(1) + "\n")
                    for line2 in subprocess.run(['samtools', 'depth', '-r', region, bamFile], stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = 'utf8').stdout.split("\n"):
                        m2 = re.search('^'+chr+'\s+(\d+)\s+(\d+)', line2)
                        if (m2 != None):
                            if int(m2.group(2)) > 0:
                                numberOfDepthCoveredBases = numberOfDepthCoveredBases + 1
                                outputDepth.write(chr + "\t" + m2.group(1) + "\n")
                    output.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" +str(length) + "\t" + str(numberOfMpileupCoveredBases) + "\t" + str(numberOfDepthCoveredBases) + "\n")

    output.close()
    outputDepth.close()
    outputMpileup.close()


if __name__ == '__main__':
    parser = ArgumentParser(description='check the base resolution CNS alignment coverage')
    parser.add_argument("-b", "--bam",
                        dest="bam",
                        type=str,
                        default="",
                        help="bam file")
    parser.add_argument("-e", "--bed",
                        dest="bed",
                        type=str,
                        default="",
                        help="bed file")
    parser.add_argument("-o", "--output",
                        dest="output",
                        type=str,
                        default="",
                        help="output file")


    args = parser.parse_args()

    if args.bam == "":
        print("Error: please specify --bam", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.bed == "":
        print("Error: please specify --bed", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.output == "":
        print("Error: please specify --output", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    runSamtools(args.bed, args.bam, args.output)
