//
// Created by bs674 on 8/9/19.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>

TEST(samFileToSeeds, c1){
    std::string samFile = "/media/bs674/2t/testPan_cns/realStory/sorghum/5.sam";
    std::map<std::string, std::string>genome;
    std::vector<std::string> seqNames;
    readFastaFile( "/media/bs674/2t/testPan_cns/realStory/maskGenomeForGenomeAlignment/masked_sorghum_k20_100.fa", genome, seqNames);

    std::map<std::string, std::string> ref_genome;
    std::vector<std::string> ref_seqNames;
    readFastaFile( "/media/bs674/2t/testPan_cns/realStory/maskGenomeForGenomeAlignment/masked_B73_v4_k20_150.fa", ref_genome, ref_seqNames);


    std::map<std::string, std::map<std::string, std::vector<Seed>>> positiveSeeds;
    std::map<std::string, std::map<std::string, std::vector<Seed>>> negativeSeeds;

    samFileToSeeds ( samFile, positiveSeeds, negativeSeeds, ref_genome, genome );

    /*for ( std::map<std::string, std::map<std::string, std::vector<Seed>>>::iterator it = positiveSeeds.begin();
            it!=positiveSeeds.end(); ++it ){
        for ( std::map<std::string, std::vector<Seed>>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
            for( Seed seed : it2->second ){
                std::cout << it->first << ":" << seed.getStart1() << "-" << seed.getEnd1() << "\t" << it2->first << ":" << seed.getStart2() << "-" << seed.getEnd2() << std::endl;
            }
        }
    }*/

    for ( std::map<std::string, std::map<std::string, std::vector<Seed>>>::iterator it = negativeSeeds.begin();
          it!=negativeSeeds.end(); ++it ){
        for ( std::map<std::string, std::vector<Seed>>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
            for( Seed seed : it2->second ){
                std::cout << it->first << ":" << seed.getStart1() << "-" << seed.getEnd1() << "\t" << it2->first << ":" << seed.getStart2() << "-" << seed.getEnd2() << std::endl;

            }
        }
    }

    ASSERT_EQ(0, 0);
}
