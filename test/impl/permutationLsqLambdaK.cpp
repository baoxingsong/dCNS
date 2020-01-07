//
// Created by bs674 on 6/8/19.
//

#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>

TEST(permutationLsqLambdaK, c1) {

    std::string referenceFasta = "/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/maize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa";
    std::vector<std::string> queryFastas;
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/Miscanthus_sinensis/Msinensis_497_v7.0.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/Setaria_italica/Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/sugarcane_tareploid/Sspon.HiC_chr_asm.fasta");

    int _open_gap_penalty = -4;
    int _extend_gap_penalty = -2;

    int matchingScore = 2;
    int mismatchingPenalty = -3;

    int32_t seed =1;
    int32_t length = 100000;
    int permutationTimes = 10000;
    bool removen = true;
    permutationLsqLambdaK( referenceFasta, queryFastas, _open_gap_penalty,  _extend_gap_penalty, matchingScore,
            mismatchingPenalty, length, permutationTimes, seed, removen);
}

TEST(permutationLsqLambdaKslow, c2) {
    std::string referenceFasta = "/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/maize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa";
    std::vector<std::string> queryFastas;
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/Miscanthus_sinensis/Msinensis_497_v7.0.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/Setaria_italica/Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/sugarcane_tareploid/Sspon.HiC_chr_asm.fasta");

    int _open_gap_penalty = -3;
    int _extend_gap_penalty = -1;

    int matchingScore = 2;
    int mismatchingPenalty = -3;

    int32_t seed =1;
    int32_t length = 1000;
    int permutationTimes = 10000;
    bool removen = true;
    permutationLsqLambdaK( referenceFasta, queryFastas, _open_gap_penalty,  _extend_gap_penalty, matchingScore,
                           mismatchingPenalty, length, permutationTimes, seed, removen);

}

TEST(permutationLsqLambdaKslow, c3) {
    std::string referenceFasta = "/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/maize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa";
    std::vector<std::string> queryFastas;
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/Miscanthus_sinensis/Msinensis_497_v7.0.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/Setaria_italica/Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa");
    queryFastas.push_back("/media/bs674/1_10t/assembly/core_Andropogoneae_genes/tabasco2.0/genomes/sugarcane_tareploid/Sspon.HiC_chr_asm.fasta");

    int _open_gap_penalty = -10;
    int _extend_gap_penalty = -1;

    int matchingScore = 4;
    int mismatchingPenalty = -5;

    int32_t seed =1;
    int32_t length = 1000;
    int permutationTimes = 10000;
    bool removen = true;
    permutationLsqLambdaK( referenceFasta, queryFastas, _open_gap_penalty,  _extend_gap_penalty, matchingScore,
                           mismatchingPenalty, length, permutationTimes, seed, removen);

}
