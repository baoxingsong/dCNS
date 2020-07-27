cmake_minimum_required(VERSION 3.0)
project(dCNS)
include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)
set(CMAKE_CXX_STANDARD 14)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
else()
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
endif()



CHECK_CXX_COMPILER_FLAG("-mavx2" CXX_COMPILER_OPT_ARCH_AVX2_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-msse2" CXX_COMPILER_OPT_ARCH_SSE2_SUPPORTED)
CHECK_C_COMPILER_FLAG("-mavx2" C_COMPILER_OPT_ARCH_AVX2_SUPPORTED)
CHECK_C_COMPILER_FLAG("-msse2" C_COMPILER_OPT_ARCH_SSE2_SUPPORTED)
if(CXX_COMPILER_OPT_ARCH_AVX2_SUPPORTED AND C_COMPILER_OPT_ARCH_AVX2_SUPPORTED)
    set(CMAKE_CXX_FLAGS "-Wall -g3 -O3 -funroll-all-loops -mavx2 -I ./ -std=gnu++14 -pthread")
    set(CMAKE_C_FLAGS "-Wall -g3 -O3 -funroll-all-loops -mavx2 -I ./ -std=gnu++14 -pthread")
elseif( CXX_COMPILER_OPT_ARCH_SSE2_SUPPORTED AND C_COMPILER_OPT_ARCH_SSE2_SUPPORTED)
    set(CMAKE_CXX_FLAGS "-Wall -msse2 -g3 -O3 -funroll-all-loops -I ./ -std=gnu++14 -pthread")
    set(CMAKE_C_FLAGS "-Wall -msse2 -g3 -O3 -funroll-all-loops -I ./ -std=gnu++14 -pthread")
else()
    set(CMAKE_CXX_FLAGS "-Wall -g3 -O3 -funroll-all-loops -I ./ -std=gnu++14 -pthread")
    set(CMAKE_C_FLAGS "-Wall -g3 -O3 -funroll-all-loops -I ./ -std=gnu++14 -pthread")
endif()


#for google unit test
add_subdirectory(./googletest/googletest)
include_directories(googletest/googletest/include googletest/googletest)

add_executable(dCNS main.cpp
        controlLayer4.cpp
        controlLayer.h
        InputParser.cpp
        InputParser.h
        impl/sequenceAlignment3.cpp
        impl/sequenceAlignment.h
        impl/fasta.cpp
        impl/fasta.h
        impl/SequenceCharToUInt8.cpp
        impl/SequenceCharToUInt8.h
        impl/impl.h
        impl/findSimilarFragmentsForPairedSequence4.cpp
        impl/findSimilarFragmentsForPairedSequence.h
        model/PairedSimilarFragment.cpp
        model/PairedSimilarFragment.h
        model/model.h
        test/impl/testExtend.cpp
        test/impl/testKmerReading.cpp
        test/impl/getCnsForMultipleSpecies.cpp
        test/impl/SequenceCharToUInt8.cpp
        test/impl/smithWaterman.cpp
        test/impl/findSimilarFragmentsForPairedSequence.cpp
        test/impl/findSimilarFragmentsForPairedSequenceControlLayer4.cpp
        test/impl/caculateK.cpp
        impl/syntenic.cpp
        impl/syntenic.h
        model/Node.cpp
        model/Node.h
        model/Score.cpp
        model/Score.h
        model/Graph.cpp
        model/Graph.h
        impl/longestPath.cpp
        impl/longestPath.h
        impl/gffToCategory.cpp
        impl/gffToCategory.h
        model/CnsMasResult.cpp
        model/CnsMasResult.h
        model/MsaGraph.cpp
        model/MsaGraph.h
        impl/getCnsForMultipleSpecies.cpp
        impl/getCnsForMultipleSpecies.h
        model/PairSequenceAlginmentResult.cpp
        model/PairSequenceAlginmentResult.h
        impl/calculateLambda.cpp
        impl/calculateLambda.h
        impl/calculateK.cpp
        impl/calculateK.h
        impl/permutationLsqLambdaK.cpp
        impl/permutationLsqLambdaK.h
        test/impl/permutationLsqLambdaK.cpp
        test/impl/samFileToSeeds.cpp
        test/impl/readFasta.cpp
        version.h
        impl/slidingWindowAlignment.cpp
        impl/slidingWindowAlignment.h
        test/impl/slidingWindowAlignment.cpp
        model/Matrix.cpp model/Matrix.h
        impl/maskGenome.cpp
        impl/maskGenome.h
        impl/getAllTheCdsPositions.cpp
        impl/getAllTheCdsPositions.h impl/samFileToSeeds.cpp impl/samFileToSeeds.h model/Seed.cpp model/Seed.h test/impl/samFileToSeeds.cpp model/FastaIndexEntry.cpp model/FastaIndexEntry.h impl/FastaIndexEntryImpl.cpp impl/FastaIndexEntryImpl.h impl/FastaWithIndex.cpp impl/FastaWithIndex.h model/utils.cpp model/utils.h impl/mapCNSToGenome.cpp impl/mapCNSToGenome.h)

target_link_libraries(dCNS gtest gtest_main)

