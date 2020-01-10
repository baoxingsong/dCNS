#include <iostream>

#include <stdlib.h>
#include "controlLayer.h"
#include "./googletest/googletest/include/gtest/gtest.h"


int main( int argc, char** argv ) {

    if( argc<=1 ){
        usage();
        return 1;
    }

    std::string program = argv[1];
    if( program.compare("-h") == 0 || program.compare("--help") == 0 ){
        usage();
        exit(1);
    }
/*
    if( program.compare("pair-cns") == 0 ) {
        return pairCns(--argc, ++argv);
    }else if ( program.compare("mult-cns") ==0){
        return multCns(--argc, ++argv);
    }else{
        usage();
    }
    */

    if( program.compare("pairCnsXExtend") == 0 ) {
        return pairCnsXExtend(--argc, ++argv);
    }else if ( program.compare("pairCns2Gaps") ==0){
        return pairCns2Gaps(--argc, ++argv);
    }else if( program.compare("weighted1Gap") == 0 ) {
        return weighted1Gap(--argc, ++argv);
    }else if ( program.compare("weighted2Gaps") ==0){
        return weighted2Gaps(--argc, ++argv);
    }else if( program.compare("cut1Gap") == 0 ) {
        return cut1Gap(--argc, ++argv);
    }else if ( program.compare("cut2Gaps") ==0){
        return cut2Gaps(--argc, ++argv);
    }else if ( program.compare("cut2Gaps2") ==0){
        return cut2Gaps2(--argc, ++argv);
    }else if ( program.compare("maskGenome") ==0 ){
        return maskGenome(--argc, ++argv);
    }else if( program.compare("ranSco") ==0 ){
        return smitherWaterManScoreOfRandomFragments(--argc, ++argv);
    }else if( program.compare("multCns") ==0 ){
        return multCns( --argc, ++argv );
    }else if( program.compare("slideWindow") ==0 ){
        return slideWindow( --argc, ++argv );
    }else if( program.compare("mapCNSToGenome") ==0 ){
        return mapCNSToGenome( --argc, ++argv );
    }else{
        usage();
    }
    return 0;
}

int mainwe( int argc, char** argv ) {
    testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    return 0;
}

// in case it is needed
// in blast the default values for lambda and K
//
// Lambda	1.33271	1.28
// K	0.620991	0.46
// H	1.12409	0.85

//TODO for those failed extension alignment result, should try to use a larger band width