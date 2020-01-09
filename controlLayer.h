/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 09:38:17
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#ifndef _CONTROLLAYER_H
#define _CONTROLLAYER_H

#include <iostream>
#include "InputParser.h"
#include <sstream>
#include "./impl/impl.h"
#include "./model/model.h"
#include <iostream>


/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/25/2017 09:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************




 ************************************************************************/




void pairCnsXExtend(std::string & _input,  std::string & _reference, std::string & _output,
                    int & _matchingScore, int & _mismatchingPenalty, int & _openGapPenalty,
                    int & _extendGapPenalty, int32_t & _seed_window_size, int32_t & _mini_cns_score,
                    int32_t & _step_size, int & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
                    const double & _kValue, const int & _w, const int & _xDrop, double & _pvalues);


void pairCns2Gaps(std::string & _input,  std::string & _reference, std::string & _output,
                  int & _matchingScore, int & _mismatchingPenalty, int & _openGapPenalty1,
                  int & _extendGapPenalty1, int & _openGapPenalty2,
                  int & _extendGapPenalty2, int32_t & _seed_window_size, int32_t & _mini_cns_score,
                  int32_t & _step_size, int & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
                  const double & _kValue, const int & _w, const int32_t & bandwidth, const int & _xDrop,
                  const int32_t & _zDrop, double & _pvalues);

void weighted1Gap(std::string & _input,  std::string & _reference, std::string & _output,
                  int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty1,
                  int32_t & _extendGapPenalty1, int32_t & _seed_window_size, int32_t & _mini_cns_score,
                  int32_t & _step_size, int32_t & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
                  const double & _kValue, int32_t & _bandwidth, int32_t & _w, const int32_t & _xDrop, const int32_t & _zDrop, const std::string & scoreFoler,
                  const std::string & refFasta, const std::string & gffFile, double & _pvalues);

void weighted2Gaps(std::string & _input,  std::string & _reference, std::string & _output,
                   int & _matchingScore, int & _mismatchingPenalty, int & _openGapPenalty1,
                   int & _extendGapPenalty1, int32_t & _seed_window_size, int32_t & _mini_cns_score,
                   int32_t & _step_size, int & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
                   const double & _kValue, int & _w, const int32_t & bandwidth, const int & _xDrop, const int32_t & _zDrop, const std::string & scoreFoler,
                   const std::string & refFasta,  const std::string & gffFile, double & _pvalues);

void cut1Gap(std::string & _input,  std::string & _reference, std::string & _output,
             int & _matchingScore, int & _mismatchingPenalty, int & _openGapPenalty1,
             int & _extendGapPenalty1, int32_t & _seed_window_size, int32_t & _mini_cns_score,
             int32_t & _step_size, int & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
             const double & _kValue, int & _w, const int & _xDrop,
             const std::string & refFasta, const std::string & queryFasta, double & _pvalues);

void cut2Gaps(std::string & _input,  std::string & _reference, std::string & _output,
              int & _matchingScore, int & _mismatchingPenalty, int & _openGapPenalty1,
              int & _extendGapPenalty1, int & _openGapPenalty2,
              int & _extendGapPenalty2, int32_t & _seed_window_size, int32_t & _mini_cns_score,
              int32_t & _step_size, int & _matrix_boundary_distance, bool & _onlySyntenic,  const double & _lambda,
              const double & _kValue, int & _w, const int32_t & bandwidth, const int & _xDrop, const int32_t & _zDrop,
              const std::string & refFasta, const std::string & queryFasta, double & _pvalues);
int pairCnsXExtend(int argc, char** argv);
int pairCns2Gaps(int argc, char** argv);
int weighted1Gap(int argc, char** argv);
int weighted2Gaps(int argc, char** argv);
int cut1Gap(int argc, char** argv);
int cut2Gaps(int argc, char** argv);
int cut2Gaps2(int argc, char** argv);
int maskGenome(int argc, char** argv);
int smitherWaterManScoreOfRandomFragments(int argc, char** argv);
int multCns( int argc, char** argv );
int slideWindow( int argc, char** argv );




void mapCNSToGenome(std::string & _input,  std::string & _reference, std::string & _output,
                    int32_t & _matchingScore, int32_t & _mismatchingPenalty, int32_t & _openGapPenalty1,
                    int32_t & _extendGapPenalty1, int32_t & _openGapPenalty2,
                    int32_t & _extendGapPenalty2);
int mapCNSToGenome(int argc, char** argv);
#endif
