/*
 * main.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: mmeier
 */

#include "log.h"
#include "Global.h"
#include "PatternCensus.h"

int main(int nargs, char **args) {
  Global::init(nargs, args);

  BackgroundModel* bgModel = new BackgroundModel(*Global::backgroundSequenceSet,
                          Global::bgModelOrder,
                          Global::bgModelAlpha,
                          Global::interpolateBG );

  PatternCensus(Global::patternLength, Global::bgModelOrder, Global::inputSequenceSet, bgModel);

  delete bgModel;
}
