//
// Created by Christian Roth on 23.09.17.
//

#include "gtest/gtest.h"
#include <math.h>

namespace {
  TEST(BasicTests, logf_safety) {
  ASSERT_FLOAT_EQ(sqrt(1024.0f), sqrt(1024.0));
  }
}