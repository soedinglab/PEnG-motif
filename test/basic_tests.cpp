//
// Created by Christian Roth on 23.09.17.
//

#include "gtest/gtest.h"
#include "../src/shared/Alphabet.h"
#include "../src/shared/SequenceSet.h"

namespace {
  class BasicTests : public ::testing::Test {
  protected:

    BasicTests() {
      std::string alphabet_string = "STANDARD";
      Alphabet::init(alphabet_string.c_str());
      std::string default_sequences_fasta = "test/test_data/default_sequence_set.fa";

      std::string empty_string = "";
      default_sequence_set = new SequenceSet{default_sequences_fasta.c_str(), true, empty_string.c_str()};
    }

    virtual ~BasicTests() {
      delete default_sequence_set;
    }

    virtual void SetUp() {}

    virtual void TearDown() {}

    SequenceSet* default_sequence_set;
    const unsigned NUMBER_OF_SEQUENCES = 3;
  };

  TEST_F(BasicTests, check_read_sequences) {
    ASSERT_EQ(default_sequence_set->getN(), NUMBER_OF_SEQUENCES);
  }
}