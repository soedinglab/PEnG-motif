//
// Created by Christian Roth on 27.09.17.
//

#include "gtest/gtest.h"
#include "../src/base_pattern.h"
#include "../src/iupac_pattern.h"


namespace {
  class BasePatternTest : public ::testing::Test {
  protected:

    BasePatternTest() {
      std::string alphabet_string = "STANDARD";
      Alphabet::init(alphabet_string.c_str());
      std::string default_sequences_fasta = "test/test_data/default_sequence_set.fa";
      std::string empty_string = "";
      default_sequence_set = new SequenceSet{default_sequences_fasta.c_str(), true, empty_string.c_str()};
      default_background_model = new BackgroundModel{*default_sequence_set, 2, Global::bgModelAlpha, Global::interpolateBG};
    }

    virtual ~BasePatternTest() {
      delete default_sequence_set;
    }

    virtual void SetUp() {}

    virtual void TearDown() {}

    SequenceSet* default_sequence_set;
    BackgroundModel* default_background_model;
  };

  TEST_F(BasePatternTest, check_kmer_extension_right) {
    BasePattern* basepattern = new BasePattern{4, Strand::BOTH_STRANDS, 2, 2, default_sequence_set, default_background_model};
    size_t pattern_CT = 1 + 3*4;
    size_t pattern_CTG = 1 + 3*4 + 2*16;
    size_t pattern_CTGA = 1 + 3*4 + 2*16 + 0*64;

    size_t exp_CTG = basepattern->add_letter_to_the_right(pattern_CT, 2, 2);
    ASSERT_EQ(exp_CTG, pattern_CTG);

    size_t exp_CTGA = basepattern->add_letter_to_the_right(pattern_CTG, 3, 0);
    ASSERT_EQ(exp_CTGA, pattern_CTGA);
    delete basepattern;
  }

  TEST_F(BasePatternTest, check_bg_kmer_conversion) {
    BasePattern* basepattern = new BasePattern{4, Strand::BOTH_STRANDS, 2, 2, default_sequence_set, default_background_model};
    // careful: we are converting between different kmer representations
    size_t pattern_CTGA = 1 + 3*4 + 2*16 + 0*64;
    size_t bg_kmer = 3*16 + 2*4 + 0;

    size_t extracted_kmer = basepattern->get_bg_id(pattern_CTGA, 4, 2);
    ASSERT_EQ(extracted_kmer, bg_kmer);
    delete basepattern;
  }

  TEST_F(BasePatternTest, iupac2base_patterns) {
    size_t pattern_length = 4;
    BasePattern* basepattern = new BasePattern{pattern_length, Strand::BOTH_STRANDS, 2, 2, default_sequence_set, default_background_model};
    size_t pattern_CTRA = 1*1 + 3*11 + 6*11*11 + 0*11*11*11;
    size_t pattern_CTAA =  1*1 + 3*4 + 0*16 + 0*64;
    size_t pattern_CTGA =  1*1 + 3*4 + 2*16 + 0*64;
    IUPACPattern* iupac = new IUPACPattern{pattern_CTRA, pattern_length};
    iupac->init(pattern_length, default_background_model->getV()[2]);

    auto results = iupac->generate_base_patterns(basepattern, pattern_CTRA);

    ASSERT_EQ(results.size(), 2);
    ASSERT_EQ(results[0], pattern_CTAA);
    ASSERT_EQ(results[1], pattern_CTGA);


    size_t pattern_NTRN = 10*1 + 3*11 + 6*11*11 + 10*11*11*11;

    std::vector<size_t> old_patterns;
    IUPACPattern::find_base_patterns(basepattern, pattern_NTRN, pattern_length, old_patterns);
    auto new_result = iupac->generate_base_patterns(basepattern, pattern_NTRN);
    ASSERT_EQ(old_patterns.size(), new_result.size());
    for(auto pattern : new_result) {
      ASSERT_TRUE(std::find(old_patterns.begin(), old_patterns.end(), pattern) != old_patterns.end());
    }

    size_t pattern_YCCT = 7 + 1*11 + 1*11*11 + 3*11*11*11;
    size_t pattern_CCCT = 1 + 1*4 + 1*16 + 3*64;
    size_t pattern_TCCT = 3 + 1*4 + 1*16 + 3*64;

    auto YCCT_result = iupac->generate_base_patterns(basepattern, pattern_YCCT);
    ASSERT_EQ(YCCT_result.size(), 2);
    ASSERT_TRUE(std::find(YCCT_result.begin(), YCCT_result.end(), pattern_CCCT) != YCCT_result.end());
    ASSERT_TRUE(std::find(YCCT_result.begin(), YCCT_result.end(), pattern_TCCT) != YCCT_result.end());

    delete basepattern;
    delete iupac;
  }
}