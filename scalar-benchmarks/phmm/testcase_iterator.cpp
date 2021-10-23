#include <istream>
#include <string>
#include <vector>
#include <limits>

#include <iostream>

#include "testcase_iterator.h"
#include "read.h"
#include "haplotype.h"

using namespace std;

inline bool end_of_file_reached(const pair<size_t, size_t>& p) {
  return p.first == 0;
}

TestcaseIterator::TestcaseIterator() :
  input_stream { nullptr },
  current      { }
{}

TestcaseIterator::TestcaseIterator(istream* const input) :
  input_stream { input },
  current      { fetch_next_element() }
{}

TestcaseIterator::TestcaseIterator(TestcaseIterator&& original) :
  input_stream { original.input_stream },
  current {move(original.current) }
{
  original.input_stream = nullptr;
}

bool TestcaseIterator::operator!=(const TestcaseIterator& rhs) {
  return input_stream != rhs.input_stream;
}

Testcase& TestcaseIterator::operator*() {
  return current;
}

Testcase& TestcaseIterator::operator++() {
  current = fetch_next_element();
  return current;
}

Testcase TestcaseIterator::fetch_next_element() {
  const auto test_size = parse_test_size();
  if (end_of_file_reached(test_size)) {
    input_stream = nullptr;
    return Testcase{};
  }
  auto reads = parse_all_reads(test_size.first);
  auto haplotypes = parse_all_haplotypes(test_size.second);
  return Testcase {move(reads), move(haplotypes)};
}

pair<size_t, size_t> TestcaseIterator::parse_test_size() {
  int n_reads, n_haplotypes;
  if (! (*input_stream >> n_reads) )
    return make_pair(0,0);
  *input_stream >> n_haplotypes; 
  return make_pair(n_reads, n_haplotypes);
}

Read<uint8_t,uint8_t> TestcaseIterator::parse_read() {
  const auto read_bases      = parse_string();
  const auto read_quals      = parse_string();
  const auto read_ins_quals  = parse_string();
  const auto read_del_quals  = parse_string();
  const auto read_gcp_quals  = parse_string();
  return Read<uint8_t,uint8_t>{read_bases, read_quals, read_ins_quals, read_del_quals, read_gcp_quals};
}

Haplotype<uint8_t> TestcaseIterator::parse_haplotype() {
  const auto haplotype_bases = parse_string();
  return Haplotype<uint8_t>{haplotype_bases};
}

string TestcaseIterator::parse_string() {
  string s;
  *input_stream >> s;
  return s;
}

vector<Read<uint8_t,uint8_t>> TestcaseIterator::parse_all_reads(size_t num_reads) {
  auto reads = vector<Read<uint8_t,uint8_t>>{};
  reads.reserve(num_reads);
  while(num_reads-- > 0) 
    reads.push_back(parse_read()); 
  return reads;
}

vector<Haplotype<uint8_t>> TestcaseIterator::parse_all_haplotypes(size_t num_haplotypes) {
  auto haplotypes = vector<Haplotype<uint8_t>>{};
  haplotypes.reserve(num_haplotypes);
  while(num_haplotypes-- > 0) 
    haplotypes.push_back(parse_haplotype());
  return haplotypes;
}
