#ifndef READ_H
#define READ_H

#include <ostream>
#include <vector>
#include <string>

#include "utils.h"

template<class QUAL_TYPE, class BASE_TYPE>
struct Read {
  size_t original_length, left_padding, right_padding;
  std::vector<BASE_TYPE> bases;
  std::vector<QUAL_TYPE> base_quals;
  std::vector<QUAL_TYPE> ins_quals;
  std::vector<QUAL_TYPE> del_quals;
  std::vector<QUAL_TYPE> gcp_quals;

  explicit Read() = default;

  Read(const std::vector<QUAL_TYPE>& bases_,
       const std::vector<QUAL_TYPE>& quals_,
       const std::vector<QUAL_TYPE>& ins_quals_,
       const std::vector<QUAL_TYPE>& del_quals_,
       const std::vector<QUAL_TYPE>& gcp_quals_) : 
    original_length{bases_.size()}, 
    left_padding{0}, 
    right_padding{0},
    bases      {bases_},
    base_quals {quals_},
    ins_quals  {ins_quals_},
    del_quals  {del_quals_},
    gcp_quals  {gcp_quals_}
  {}

  Read(const std::string& bases_,
       const std::string& base_quals_,
       const std::string& ins_quals_,
       const std::string& del_quals_,
       const std::string& gcp_quals_,
       const size_t original_length_ = 0) :
    original_length {original_length_},
    bases     {convert_bytes<std::vector<BASE_TYPE>>(bases_)},
    base_quals{convert_bytes<std::vector<QUAL_TYPE>>(base_quals_, -QUAL_OFFSET)},
    ins_quals {convert_bytes<std::vector<QUAL_TYPE>>(ins_quals_,  -QUAL_OFFSET)},
    del_quals {convert_bytes<std::vector<QUAL_TYPE>>(del_quals_,  -QUAL_OFFSET)},
    gcp_quals {convert_bytes<std::vector<QUAL_TYPE>>(gcp_quals_,  -QUAL_OFFSET)}
  {}

};

template<class QUAL_TYPE, class BASE_TYPE>
std::ostream& operator<<(std::ostream& out, const Read<QUAL_TYPE, BASE_TYPE>& read) {
  out << convert_bytes<std::string>(read.bases) << " ";
  out << convert_bytes<std::string>(read.base_quals, QUAL_OFFSET) << " ";
  out << convert_bytes<std::string>(read.ins_quals, QUAL_OFFSET) << " ";
  out << convert_bytes<std::string>(read.del_quals, QUAL_OFFSET) << " ";
  out << convert_bytes<std::string>(read.gcp_quals, QUAL_OFFSET) << " ";
  return out;
}

#endif
