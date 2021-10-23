#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

static constexpr uint8_t QUAL_OFFSET = 33;

template<typename Result, typename Input>
static inline Result convert_bytes(const Input& s, const uint8_t offset = 0) {
  auto result = Result {};
  result.reserve(s.size());
  for (int i = 0; i != s.size(); ++i)
    result.push_back(s[i] + offset);
  return result;
}

#endif

