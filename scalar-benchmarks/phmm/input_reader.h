#ifndef INPUT_READER_H
#define INPUT_READER_H

#include <iostream>
#include <fstream>

template <typename Iterator>
class InputReader {
  std::ifstream file_stream;             /// < a file stream (in case we are reading from a file)
  std::istream* input_stream{&std::cin}; /// < a pointer to the input stream

 public:
  
  InputReader() = default;
  ~InputReader() { file_stream.close(); }
  InputReader(const InputReader&) = delete;

  Iterator begin() { return Iterator{input_stream}; }
  Iterator end()   { return Iterator{}; }

  void from_file(const std::string& filename) {
    file_stream.open(filename);
    input_stream = &file_stream;
  }
};

#endif
