#ifndef HOST_INPUT_H
#define HOST_INPUT_H

#include <vector>
#include <cstdint>

typedef int64_t anchor_idx_t;
typedef uint32_t tag_t;
typedef int32_t loc_t;
typedef int32_t loc_dist_t;
typedef int32_t score_t;
typedef int32_t parent_t;
typedef int32_t target_t;
typedef int32_t peak_score_t;

#define ANCHOR_NULL (anchor_idx_t)(-1)

// struct anchor_t {
//     uint64_t x;
//     uint64_t y;
// };

// struct call_t {
//     anchor_idx_t n;
//     float avg_qspan;
//     int max_dist_x, max_dist_y, bw, n_segs;
//     std::vector<anchor_t> anchors;
// };

struct call_t {
    anchor_idx_t n;
    float avg_qspan;
    int max_dist_x, max_dist_y, bw, n_segs;
    std::vector<uint64_t> anchors_x;
    std::vector<uint32_t> anchors_x32; // WRONG in some cases.

    std::vector<uint64_t> anchors_y;
    std::vector<uint32_t> anchors_y32; // WRONG in some cases.

    std::vector<int32_t> q_spans;
};

struct return_t {
    anchor_idx_t n;
    std::vector<score_t> scores;
    std::vector<parent_t> parents;
    std::vector<target_t> targets;
    std::vector<peak_score_t> peak_scores;
};

#endif // HOST_INPUT_H
