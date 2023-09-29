#include "host_data_io.h"
#include "host_data.h"

void skip_to_EOR(FILE *fp) {
    const char *loc = "EOR";
    while (*loc != '\0') {
        if (fgetc(fp) == *loc) {
            loc++;
        }
    }
}

call_t read_call(FILE *fp) {
    call_t call;

    long long n;
    float avg_qspan;
    int max_dist_x, max_dist_y, bw, n_segs;

    int t = fscanf(fp, "%lld%f%d%d%d%d",
            &n, &avg_qspan, &max_dist_x, &max_dist_y, &bw, &n_segs);
    // fprintf(stderr, "read %d arguments\n", t);
    if (t != 6) {
        call.n = ANCHOR_NULL;
        call.avg_qspan = .0;
        return call;
    }

    call.n = n;
    call.avg_qspan = avg_qspan;
    call.max_dist_x = max_dist_x;
    call.max_dist_y = max_dist_y;
    call.bw = bw;
    call.n_segs = n_segs;
    // fprintf(stderr, "%lld\t%f\t%d\t%d\t%d\t%d\n", n, avg_qspan, max_dist_x, max_dist_y, bw, n_segs);

    // call.anchors.resize(call.n);

    // for (anchor_idx_t i = 0; i < call.n; i++) {
    //     uint64_t x, y;
    //     fscanf(fp, "%llu%llu", &x, &y);

    //     anchor_t t;
    //     t.x = x; t.y = y;

    //     call.anchors[i] = t;
    // }

    // Some extra space for vectorization with intrinsics.
    call.anchors_x.resize(call.n + 64);
    call.anchors_x32.resize(call.n + 64);

    call.anchors_y.resize(call.n + 64);
    call.anchors_y32.resize(call.n + 64);

    call.q_spans.resize(call.n + 64);

    // Add padding before and after the data.
    for (anchor_idx_t i = 32; i < call.n + 32; i++) {
        uint64_t x, y;
        fscanf(fp, "%llu%llu", &x, &y);

        call.anchors_x[i] = x;
        call.anchors_x32[i] = x;

        call.anchors_y[i] = y;
        call.anchors_y32[i] = y;

        call.q_spans[i] = y >> 32 & 0xff;
    }

    skip_to_EOR(fp);
    return call;
}

void print_return(FILE *fp, const return_t &data)
{
    fprintf(fp, "%lld\n", (long long)data.n);
    // Keep padding in mind (32 before and 32 after)
    for (anchor_idx_t i = 32; i < data.n + 32; i++) {
        fprintf(fp, "%d\t%d\n", (int)data.scores[i], (int)data.parents[i]);
    }
    fprintf(fp, "EOR\n");
}
