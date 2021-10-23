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

    call.anchors.resize(call.n);

    for (anchor_idx_t i = 0; i < call.n; i++) {
        uint64_t x, y;
        fscanf(fp, "%llu%llu", &x, &y);

        anchor_t t;
        t.x = x; t.y = y;

        call.anchors[i] = t;
    }

    skip_to_EOR(fp);
    return call;
}

void print_return(FILE *fp, const return_t &data)
{
    fprintf(fp, "%lld\n", (long long)data.n);
    for (anchor_idx_t i = 0; i < data.n; i++) {
        fprintf(fp, "%d\t%d\n", (int)data.scores[i], (int)data.parents[i]);
    }
    fprintf(fp, "EOR\n");
}
