```shell
# Create libnanopolish.a
cd ../../tools/nanopolish
make
# Compile run_align
cd ../../long-reads/nanopolish
make
./run_align -p data/align_parameter_nanopolish -r data/align_result_nanopolish
```