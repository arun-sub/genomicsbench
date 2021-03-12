First, go to tools/Clair and create envirnment with instructions.

Second, come back and run:

python3 my_prediction.py --chkpnt_fn ./model --sampleName HG001 --threads 1 --qual 100 --input_fn <prediction_input.h5> --output_fn <prediction_output.h5>