To generate config files, at first the required parameters are written in sweep.json.
Next, the code can be run in the terminal with the command python parameter_sweep.py \
    --base-config ../inputs/run_config.txt \
    --sweep sweep.json \
    --benchmark-root benchmark_runs
This will produce the folder benchmark_runs

Currently, slurm script is still under construction.
    