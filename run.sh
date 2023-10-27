# snakemake --dag  | dot -Tsvg > dag.svg
nohup snakemake -s Snakefile --cores 50 --latency-wait 100 &
