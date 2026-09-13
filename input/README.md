# Input files

The tutorial requires four input files. By default, the command-line script looks for:

```text
input/example_pixy_fst_data.csv
input/example_pixy_pi_data.csv
input/example_comparison_map.txt
input/example_chromosome_map.txt
```

The Fst and π examples are comma-separated. The two map files are tab-separated. The command-line script also accepts tab-separated Fst and π files when their filenames do not end in `.csv`.

## Fst output

Required columns:

```text
pop1  pop2  chromosome  window_pos_1  window_pos_2  avg_wc_fst
```

## π output

Required columns:

```text
pop  chromosome  window_pos_1  window_pos_2  avg_pi
```

## Comparison map

Each row defines one pairwise comparison:

```text
comp  population_1  population_2  plot_order
```

Population names must exactly match the names in the pixy files. `population_1` is the numerator and `population_2` is the denominator of the π ratio. Comparison names must use the format `population_1.v.population_2`.

## Chromosome map

Each row connects a chromosome name in the pixy output to its numeric plot label and order:

```text
chromosome  chromo.num
```

## Custom filenames

Supply custom filenames after the project directory. Relative paths are resolved from the project directory:

```bash
Rscript scripts/meadow_plot_cli.R PROJECT_DIR FST_FILE PI_FILE COMPARISON_MAP_FILE CHROMOSOME_MAP_FILE [MAX_INTERVENING_NONOUTLIERS]
```

For example:

```bash
Rscript scripts/meadow_plot_cli.R /path/to/project input/my_fst.tsv input/my_pi.tsv input/my_comparisons.txt input/my_chromosomes.txt 3
```
