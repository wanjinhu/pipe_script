#! /bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate lefse
lefse_format_input.py $1 $2 -c 1 -s -1 -u 2 -o 1000000
lefse_run.py $2 $3 -l $7
lefse_plot_res.py $3 $4 --format pdf --width 15 --height 35
lefse_plot_cladogram.py $3 $5 --format pdf --clade_sep 0.05
mkdir $6
lefse_plot_features.py $2 $3 $6/ --format pdf --title_font_size 8
# cladogram的结果用AI调整下画布大小

## demo
# conda activate lefse
# lefse_format_input.py lefse_input.txt lefse_input.in -c 1 -s -1 -u 2 -o 1000000
# lefse_run.py lefse_input.in lefse_input.res -l 2.0
# lefse_plot_res.py lefse_input.res lefse_diff.pdf --format pdf --width 15 --height 35
# lefse_plot_cladogram.py lefse_input.res lefse_cladogram.pdf --format pdf --clade_sep 0.05
# mkdir biomarkers_raw_images
# lefse_plot_features.py lefse_input.in lefse_input.res biomarkers_raw_images/ --format pdf --title_font_size 8