# MD7115-NTU

Teaching materials for `MD7115: R for Graph, Network and Trees` at NTU Singapore.

## Overview

This repository contains slides, example scripts, and demo datasets for a 3-hour teaching session on graph, network, and tree analysis in R, with applications in biomedical research.

The course covers:

- basic concepts of graphs, networks, and trees
- biological examples of network structures
- network analysis and visualization in R
- tree visualization for microbiome and related biological data
- interpretation, limitations, and common pitfalls

## Main teaching material

The primary slide deck is:

- `3_materials/MD7115-Xiaotao.Rmd`

The rendered slide deck is:

- `3_materials/MD7115-Xiaotao.html`

## Repository structure

```text
MD7115-NTU/
├── 1_code/
│   ├── 0_install_packages.R
│   ├── 1_example1.R
│   ├── 2_example2.R
│   └── 3_example.R
├── 2_demo_data/
│   ├── example_edge_data.csv
│   ├── example_node_data.csv
│   ├── microbiome_data.RData
│   ├── pathway_result.csv
│   ├── pathway_result.rda
│   └── personalized_score.csv
├── 3_materials/
│   ├── MD7115-Xiaotao.Rmd
│   ├── MD7115-Xiaotao.html
│   ├── images/
│   ├── libs/
│   ├── my-theme.css
│   └── xaringan-themer.css
└── MD7115-NTU.Rproj
```

## Files for students

Before class, students should mainly use:

- `1_code/0_install_packages.R` to install required packages
- `MD7115-NTU.Rproj` to open the project in RStudio
- `3_materials/MD7115-Xiaotao.html` or `3_materials/MD7115-Xiaotao.Rmd` for the teaching slides

## Example scripts

The `1_code/` folder contains scripts used for demonstration and follow-up practice:

- `1_example1.R`: network visualization example
- `2_example2.R`: GO similarity network example
- `3_example.R`: tree visualization example

## Demo data

The `2_demo_data/` folder contains local data used in class so that examples can be run without relying on internet access during teaching.

## Before class

Recommended preparation:

1. Install R
2. Install RStudio
3. Download or clone this repository
4. Open `MD7115-NTU.Rproj`
5. Run `1_code/0_install_packages.R`

## Rendering slides

To render the main slide deck in R:

```r
rmarkdown::render("3_materials/MD7115-Xiaotao.Rmd")
```

## Contact

Xiaotao Shen  
Assistant Professor  
Lee Kong Chian School of Medicine, NTU Singapore  
Website: http://shen-lab.org  
Email: xiaotao.shen@ntu.edu.sg
