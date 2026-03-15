# MD7115-NTU

Teaching materials for `MD7115: R for Graph, Network and Trees` at NTU Singapore.

## Course overview

This course introduces graph, network, and tree concepts with a focus on biomedical research applications. The class is designed for a 3-hour session and combines conceptual explanations with hands-on R examples.

Main topics:

- graph, network, and tree concepts
- biomedical examples of network analysis
- network visualization in R
- tree visualization in biology
- interpretation, common pitfalls, and take-home messages

## Course structure

The teaching materials are organized into five parts.

1. `Part 1: Introduction and setup` (15-20 min)
2. `Part 2: What are graph, network, and tree?` (25-30 min)
3. `Part 3: Network analysis in R` (60-70 min)
4. `Part 4: Tree visualization in biology` (35-40 min)
5. `Part 5: Interpretation, pitfalls, and take-home messages` (10-15 min)

## Main slide files

The new modular course slides are in `3_materials/`:

- `MD7115-Part1-Introduction-Setup.Rmd`
- `MD7115-Part2-Graph-Network-Tree.Rmd`
- `MD7115-Part3-Network-Analysis-R.Rmd`
- `MD7115-Part4-Tree-Visualization-Biology.Rmd`
- `MD7115-Part5-Interpretation-Pitfalls.Rmd`

The original full slide deck is also retained:

- `MD7115-Xiaotao.Rmd`

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
│   └── personalized_score.csv
├── 3_materials/
│   ├── MD7115-Part1-Introduction-Setup.Rmd
│   ├── MD7115-Part2-Graph-Network-Tree.Rmd
│   ├── MD7115-Part3-Network-Analysis-R.Rmd
│   ├── MD7115-Part4-Tree-Visualization-Biology.Rmd
│   ├── MD7115-Part5-Interpretation-Pitfalls.Rmd
│   └── MD7115-Xiaotao.Rmd
└── MD7115-NTU.Rproj
```

## Before class

Students should complete the following steps before attending class:

1. Install R and RStudio
2. Download or clone this repository
3. Open `MD7115-NTU.Rproj`
4. Run `1_code/0_install_packages.R`

## Code and data

The course uses local data files stored in `2_demo_data/` so that class examples can run without depending on external downloads during teaching.

The code files in `1_code/` provide standalone examples that can be used before, during, or after class.

## Rendering slides

To render a slide deck in R:

```r
rmarkdown::render("3_materials/MD7115-Part1-Introduction-Setup.Rmd")
```

Replace the file name with any of the other `.Rmd` files in `3_materials/`.

## Contact

Xiaotao Shen  
Assistant Professor  
Lee Kong Chian School of Medicine, NTU Singapore  
Website: http://shen-lab.org  
Email: xiaotao.shen@ntu.edu.sg
