# Cancer prognosis with hiPathia

This repo contains data and R codes for reproducing numerical results of the following paper:

> Yunlong Jiao, Marta Hidalgo, Cankut Cubuk, Alicia Amadoz, Jose Carbonell-Caballero, Jean-Philippe Vert, and Joaquin Dopazo. "Signaling Pathway Activities Improve Prognosis for Breast Cancer." bioRxiv preprint bioRxiv-132357, 2017. [bioRxiv-132357](https://doi.org/10.1101/132357).

## Quick start

See the notebook in `results/notebook.[md|html]` for the detailed pipeline, codes and results for the numerical experiments of this study.

## Directories

The top level structure is as follows:

* `data/` - static RData to be used in running experiments, including
  - `path.genes.vals.RData`, `other.genes.vals.pt[1|2].RData` are gene-level profiles of breast tumors and `surv.grps.RData` is donor vital outcome, all downloaded from TCGA-ICGC data portal [release No.20](https://dcc.icgc.org/releases/release_20/Projects/BRCA-US) and further processed as described in paper;
  - `eff.vals.RData`, `path.vals.RData` are pathway-level profiles, processed from gene expression with pathway analysis tool [hiPathia](http://hipathia.babelomics.org/);
  - `fun.vals.RData`, `go.vals.RData` are function-level profiles, processed from pathway activities with [UniProt](http://www.uniprot.org/) or [GO](http://www.geneontology.org/) annotations;
  - `fpgs.RData` contains detailed info of KEGG pathways modeled in paper.
* `results/` - results of numerical experiments, including
  - `notebook.html` is the project notebook with raw codes found in `notebook.Rmd`, markdown version in `notebook.md` and figures saved in `notebook_files/`;
  - `results.scores.txt` contain evaluation scores of prediction performance and `results.pathways.txt`, `results.othergenes.txt` contain selected top features in each type of profile;
  - `runPredict.R` is the R script for running using different profiles to make prediction with an example shell script to submit parallelized tasks to SGE cluster in `runPredict.sh` on parameters defined in `runPredict.param.txt`.
* `src/` - source code and general purpose scripts, including
  - `func.R` implements general functions and classifiers for the entire study.

## Hands-on

To build the project notebook `results/notebook.html` locally, first make sure your local machine has installed the following R packages (or run the corresponding commands to install)

```r
> require(rmarkdown)    # install.packages('rmarkdown')
> require(knitr)        # install.packages('knitr')
> require(devtools)     # install.packages('devtools')
> require(ggplot2)      # install.packages('ggplot2')
> require(igraph)       # install.packages('igraph')
> require(org.Hs.eg.db) # source('https://bioconductor.org/biocLite.R'); biocLite('org.Hs.eg.db')
> require(GO.db)        # source('https://bioconductor.org/biocLite.R'); biocLite('GO.db')
```

then run in shell, which should take seconds to get the project notebook `results/notebook.html`

```sh
$ git clone git@github.com:YunlongJiao/hipathiaCancerPrognosis.git
$ cd hipathiaCancerPrognosis/results/
$ Rscript -e "rmarkdown::render('notebook.Rmd', output_format = 'all')"
```

Note that in order to build the project from scratch, one needs to run the predictions to obtain the three results files `results/results.[scores|pathways|othergenes].txt`, which requires to run `runPredict.R` before building the project notebook. See `results/notebook.html` for detail.

## Built with

* [hiPathia](https://github.com/babelomics/hipathia) - Signaling pathway model

## Authors

* **Yunlong Jiao** - main contributor
* **Marta Hidalgo** - data processing with hiPathia

