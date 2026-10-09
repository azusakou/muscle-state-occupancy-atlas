# Reproducibility package

Code and processed inputs for the manuscript analyses and five figures.

## Environment and execution

Use Python 3.12, R 4.6.1 and Bioconductor 3.23. Python requirements are listed
in `requirements.txt`; tested direct R package versions are in `R-packages.tsv`.
Install these dependencies in your own environment, then run from this directory:

```sh
python3 -m pip install -r requirements.txt
python3 run.py
```

`Rscript` must be available on PATH; alternatively set `RSCRIPT` to its executable.
For a separate R library, set `R_LIBS_USER`. Results and run logs are written to
`results/`; the five figures are written to `figures/`. The workflow verifies
release checksums and numerical results and stops on failure. Software versions
can affect numerical agreement; R session information is recorded with results.

## Contents and scope

`code/` contains analyses and validation; `data/` contains the required inputs.
`DATA_PROVENANCE.json` records data provenance. Model settings and evidence tiers are supplied in the corresponding TSV files.

Reproduction starts from processed inputs. Three lineage DESeq2 contrasts are
refitted from supplied counts. Spatial and denervation gene-level statistics,
fixed-threshold annotation counts and NMJ scores are processed inputs; their upstream fitting,
raw-read alignment, full quality control and annotation training are outside
this package. Individual annotation scores are omitted, so those filtering
thresholds cannot be changed. Unused metadata and all-zero lineage genes are
omitted. No precomputed output figures or analysis results are bundled.

`NA` denotes unavailable age, a standard deviation with only one baseline, or
an undefined/filtered adjusted P value; it must not be replaced with zero.
