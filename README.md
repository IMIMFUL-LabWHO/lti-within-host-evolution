# lti-within-host-evolution

Code, processed data, and workflow files accompanying the study of longitudinal SARS-CoV-2 within-host evolution during prolonged infection.

This repository supports two complementary levels of reproducibility:

| Level | Starting point | Purpose |
| --- | --- | --- |
| **Level 1: processed-data analysis** | Supplied metadata and variant tables | Reproduce the statistical analyses and figures reported in the manuscript without reprocessing raw sequencing reads. |
| **Level 2: full sequencing pipeline** | Raw FASTQ files | Reprocess the sequencing reads with LTIseek, regenerate the variant dataset, and then repeat the statistical analyses. |

## Repository contents

| Path | Description |
| --- | --- |
| `data/metadata.csv` | Sample- and patient-level metadata used in the statistical analyses. |
| `data/variants.csv` | Processed SARS-CoV-2 variant dataset used in the statistical analyses. |
| `statistical_analysis/lti-within-host-evolution.R` | R script used to reproduce the statistical analyses and manuscript figures. |
| `scripts/download_fastq.sh` | Script for downloading the raw FASTQ files. |
| `LTIseek/20250313_lti.def` | Singularity definition file for building the analysis container. |
| `LTIseek/config.sh` | LTIseek input paths and pipeline settings. |
| `LTIseek/lti_main.sh` | Main LTIseek workflow entry point. |
| `LTIseek/reference` | LTIseek reference genome, GFF file and primer info. |

## Requirements

### Level 1

- [Git](https://git-scm.com/)
- [R](https://www.r-project.org/)
- The R packages loaded by `statistical_analysis/lti-within-host-evolution.R`

### Level 2

- All Level 1 requirements
- A Linux environment with [SingularityCE](https://docs.sylabs.io/guides/latest/admin-guide/installation.html)
- Sufficient storage for the raw FASTQ files and pipeline outputs
- Internet access for downloading the sequencing data and building the container

## Get the repository

Clone the repository and enter its root directory:

```bash
git clone https://github.com/IMIMFUL-LabWHO/lti-within-host-evolution.git
cd lti-within-host-evolution
```

The commands below assume that they are run from the repository root.

## Level 1: reproduce the statistical analyses and figures

Level 1 starts from the processed files supplied in `data/` and does not require the raw sequencing reads.

Run the analysis script from the repository root:

```bash
Rscript statistical_analysis/lti-within-host-evolution.R
```

Alternatively, open the script in RStudio or another R environment, set the repository root as the working directory, and run it interactively. If R reports a missing package, install that package and rerun the script.

## Level 2: reproduce the full sequencing workflow

Level 2 starts from the raw FASTQ files, processes them with LTIseek, regenerates the variant dataset, and then repeats the Level 1 analysis.

### 1. Install SingularityCE

Install SingularityCE by following the [official installation guide](https://docs.sylabs.io/guides/latest/admin-guide/installation.html). On a shared computing system, SingularityCE may already be available as an environment module.

### 2. Build the container

From the repository root, build `lti.sif` from the supplied definition file.

With administrator privileges:

```bash
sudo singularity build lti.sif 20250313_lti.def
```

Without administrator privileges, use `--fakeroot` if it is enabled and configured on the system:

```bash
singularity build --fakeroot lti.sif 20250313_lti.def
```

These are alternative build commands; only one is required.

### 3. Download the raw FASTQ files

Assuming that `lti.sif` was built in the repository root, run:

```bash
singularity exec "$PWD/lti.sif" bash scripts/download_fastq.sh
```

Raw-read availability and any applicable access conditions are described in the Data Availability statement of the accompanying paper.

### 4. Configure LTIseek

Before starting the workflow, inspect `config.sh` and verify:

- the path to the `fastq_data` directory;
- input and output paths; and
- all other pipeline settings relevant to the local computing environment.

Absolute paths are recommended, particularly when running the workflow on a high-performance computing system.

### 5. Run LTIseek

```bash
singularity exec "$PWD/lti.sif" bash lti_main.sh
```

If the container was built elsewhere, replace `"$PWD/lti.sif"` with its absolute path.

### 6. Repeat the statistical analysis

After LTIseek completes:

1. confirm that the regenerated variant table is available at the location expected by `statistical_analysis/lti-within-host-evolution.R`;
2. retain the supplied `data/variants.csv` if you want to compare the original and regenerated results; and
3. rerun the analysis:

```bash
Rscript statistical_analysis/lti-within-host-evolution.R
```

## Reproducibility notes

- Use the repository release or Git commit associated with the paper.
- Run commands from the repository root unless a script explicitly states otherwise.
- The Singularity definition file records the software environment used by the sequencing workflow.
- Minor numerical or graphical differences can arise when the R analysis is run with different R or package versions.

## Citation

If you use this repository, its code, or its processed data, please cite the accompanying paper. Full citation details will be added once the article is published.

## Contact

For questions about the repository or workflow, contact [Alen Suljič](mailto:alen.suljic@mf.uni-lj.si).
