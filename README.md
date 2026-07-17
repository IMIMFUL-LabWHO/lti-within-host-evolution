# lti-within-host-evolution
Repository of code and metadata used to analyse the SARS-CoV-2 sequencing data for LTI within-host evolution study.

When in doubt, write alen.suljic@mf.uni-lj.si

Two levels of reproducibilty are provided.

**Level 1 (faster):**

Use the supplied metadata.csv and variants.csv files in data directory to reproduce all statistical analyses and figures in the manuscript using R script lti-within-host-evolution.R, provided in statistical_analysis directory.

**Level 2 (full pipeline):**

Download raw FASTQ files, process them using the included LTIseek workflow, and regenerate the variants dataset before repeating the statistical analyses.

Level 2 workflow steps:

1. Clone the repository in the desired directory
```
git clone https://github.com/IMIMFUL-LabWHO/lti-within-host-evolution.git
```

2. Install Singularity (https://singularity-tutorial.github.io/01-installation/)

3. Build the Singularity container:
```
sudo singularity build lti.sif 20250313_lti.def
singularity build --fakeroot lti.sif 20250313_lti.def #if no sudo privileges
```

4. Download FASTQ files using the download_fastq.sh file provided in scripts directory
```
singularity exec /path/to/lti.sif bash download_fastq.sh
```

5. Run LTIseek (check config.sh file for correct path to fastq_data directory and pipeline settings):
```
singularity exec /path/to/lti.sif bash lti_main.sh
```
