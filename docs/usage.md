# nf-core/chipseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/chipseq/usage](https://nf-co.re/chipseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple replicates

The `sample` identifier should be identical when you have multiple replicates from the same experimental group; just increment the `replicate` identifier appropriately. The first replicate value for any given experimental group must be 1.

The `antibody` column is required to separate the downstream consensus peak merging for different antibodies. It is not advisable to generate a consensus peak set across different antibodies especially if their binding patterns are inherently different e.g. narrow transcription factors and broad histone marks.

The `control` column should be the `sample` identifier for the controls for any given IP. This column together with the `control_replicate` column will set the corresponding control for each of the samples in the table.

```console
group,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,BLA203A1_S27_L006_R1_001.fastq.gz,,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,BLA203A25_S16_L002_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A49_S40_L001_R1_001.fastq.gz,,3,BCATENIN,WT_INPUT,3
WT_INPUT,BLA203A6_S32_L006_R1_001.fastq.gz,,1,,,
WT_INPUT,BLA203A30_S21_L002_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A31_S21_L003_R1_001.fastq.gz,,3,,,
```

### Multiple runs of the same library

Both the `sample` and `replicate` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will perform the alignments in parallel, and subsequently merge them before further analysis. Below is an example where the samples called `WT_BCATENIN_IP` and `WT_INPUT` have been re-sequenced multiple times:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,BLA203A1_S27_L006_R1_001.fastq.gz,,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,BLA203A25_S16_L001_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A25_S16_L002_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A25_S16_L003_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A49_S40_L001_R1_001.fastq.gz,,3,BCATENIN,WT_INPUT,3
WT_INPUT,BLA203A6_S32_L006_R1_001.fastq.gz,,1,,,
WT_INPUT,BLA203A30_S21_L001_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A30_S21_L002_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A31_S21_L003_R1_001.fastq.gz,,3,,,
```

### Full design

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 7 columns to match those defined in the table below.

A final design file may look something like the one below. This is for two antibodies and associated controls, where the second replicate of the `WT_BCATENIN_IP` and `NAIVE_BCATENIN_IP` samples have been sequenced twice:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,BLA203A1_S27_L006_R1_001.fastq.gz,,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,BLA203A25_S16_L001_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A25_S16_L002_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A49_S40_L001_R1_001.fastq.gz,,3,BCATENIN,WT_INPUT,3
NAIVE_BCATENIN_IP,BLA203A7_S60_L001_R1_001.fastq.gz,,1,BCATENIN,NAIVE_INPUT,1
NAIVE_BCATENIN_IP,BLA203A43_S34_L001_R1_001.fastq.gz,,2,BCATENIN,NAIVE_INPUT,2
NAIVE_BCATENIN_IP,BLA203A43_S34_L002_R1_001.fastq.gz,,2,BCATENIN,NAIVE_INPUT,2
NAIVE_BCATENIN_IP,BLA203A64_S55_L001_R1_001.fastq.gz,,3,BCATENIN,NAIVE_INPUT,3
WT_TCF4_IP,BLA203A3_S29_L006_R1_001.fastq.gz,,1,TCF4,WT_INPUT,1
WT_TCF4_IP,BLA203A27_S18_L001_R1_001.fastq.gz,,2,TCF4,WT_INPUT,2
WT_TCF4_IP,BLA203A51_S42_L001_R1_001.fastq.gz,,3,TCF4,WT_INPUT,3
NAIVE_TCF4_IP,BLA203A9_S62_L001_R1_001.fastq.gz,,1,TCF4,NAIVE_INPUT,1
NAIVE_TCF4_IP,BLA203A45_S36_L001_R1_001.fastq.gz,,2,TCF4,NAIVE_INPUT,2
NAIVE_TCF4_IP,BLA203A66_S57_L001_R1_001.fastq.gz,,3,TCF4,NAIVE_INPUT,3
WT_INPUT,BLA203A6_S32_L006_R1_001.fastq.gz,,1,,,
WT_INPUT,BLA203A30_S21_L001_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A31_S21_L003_R1_001.fastq.gz,,3,,,
NAIVE_INPUT,BLA203A12_S3_L001_R1_001.fastq.gz,,1,,,
NAIVE_INPUT,BLA203A48_S39_L001_R1_001.fastq.gz,,2,,,
NAIVE_INPUT,BLA203A49_S1_L006_R1_001.fastq.gz,,3,,,
```

| Column              | Description                                                                                                                                                                            |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`            | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`           | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`           | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `replicate`         | Integer representing replicate number. This will be identical for re-sequenced libraries. Must start from `1..<number of replicates>`.                                                 |
| `antibody`          | Antibody name. This is required to segregate downstream analysis for different antibodies. Required when `control` is specified.                                                       |
| `control`           | Sample name for control sample.                                                                                                                                                        |
| `control_replicate` | Integer representing replicate number for control sample.                                                                                                                              |

Example design files have been provided with the pipeline for [paired-end](../assets/samplesheet_pe.csv) and [single-end](../assets/samplesheet_se.csv) data.

> **NB:** The `group` and `replicate` columns were replaced with a single `sample` column as of v2.0 of the pipeline. The `sample` column is essentially a concatenation of the `group` and `replicate` columns. If all values of `sample` have the same number of underscores, fields defined by these underscore-separated names may be used in the PCA plots produced by the pipeline, to regain the ability to represent different groupings.

## Reference genome files

The minimum reference genome requirements are a FASTA and a GTF file, all other files required to run the pipeline can be generated from these files. However, it is more storage and compute friendly if you are able to re-use reference genome files as efficiently as possible. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices (e.g. those unavailable on [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space. You can then either provide the appropriate reference genome files on the command-line via the appropriate parameters (e.g. `--bwa_index '/path/to/bwa/index/'`) or via a [custom config file](https://nf-co.re/usage/configuration#custom-configuration-files).

- If `--genome` is provided then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`.
- If `--gene_bed` is not provided then it will be generated from the GTF file.

> **NB:** Compressed reference files are also supported by the pipeline i.e. standard files with the `.gz` extension and indices folders with the `tar.gz` extension.

## Blacklist bed files

The blacklist bed files where obtained using the commands below:

```console
cd ..
mkdir -p v1.0
cd v1.0
wget -L https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz && gunzip ENCFF001TDO.bed.gz && mv ENCFF001TDO.bed hg19-blacklist.v1.bed

mkdir -p assets/blacklists/v2.0/
cd assets/blacklists/v2.0/
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.v2.bed.gz && gunzip ce10-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.v2.bed.gz && gunzip ce11-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.v2.bed.gz && gunzip dm3-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz && gunzip dm6-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz && gunzip hg19-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz && gunzip hg38-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz && gunzip mm10-blacklist.v2.bed.gz

cd ..
mkdir -p v3.0
cd v3.0
wget -L https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz && gunzip ENCFF356LFX.bed.gz && mv ENCFF356LFX.bed hg38-blacklist.v3.bed
```

> **NB:** A detailed description of the different versions of the files can be found [here](https://github.com/Boyle-Lab/Blacklist/blob/master/README.md). Also, to to see which blacklist bed files are assigned by default to the respective reference genome check the [igenomes.config](https://github.com/nf-core/chipseq/blob/master/conf/igenomes.config).

### Renaming chromosomes (`--chrom_alias`)

If the genome FASTA uses one naming convention (e.g. RefSeq `NC_060925.1`) but you want every output, and matching with external tracks such as UCSC RepeatMasker/CenSat BigBeds, to use another (e.g. UCSC `chr1`), provide a chromosome alias table:

```bash
nextflow run nf-core/chipseq --fasta GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz --gtf genes.gtf.gz \
    --chrom_alias hs1.chromAlias.txt ...
```

The FASTA, GTF (also when converted from `--gff`), gene BED and blacklist are renamed before chromosome sizes and aligner indices are built, so BAM/bigWig/peak files, IGV sessions and all annotations use the new names. Every name on a line of the table is treated as equivalent, and sequences are renamed to the name in `--chrom_alias_column` (default `ucsc`, the first column of UCSC `<assembly>.chromAlias.txt` files, e.g. `https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chromAlias.txt`; a 1-based column number also works). Sequences missing from the table keep their names; per-file reports are written to `genome/chrom_rename/`, and the renamed reference files to `genome/` with `--save_reference`. A pre-built aligner index cannot be combined with `--chrom_alias` because it contains the original names.

## Optional analyses added to this version

All four additions are off by default.

### Multi-mapped read allocation with Allo

```bash
nextflow run nf-core/chipseq --aligner bowtie2 --with_allo ...
```

Bowtie2 is run with `-k 25` (`--allo_max_alignments`; paired-end data also get `--no-mixed --no-discordant`), the raw alignments are grouped by read name with `samtools collate` and written as SAM, and [Allo](https://github.com/seqcode/allo) is run on them with `--mixed`. Allo keeps the SAM record of the location it chooses unchanged, so an allocated read that Bowtie2 reported as a secondary alignment keeps FLAG 0x100 and would be dropped by Picard, MACS3 and featureCounts; the pipeline therefore clears 0x100 on reads carrying Allo's `ZA`/`ZZ` tags. Every downstream step (merging, duplicate marking, filtering, peak calling, ...) then uses the Allo alignments. Because allocated reads keep Bowtie2's low MAPQ, the BAMTools `-q 1` filter is disabled when `--with_allo` is set.

- **Container / environment.** The public biocontainer (`allo:1.2.0--pyhdfd78af_0`) cannot be used: it lacks the `keras` package that TensorFlow >= 2.16 needs, and Allo fails with `ModuleNotFoundError: No module named 'tensorflow.keras'`. Recipes for a working environment (Allo 1.2.0, TensorFlow 2.21.0, Keras 3.15.1) are shipped in `containers/allo/`:

  ```bash
  # Docker / Podman
  docker build -t allo:1.2.0-keras containers/allo
  nextflow run ... -profile docker --with_allo --allo_container allo:1.2.0-keras

  # Singularity / Apptainer
  singularity build allo-1.2.0-keras.sif containers/allo/allo.def
  nextflow run ... -profile singularity --with_allo --allo_container $PWD/allo-1.2.0-keras.sif
  ```

  With `-profile conda`/`mamba`, `containers/allo/environment.yml` (which adds `keras`) is used automatically. With a container engine and neither `--allo_container` nor `--allo_use_local`, the pipeline stops at start-up with these instructions.
- `--allo_use_local`: no container or conda environment is used for the Allo task; the `allo` executable on the host `$PATH` is run instead (`pip install bio-allo keras`). The task fails with a clear message if `allo` is not found or `keras` is missing.
- `--allo_max_workers` (default `4`): number of parallel Allo workers (`-p`). Every worker loads TensorFlow and the CNN and gets its own copy of the genome-wide map of uniquely mapped reads, so memory grows roughly linearly with this value (about 0.5 GB per worker before any read data, several GB per worker on deep human data). If a worker crashes (`TerminatedWorkerError`, SIGSEGV or SIGKILL), the task exits with 139 and is retried up to twice with half as many workers and more memory. Set `0` to use all task CPUs.
- `--allo_args`: extra Allo arguments, e.g. `'--remove-zeros'`, `'-max 10'` or `'--keep-unmap'`.
- `--allo_random_control`: also pass `--random` for control samples, as suggested by the Allo authors.

Note that Bowtie2 with `-k` is considerably slower, that Allo removes unmapped reads by default (so library-level mapping rates in MultiQC appear close to 100%), and that an "LIB: Allo multi-mapped read allocation" table is added to MultiQC.

### Motif discovery with XSTREME

```bash
nextflow run nf-core/chipseq --run_xstreme [--xstreme_motif_db JASPAR.meme] ...
```

XSTREME is run on the sequences of every per-sample peak set and every consensus peak set with default settings, except `--maxw`, which is set to the length of the shortest read in the final filtered, sorted BAM file (`*.mLb.clN.sorted.bam`) used to call the peaks (for consensus peaks: the shortest across that antibody's samples), read from the samtools stats read-length histogram. STREME, which XSTREME runs, only supports motif widths up to 30, so by default the value is capped at 30 with a warning; use `--xstreme_maxw_cap 0` to pass the read length unchanged. `--xstreme_motif_db` optionally supplies known motifs (`--m`) so discovered motifs are compared with them.

### RepeatMasker / CenSat peak annotation

```bash
# First run: build the BigBeds from source annotations (or with --run_repeatmasker --repeatmasker_species human)
nextflow run nf-core/chipseq --annotate_features \
    --repeatmasker_annotation rmsk.txt.gz --censat_bed censat.bed ...

# Later runs: reuse the BigBeds published in <outdir>/genome/annotation/
nextflow run nf-core/chipseq --annotate_features \
    --repeatmasker_bigbed repeatmasker.bb --censat_bigbed censat.bb ...
```

Sources, in order of precedence: `--repeatmasker_bigbed` > `--repeatmasker_annotation` (RepeatMasker `.out`, UCSC `rmsk.txt[.gz]` or BED) > `--run_repeatmasker` (optionally `--repeatmasker_species`, `--repeatmasker_lib`); `--censat_bigbed` > `--censat_bed`. Remote URLs (e.g. UCSC hub BigBeds) can be passed directly. Chromosome names must match the genome FASTA; features on other sequences are dropped. For BigBeds not built by the pipeline, the feature class is read from an autoSql field named `repClass`/`class`/`type` if present, from RepeatMasker-style names (`AluY#SINE/Alu`), or derived from the name (CenSat: `hsat2_2(...)` -> `hsat2`).

Peaks are reported as covering a feature when they overlap at least `--feature_min_overlap` (default `0.8`, i.e. 80%) of the feature's length.

**Chromosome names** of the RepeatMasker/CenSat features must match the genome FASTA. Exact matches are used first, then an optional alias table (`--feature_chrom_alias`, e.g. UCSC `hs1.chromAlias.txt` to map T2T GenBank/RefSeq accessions such as `CP068277.2` to `chr1`), then `chr` prefix differences (`1` <-> `chr1`) and `chrM` <-> `MT` are fixed automatically. Features on chromosomes that still cannot be matched are dropped, and the run stops with an error if none match. A per-chromosome report (`genome/annotation/<feature_set>.chrom_report.tsv`) shows how each name was matched.

**Expectation for CenSat:** CenSat arrays (HOR, HSat, ...) are typically tens of kb to several Mb long, so a peak rarely covers 80% of one. Empty `*.censat.feature_covered.tsv` tables are therefore expected; use `*.censat.all_overlaps.tsv` (which reports `peak_pct_covered`) to see which peaks fall inside satellite arrays, or lower `--feature_min_overlap`.

### Differential binding and visualisation

```bash
nextflow run nf-core/chipseq --run_differential [--differential_fdr 0.05] [--differential_lfc 1] ...
```

For each antibody with at least two sample groups and replicates, DESeq2 is run on the consensus-peak counts and every pair of groups (sample name without `_REP<n>`) is compared. Results are plotted together with the HOMER annotation of the consensus peaks and, when `--annotate_features` is also set, the RepeatMasker / CenSat classes those peaks cover.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/chipseq --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/chipseq -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/chipseq
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/chipseq releases page](https://github.com/nf-core/chipseq/releases) and find the latest pipeline version - numeric only (eg. `2.1.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.1.0`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
