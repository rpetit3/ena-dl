[![Anaconda-Server Badge](https://anaconda.org/bioconda/ena-dl/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ena-dl/badges/downloads.svg)](https://anaconda.org/bioconda/ena-dl)

# Please use fastq-dl
ena-dl has run its course, while it should still does what it says it does (haha at least at the time of writing this), it has been replaced by [fastq-dl](https://github.com/rpetit3/fastq-dl). fastq-dl has the same support as ena-dl but also includes the ability to download from SRA.

# ena-dl
ena-dl, short for *'ENA download'*, is a tool for downloading FASTQ files from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena). ENA was chosen over the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) due to the availability of FASTQ files without the need for conversion (e.g. SRA format to FASTQ).

ena-dl takes an ENA/SRA accession (Study, Experiment, or Run) and queries ENA (via [Data Warehouse API](https://www.ebi.ac.uk/ena/browse/search-rest)) to determine the Aspera/FTP FASTQ paths. These FASTQ files are then downloaded for each Run. For Samples or Experiments with multiple Runs, users can optionally merge the runs.

### Alternatives
ena-dl was originally developed for [Staphopia](https://staphopia.emory.edu) and is now being used for [Bactopia](https://github.com/bactopia/bactopia), as you might imagine some of the functions are more geared towards these projects. There are currently no plans to add any other functions to ena-dl. Although don't let that stop you from using it, it still does what it says it does!

If you are looking for more rounded tools for downloading FASTQs, may I suggest [enaBrowserTools](https://github.com/enasequence/enaBrowserTools) or [SRA Toolkit](https://github.com/ncbi/sra-tools). These tools are built and maintained by EBI/NCBI and provide are more extensive access to their databases.

# Installation
### Bioconda
ena-dl is available from [Bioconda](https://bioconda.github.io/) and I highly recommend you go this route to for installation.
```
conda install -c conda-forge -c bioconda ena-dl
```

### From Source
However, being a single Python script, with a single dependency, ena-dl is pretty straight forward to set up from source.
```
git clone git@github.com:rpetit3/ena-dl.git
cd ena-dl
pip3 install -r requirements.txt
python3 ena-dl.py
```

Feel free to drop it in your $PATH somewhere!

# Usage
*ena-dl* requires a single ENA/SRA Study, Experiment, or Run accession and FASTQs for all Runs that fall under the given accession will be downloaded. For example, if a Study accession is given all Runs under that studies umbrella will be downloaded.

### Usage Output
```
usage: ena-dl [-h] [--aspera STRING] [--aspera_key STRING]
              [--aspera_speed STRING] [--is_study] [--is_experiment]
              [--is_run] [--group_by_experiment] [--group_by_sample]
              [--outdir OUTPUT_DIR] [--max_retry INT] [--ftp_only] [--silent]
              [--debug] [--version]
              ACCESSION

ena-dl - Download FASTQs from ENA

optional arguments:
  -h, --help            show this help message and exit

Required Options:

  ACCESSION             ENA accession to query. (Study, Experiment, or Run
                        accession)

Aspera Connect Options:
  --aspera STRING       Path to the Aspera Connect tool "ascp" (Default:
                        "which ascp")
  --aspera_key STRING   Path to Aspera Connect private key, if not given,
                        guess based on ascp path
  --aspera_speed STRING
                        Speed at which Aspera Connect will download. (Default:
                        100M)

Query Related Options:
  --is_study            Query is a Study.
  --is_experiment       Query is an Experiment.
  --is_run              Query is a Run.
  --group_by_experiment
                        Group Runs by experiment accession.
  --group_by_sample     Group Runs by sample accession.

Helpful Options:
  --outdir OUTPUT_DIR   Directory to output downloads to. (Default: ./)
  --max_retry INT       Maximum times to retry downloads (Default: 10)
  --ftp_only            FTP only downloads.
  --silent              Only critical errors will be printed.
  --debug               Skip downloads, print what will be downloaded.
  --version             show program's version number and exit
```

### Example Usage
#### Download a Study
```
ena-dl.py PRJNA248678
```

The above command would download 3 runs that fall under Study PRJNA248678. The relationship of Study to Run is a 1-to-many relationship, or there can be many Run accessions associated with a single Study Accession. You could use `--group_by_experiment` to group these runs by Experiment accession (or `--group_by_sample` for Sample accession).

#### Download an Experiment
```
ena-dl.py SRX477044
```

The above command would download the single run that fall under Experiment SRX477044. The relationship of Experiment to Run is a 1-to-many relationship, or there can be many Run accessions associated with a single Experiment Accession (e.g. resequencing the same sample). Although in most cases, especially for bacterial samples, it is a 1-to-1 relationship. In any case, you can use `--group_by_experiment` to merge multiple runs associated with an Experiment accession into a single FASTQ file (or `--group_by_sample` for Sample accession).

#### Download a Run
```
ena-dl.py SRR1178105
```

The above command would download the Run SRR1178105. Run accessions are the end of the line (1-to-1 relationship)

#### Using Aspera Connect
Users can use [Aspera Connect](https://downloads.asperasoft.com/connect2/) to speed up the download of FASTQs from ENA. Installation and setup of Aspera Connect is out of the scope of this documentation, but I can assure you its a simple installation.

```
ena-dl.py SRR1178105 --aspera /path/to/ascp
```

The above command will attempt to download SRR1178105 using the `ascp` tool. By default it will try to use the private key file that is included during the Aspera Connect installation. If it is not found you will need to use the `--aspera_key` parameter to specify its path.
