# `ftclust`
Functional Tight Clustering

Dependencies: `git`, `sbt`, `java >= 8`, `bash`

1. First download the software:

   ```bash
   $ git clone --branch ftclust https://github.com/scalation/fda.git
   ```

2. Place the software package in an appropriate directory, for instance:

   ```bash
   /home/fda$
   ```

   If you downloaded the software package to your `/Downloads/` directory you can also run the software directly from there:

   ```bash
   /Downloads/fda$
   ```

   **NOTE:** From this point forward, this directory will be abbreviated simply as `$` in the instructions.

3. **[PREPARING]** To prepare the software, navigate to the `scalation_mathstat` directory with the following command:

   ```bash
   $ cd scalation_1.3/scalation_mathstat
   ```
   
   From here type the following command:

   ```bash
   $ sbt publishLocal
   ```

   Next, navigate to the `scalation1_3/scalation_modeling` directory with the following command:

   ```bash
   $ cd ../scalation_modeling
   ```

4. **[RUNNING]** The software is run through the `bash` script from this directory with the following command: 

   ```bash
   $ ./gene_analysis.sh <INFILE> <OUTFILE> [OPTIONS]
   ```

   An `INFILE` and an `OUTFILE` argument must be provided. Please read further for instructions on the `INFILE` format as well as the various options you can pass. 

## File Format

The data you wish to cluster should be placed in a CSV file in the same directory as the `gene_analysis.sh` file. 
This is the `INFILE`.
Each row of the input file should contain observation information about a single sample.
Each column of a row should be a separate observation for that particular sample.
The input file must have a header row labeling the observations (columns) as well as a leading column labeling the samples (rows).
Both the header row and the leading label column may be left blank, but they must at least exist.
In which case the leading label column is left blank, it will be filled in sequentially for purposes of the dataoutput files.

## Software Output

The output files will be placed in a sibling directory to the directory containing the bash script: 

```bash
scalation1_3/scalation_modeling/data/
```

Depending on your choice of options, there will be up to 9 output files. 


## Example

As a simple example for illustrative purposes, you may run the provided simulated dataset, `simu.csv`, as follows:

   ```bash
   $ cd scalation_1.3/scalation/modeling
   $ ./gene_analysis.sh simu.csv simuOut -r \-100
   ```

## Software Options

The following is a comprehensive list of options for running the software:

| Option | Description |
| --- | --- |
| `-r <OPTARG>` | Specify the row sum as `OPTARG` (for filtering the results). DEFAULT is 100. This filters out "housekeeping" or inconsequential samples whose gene expression levels across all observations do not warrant scrutiny. NOTE that a negative number may be passed here as a parameter, however you must use an escape sequence to pass it (i.e. - to pass -1 you must use \-1). |
| `-m <OPTARG>` | Specify k-max as `<OPTARG>`, the starting point for `k0` in the tight clustering algorithm. DEFAULT is `6`. |
| `-s`          | Do not to use Singular Value Decomposition (SVD) for Gap Statistic clustering. DEFAULT is to use SVD. |
| `-p`          | Generate plots for raw and smoothed data points. DEFAULT is to not generate plots. |
| `-k <OPTARG>` | Specify k-est as as `OPTARG`, the estimated number of clusters to try to find in the data for Tight Clustering. DEFAULT is 1. |
| `-a <OPTARG>` | Specify the alpha level as `OPTARG` for Tight Clustering. DEFAULT is `.2`. | 
| `-B <OPTARG>` | Specify the beta level as `OPTARG` for Tight Clustering. DEFAULT is `.9`. |
| `-q <OPTARG>` | Specify top-can as as `OPTARG` the number of top-candidates to chose for each `k0` value for Tight Clustering. DEFAULT is `7`. | 
| `-b <OPTARG>` | Specify `resamp.num` as as `OPTARG` the number of times to resample for Tight Clustering. DEFAULT is `10`. |
| `-w`          | Do not cluster the raw data. DEFAULT is to cluster raw data. |  
| `-t`          | Do not cluster smoothed data. DEFAULT is to cluster smoothed data. |
| `-c`          | Do not cluster the data by coefficients. DEFAULT is to cluster by coefficients. |
| `-l`          | Do not cluster the data with normal k-means++ clustering (i.e., non-tight or loose clustering). DEFAULT is to cluster loose. | 
| `-g`          | Do not cluster the data with tight clustering. DEFAULT is to cluster tight. |
| `-x`          | Cluster the data using the Gap statistic. DEFAULT is to not cluster with the Gap statistic. BEWARE - Gap statistic runs in super exponential time. |
| `-R <OPTARG>` | Specify `samp.p` as `<OPTARG>` the sample ratio for repeated subsamplings in Tight Clustering. DEFAULT is `.7`. |
| `-L <OPTARG>` | Specify `resamp.num` as `<OPTARG>` the number of sequential `k0` to try to find tight and stable clusters. DEFAULT is `3`.
| `-h`          | Ask for help. |

>>>>>>> 0cdd2e83ac02af25eade09e95bf3f14838d76d48
