# fda
Functional Data Analysis Group

Dependencies: git, sbt, java 8 jdk, bash shell (and with it getopts)

First download the software:

   >$ git clone --branch ftclust https://github.com/scalation/fda.git

Place the software package in an appropriate directory, for instance:

   >/home/fda$

If you downloaded the software package to your /Downloads/ directory you can also run the software directly from there:

   >/Downloads/fda$

To run the software navigate to the scalation_mathstat directory with the following command:

   >/home/fda$ cd scalation_1.3/scalation_mathstat

   >/home/fda/scalation_1_3/scalation_mathstat$
   
From here type the following command:

   >/home/fda/scalation_1_3/scalation_mathstat$ sbt publishLocal

Next, navigate to the /home/fda/scalation1_3/scalation_modeling directory with the following command:

   >/home/fda/scalation_1_3/scalation_mathstat$ cd ../scalation_modeling

   >/home/fda/scalation_1_3/scalation_modeling$ 

The software is run through the bash script from this directory with the following command: 

   >/home/fda/scalation_1_3/scalation_modeling$ ./gene_analysis.sh INFILE OUTFILE [OPTIONS]

An INFILE and an OUTFILE argument must be provided. Please read further for instructions on the INFILE format as well as the various options you can pass. 

The data you wish to cluster should be placed in a CSV file in the same directory as the gene_analysis.sh file. This is the INFILE.
Each row of the input file should contain observation information about a single tissue sample.
Each column of a row should be a separate observation for that particular tissue sample.
The input file must have a header row labeling the observations (columns) as well as a leading column labeling the tissue samples (rows).
Both the header row and the leading label column may be left blank, but they must at least exist.
In which case the leading label column is left blank, it will be filled in sequentially for purposes of the dataoutput files.

The output files will be placed in a sibling directory to the directory containing the bash script: 

   >/home/fda/scalation1_3/scalation_modeling/data/

Depending on your choice of options, there will be up to 9 output files. 

As a simple example for illustrative purposes, you may run the provided simulated dataset, 'simu.csv', as follows:

   >/home/fda/scalation_1.3/scalation/modeling$ ./gene_analysis.sh simu.csv simuOut -r \-10

The following is a comprehensive list of options for running the software:

	-r <OPTARG>
		specify the row sum as OPTARG (for filtering the results). DEFAULT is 100. This filters out "housekeeping" or inconsequential samples
		whose gene expression levels across all observations do not warrant scrutiny. NOTE that a negative number may be passed here as a
		parameter, however you must use an escape sequence to pass it (i.e. - to pass -1 you must use \-1)
	-m <OPTARG>
		specify k-max as <OPTARG>, the starting point for k0 in the tight clustering algorithm. DEFAULT is 6
	-s
		do not to use SVD for Gap Statistic clustering. DEFAULT is to use SVD
	-p
		plot the data points raw and smoothed data. DEFAULT is not plots. 
	-k <OPTARG>
		specify k-est as as OPTARG, the estimated number of clusters to try to find in the data for Tight Clustering. DEFAULT is 1.  
	-a <OPTARG>
		specify the alpha level as OPTARG for Tight Clustering. DEFAULT is .2
	-B <OPTARG>
		specify the beta level as OPTARG for Tight Clustering. DEFAULT is .9
	-q <OPTARG>
		specify top-can as as OPTARG the number of top-candidates to chose for each k0 value for Tight Clustering. DEFAULT is 7 
	-b <OPTARG>
		specify resamp.num as as OPTARG the number of times to resample for Tight Clustering. resamp.num. DEFAULT is 10. 
	-w
		do not cluster the raw data. DEFAULT is to cluster raw data. 
	-t
		do not cluster smoothed data. DEFAULT is to cluster smoothed data. 
	-c
		do not cluster the data by coefficients. DEFAULT is to cluster by coefficients. 
	-l
		do not cluster the data with normal KMeans++ clustering (i.e.-non-tight or loose clustering). DEFAULT is to cluster loose. 
	-g
		do not cluster the data with tight clustering. DEFAULT is to cluster tight. 
	-x
		cluster the data using the Gap statistic. DEFAULT is to not cluster with the Gap statistic. BEWARE - Gap statistic runs in super exponential time.
	-R <OPTARG>
		specify samp.p as <OPTARG> the sample ratio for repeated subsamplings in Tight Clustering. DEFAULT is .7
	-L <OPTARG>
		specify resamp.num as <OPTARG> the number of sequential k0 to try to find tight and stable clusters. DEFAULT is 3
	-h
		ask for help. 
	
