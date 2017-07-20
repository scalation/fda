#!/usr/bin/env bash

USAGE='$0 <infile> <outfile> [options]'

HELP='Help'

if [ "$#" -lt "2" ]; then
    echo "usage: $USAGE"
    exit 1
fi

ARG_1="$1"        #infile
ARG_2="$2"        #outfile
shift 2

mkdir -p data

ARG_3="100"       #rowSum
ARG_4="6"         #k-max (called k-min in the tight clustering class...)
ARG_5="1"         #useSVD
ARG_6="0"        #plots
ARG_7="1"         #k-est (default value used to be 6, make sure to specify for consistency) same thing that we have as k0 in TightCluster but kMin in the gene analysis app...
ARG_8="0.2"       #alpha
ARG_9="0.9"       #beta
ARG_10="7"        #top-can (q in the tight clustering class)
ARG_11="10"       #resamp.num (b) 
ARG_12="1"       #clusterRawData
ARG_13="1"       #clusterSmoothedData
ARG_14="1"       #clusterCoefficients
ARG_15="1"       #clusterLoose
ARG_16="1"       #clusterTight
ARG_17="0"      #clusterGap
ARG_18="0.7"     #sample ratio
ARG_19="3"       #levels

while getopts "R:k:r:m:spa:b:B:q:t:ahclgxL:" opt ; do
    case $opt in
	r) ARG_3="$OPTARG" ;;
	    #echo "row sum"
	    #case "$2" in
		#"") ARG_3="100" ; echo "empty row sum arg found" ; shift 2 ;;      #rowSum
		#*) ARG_3="$2" ; shift 2 ;;
	    #esac ;;
	m) ARG_4="$OPTARG" ;;
	    #case "$2" in
		#"") ARG_4="6" ; shift 2 ;;        #k-max
		#*) ARG_4="$2" ; shift 2 ;;
	    #esac ;;			  
	s) ARG_5="0" ;; #shift ;;		   
	p) ARG_6="1" ;; #shift ;;
	k) ARG_7="$OPTARG" ;;                      #k-est
	    #case "$2" in
		#"") ARG_7="6" ; shift 2 ;;
		#*) ARG_7="$2" ; shift 2 ;;
	    #esac ;;
	a) ARG_8="$OPTARG" ;;
	    #case "$2" in
		#"") ARG_8="0.2" ; shift 2 ;;      #alpha
		#*) ARG_8="$2" ; shift 2 ;;
	#esac ;;
	B) ARG_9="$OPTARG" ;;
       	    #case "$2" in
		#"") ARG_9="0.9" ; shift 2 ;;      #beta
		#*) ARG_9="$2" ; shift 2 ;;
	    #esac ;;
	q) ARG_10="$OPTARG" ;;
	    #case "$2" in
		#"") ARG_10="7" ; shift 2 ;;       #top-can (i.e. - q) 
		#*) ARG_10="$2" ; shift 2 ;;
	    #esac ;;
	b) ARG_11="$OPTARG" ;;
	    #case "$2" in
		#"") ARG_11="6" ; shift 2 ;;        #target
		#*) ARG_11="$2" ; shift 2 ;;
	    #esac ;;
	w) ARG_12="0" ;;
	t) ARG_13="0" ;;
	c) ARG_14="0" ;;
	l) ARG_15="0" ;;
	g) ARG_16="0" ;;
	x) ARG_17="1" ;;
	R) ARG_18="$OPTARG" ;;
	L) ARG_19="$OPTARG" ;;
	h) echo "$HELP"; exit 0 ;;
	--) shift ; break ;;
	*)  echo "$USAGE" ; exit 1 ;;     
   esac
done
echo "CLI options: $ARG_1 $ARG_2 $ARG3 $ARG_4 $ARG_5 $ARG_6 $ARG_7 $ARG_8 $ARG_9 $ARG_10 $ARG_11 $ARG_12 $ARG_13 $ARG_14 $ARG_15 $ARG_16 $ARG_17 $ARG_18 $ARG_19"
sbt "run-main apps.GeneAnalysis2 $ARG_1 $ARG_2 $ARG_3 $ARG_4 $ARG_5 $ARG_6 $ARG_7 $ARG_8 $ARG_9 $ARG_10 $ARG_11 $ARG_12 $ARG_13 $ARG_14 $ARG_15 $ARG_16 $ARG_17 $ARG_18 $ARG_19"

