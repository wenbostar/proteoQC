#!/bin/bash

outdir='./qc'
mode='identification'
miss=0
enzyme=1
varmod=2
fixmod=1
tol=10
itol=0.6
cpu=2

subcommand=$1; shift
if [ "$subcommand" = "report" ]; then
  while true; do
    if [ "$1" = "-h" ]; then
        shift;
        echo -e "\nUsage: proteoQC [module] -s [spectralist] -f [fasta] -[options] [values]\n";
        echo -e "Options:\n"
        echo -e "-s    A file contains the experiment design or a single mgf file"
        echo -e "-f    Database file, must contain decoy sequences"
        echo -e "-o    Output directory | default: ./result"
        echo -e "-d    Identification or quantification | default: identification"
        echo -e "-m    Max miss clevages | default: 0"
        echo -e "-e    Enzyme | default: 1"
        echo -e "-v    Variable modifications are those which may or may not be present | default: 2"
        echo -e "-x    Fixed modifications are applied universally, to every instance of the specified residue(s) or terminus | default: 1"
        echo -e "-t    The error window on experimental peptide mass values | default: 10"
        echo -e "-i    Error window for MS/MS fragment ion mass values | default: 0.6"
        echo -e "-c    Max number of cpu used | default: 2"
    else
      case $1 in
        -s  )
          shift; spectralist=$1; shift;;
        -f  )
          shift; fasta=$1; shift;;
        -o  )
          shift; outdir=$1; shift ;;
        -d  )
          shift; mode=$1; shift ;;
        -m  )
          shift; miss=$1; shift ;;
        -e  )
          shift; enzyme=$1; shift ;;
        -v  )
          shift; varmod=$1; shift ;;
        -x  )
          shift; fixmod=$1; shift ;;
        -t  )
          shift; tol=$1; shift ;;
        -i  )
          shift; itol=$1; shift ;;
        -c  )
          shift; cpu=$1; shift ;;
        -*  )
          echo "$0: Unrecognized option $1" >&2
          exit 2;;
        *) break ;;
      esac

      echo "$spectralist $fasta $miss"

      #INSTALL
      # R -e "if (!require('proteoQC')) {source('https://bioconductor.org/biocLite.R');\
      # biocLite('proteoQC');}"

      #RUN
      R -e "library(proteoQC);\
            design <- system.file('extdata/$spectralist-design.txt', package='proteoQC');\
            fas <- unzip('$fasta');\
            qcres <- msQCpipe(spectralist=design,fasta=fas,outdir ='$outdir',\
            miss=$miss,enzyme=$enzyme,varmod=$varmod,fixmod=$fixmod,\
            tol=$tol,itol=$itol,cpu=$cpu,mode='$mode');\
            zpqc <- system.file('extdata/qc.zip', package='proteoQC');\
            unzip(zpqc);\
            qcres <- loadmsQCres('./qc');\
            html <- reportHTML(qcres)"
    fi
  done
else
  echo "$subcommand is not a valid module"
fi
