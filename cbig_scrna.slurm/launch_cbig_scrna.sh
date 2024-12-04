#!/bin/bash

VERSION=0.0.1
DATE0=2024-12-04
echo -e "C-BIG scRNAseq pipeline version $VERSION"

# ===============================================
# default variables values
# ===============================================
unset SAMPLE OUTPUT_DIR PIPELINE_HOME QUEUE ACCOUNT 

#Assuming script is in root of PIPELINE_HOME
#set queue-specific values, depending on system check
# submit_cmd="bash"
# if [  -z "${submit_cmd}"  ]; then
#    submit_cmd="bash"
# fi
# # echo $submit_cmd

# if [[ $submit_cmd =~ sbatch ]]; then
#    QUEUE="sbatch"          # default job scheduler: qsub
# elif [[ $submit_cmd =~ bash ]]; then
#     QUEUE="bash"
# else
#     echo "The pipeline is testing under slurm system,"
#     echo "please choose sbatch for slurm."   
#     exit 42
# fi


PIPELINE_HOME0=`realpath ${BASH_SOURCE[0]}`
export PIPELINE_HOME=$(cd $(dirname $PIPELINE_HOME0) && pwd -P)

TIMESTAMP=`date +%FT%H.%M.%S`

# create function to handle error messages
# ===============================================
Usage() {
	echo
  echo "------------------- " 
	echo -e "Usage:\t$0 [arguments]"
	echo -e "\tmandatory arguments:\n" \
          "\t\t-d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)\n" \
          "\t\t--steps  =  Specify what steps, e.g., 1 to run step 1; 2 to run step 2\n" 
	echo -e "\toptional arguments:\n " \
          "\t\t-h  (--help)  = See helps regarding the pipeline arguments. \n" \
          "------------------- \n" \
          "For a comprehensive help, please contact Michael Fiorini, Farhan Lab (michael.fiorini@mail.mcgill.ca)."

echo 
}

# ===============================================
# PARSING ARGUMENTS
# ===============================================
if ! options=$(getopt --name pipeline --unquoted --options d:h --longoptions dir:,steps:,help -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    Usage
    exit 42
fi

# ===============================================
# LOAD & OVERRIDE EXTRA CONFIG FILE FOR PROJECT
# ===============================================
set -- $options

while [ $# -gt 0 ]

do
    case $1 in
    -x| --extra) 
      EXTRA_CONF="$2" ;
      if [ -f $EXTRA_CONF ]; then
        echo "* LOADING EXTRA CONFIG FILE $EXTRA_CONF";
        . $EXTRA_CONF
      else
        echo "ERROR: invalid EXTRA CONFIG file: $EXTRA_CONF";
        echo "Please check options and try again"; exit 42;
      fi
    esac
    shift
done

if [[ -n $ACCOUNT ]]; then ACCOUNT="-A $ACCOUNT"; fi
if [[ $MODULEUSE ]]; then module use $MODULEUSE ; fi
if [[ ${CELLRANGER} ]]; then module load ${CELLRANGER}/${CELLRANGER_VERSION} ; fi


# ===============================================
# LOAD ALL OTHER OPTIONS
# ===============================================
set -- $options

while [ $# -gt 0 ]
do
    case $1 in
    -h| --help) Usage; exit 0;;
    -d| --dir) OUTPUT_DIR="$2" ; shift ;;
    -v| --verbose) VERBOSE=1 ;; 
    --steps) MODE="$2"; shift ;;   
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 42;;
    (*) break;;
    esac
    shift
done

JOB_MODE=slurm

source $PIPELINE_HOME/tools/utils.sh
FOUND_ERROR=0

# ===============================================
# CHECKING VARIABLES
# ===============================================
#check to ensure all mandatory arguments have been entered

if [ -z $OUTPUT_DIR ]; then echo "ERROR: missing mandatory option: -d (--dir) must be specified"; FOUND_ERROR=1; fi

if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi

if [[ -z $SINFO ]]; then  SINFO="F"; fi

# if [ $MODE == 'ALL' ]; then  MODE00=`echo {2..10}`; MODE0=`eval echo $MODE00`; fi 
if [[ "$MODE" == *"-"* ]]; then
  MODE00=`echo $MODE | sed  "s/-/../g"  | awk '{print "{"$0"}"}'`
  MODE0=`eval echo $MODE00`
else 
  MODE0=$MODE
fi

JOB_MODE=slurm

# STEP 0: RUN setting 
# ===============================================
#


if [[ ${MODE0[@]} == 0 ]]; then 
  JOB_OUTPUT_DIR=$OUTPUT_DIR/job_info
  if [[ -s $OUTPUT_DIR/job_info/configs ]]; then
    # echo "config file already exists in $OUTPUT_DIR. Overwrite? y|n"
    # read answer
    read -p "config file already exists in $OUTPUT_DIR. Overwrite? y|n: "$'\n'  answer
    # answer=${answer,,}; answer=${answer:0:1}
    # echo $answer
    if [[ $answer =~ y ]]; then 
      echo "NOTE: the config file and parameters are Overwritten."
        if [ ! -d $JOB_OUTPUT_DIR ]; then   
        mkdir -p $OUTPUT_DIR/job_info
        mkdir -p $OUTPUT_DIR/job_info/parameters
        mkdir -p $OUTPUT_DIR/job_info/logs
        mkdir -p $OUTPUT_DIR/job_info/configs
        mkdir -p $OUTPUT_DIR/job_info/.tmp        
        fi
        if [ -d $OUTPUT_DIR/job_info ]; then
            rm -r $OUTPUT_DIR/job_info
            mkdir -p $OUTPUT_DIR/job_info
            mkdir -p $OUTPUT_DIR/job_info/parameters
            mkdir -p $OUTPUT_DIR/job_info/logs
            mkdir -p $OUTPUT_DIR/job_info/configs
            mkdir -p $OUTPUT_DIR/job_info/.tmp
        fi
        touch $JOB_OUTPUT_DIR/summary_report.txt
        cp $PIPELINE_HOME/scrna/configs/cbig_scrna_config.ini $OUTPUT_DIR/job_info/configs/
        cp $PIPELINE_HOME/scrna/pars/* $OUTPUT_DIR/job_info/parameters/
    else
     echo "NOTE: the pipeline is using the existing config file and parameters."
    fi
  else 
    echo "The configuration files do not exist, the pipeline will create them during execution."
        if [ ! -d $JOB_OUTPUT_DIR ]; then 
        mkdir -p $OUTPUT_DIR/job_info
        mkdir -p $OUTPUT_DIR/job_info/parameters
        mkdir -p $OUTPUT_DIR/job_info/logs
        mkdir -p $OUTPUT_DIR/job_info/configs
        mkdir -p $OUTPUT_DIR/job_info/.tmp                
        fi
        touch $JOB_OUTPUT_DIR/summary_report.txt
        cp $PIPELINE_HOME/scrna/configs/cbig_scrna_config.ini $OUTPUT_DIR/job_info/configs/
        cp $PIPELINE_HOME/scrna/pars/* $OUTPUT_DIR/job_info/parameters/
  fi
  if [[   ${SINFO}  =~  T ]]; then
    if  [ ! -d $OUTPUT_DIR/samples_info ]; then   echo 'ERROR: The pipeline can not find samples_info directory,  '; exit 42 ; fi 
    if  [ -d ${OUTPUT_DIR}/step1 ]; then 
      rm ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
      rm ${OUTPUT_DIR}/job_info/.tmp/sample.list
      search_dir=${OUTPUT_DIR}/step1
      for entry in "$search_dir"/*
        do
          echo "$entry" >> ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
          echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_info/.tmp/sample.list
        done
   fi
  fi 
        EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/summary_report.txt
        echo -e  "------------------------------------------------------" >> $EXPECTED_DONE_FILES  
        echo -e  "----------------- Pipeline is set up -----------------" $VERSION >> $EXPECTED_DONE_FILES
        echo "----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
        echo "The Output is under ${OUTPUT_DIR}/job_info/" >> $EXPECTED_DONE_FILES
#  exit 0
fi 

declare -A THREADS_ARRAY
declare -A  WALLTIME_ARRAY
declare -A  MEM_ARRAY
source $OUTPUT_DIR/job_info/configs/cbig_scrna_config.ini

export JOB_OUTPUT_DIR=$OUTPUT_DIR/job_info
export EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/summary_report.txt
chmod 775 $EXPECTED_DONE_FILES
if [[ $MODULEUSE ]]; then export MODULEUSE=$MODULEUSE ; fi
export CELLRANGER=${CELLRANGER}
export ACCOUNT=$ACCOUNT
export OUTPUT_DIR=$OUTPUT_DIR

echo -e "NOTE: the pipeline is running"

TEMPCONFIG=$OUTPUT_DIR/job_info/.tmp/temp_config.ini
touch $TEMPCONFIG
if [[ -f "TEMPCONFIG" ]]; then 
  rm $TEMPCONFIG
fi

if [  -z "${JOB_MODE}"  ]; then
   JOB_MODE=local
fi

if [[ $JOB_MODE =~ slurm ]]; then
   QUEUE="sbatch"
elif [[ $JOB_MODE =~ local ]]; then
    QUEUE="bash"
else
    echo "The pipeline is tested under slurm system, and linux"
    echo "please choose choose slurm or local"   
    exit 42
fi

echo " # IT IS A temp FILE. DO NOT EDIT THIS FILE DIRECTLY."  > $TEMPCONFIG
echo SINFO=$SINFO  >> $TEMPCONFIG
echo OUTPUT_DIR=$OUTPUT_DIR >> $TEMPCONFIG
echo JOB_OUTPUT_DIR=$JOB_OUTPUT_DIR  >> $TEMPCONFIG
echo EXPECTED_DONE_FILES=$EXPECTED_DONE_FILES   >> $TEMPCONFIG
echo MODE=$MODE >> $TEMPCONFIG
echo QUEUE=$QUEUE >> $TEMPCONFIG
echo VERSION=$VERSION >> $TEMPCONFIG
if [[ $MODULEUSE ]]; then 
    echo MODULEUSE=$MODULEUSE >> $TEMPCONFIG
fi
if [[ $CELLRANGER ]]; then 
    echo CELLRANGER=$CELLRANGER >> $TEMPCONFIG
fi
if [[ ${CELLRANGER_VERSION} ]]; then 
    echo CELLRANGER_VERSION=${CELLRANGER_VERSION} >> $TEMPCONFIG
fi

bash ${PIPELINE_HOME}/launch/launch_cbig_scrna_slurm.sh

echo -e " \n"
exit 0



