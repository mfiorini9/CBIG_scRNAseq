#!/bin/bash

cbig_scrna_config.ini

declare -A THREADS_ARRAY
declare -A  WALLTIME_ARRAY
declare -A  MEM_ARRAY
source $OUTPUT_DIR/job_info/configs/cbig_scrna_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini
source $PIPELINE_HOME/tools/utils.sh

# if [ $MODE == 'ALL' ]; then  MODE00=`echo {2..8}`; MODE0=`eval echo $MODE00`; fi 
if [[ "$MODE" == *"-"* ]]; then
  MODE00=`echo $MODE | sed  "s/-/../g"  | awk '{print "{"$0"}"}'`
  MODE0=`eval echo $MODE00`
else 
  MODE0=$MODE
fi


export QUEUE=$QUEUE  
# ===============================================
# STEP 1: 
# ===============================================
#
STEP=step_1

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  1  ]]; then
  export ACCOUNT=$ACCOUNT
  if [ ! -d $OUTPUT_DIR/step1 ]; then 
    mkdir -p $OUTPUT_DIR/step1 
  fi
  step1_par_auto=`grep 'par_automated_library_prep=' ${OUTPUT_DIR}/job_info/parameters/parameters.txt`
  step1_par_auto1=`echo ${step1_par_auto//[[:blank:]]/} | tr 'A-Z' 'a-z'`
  if [[ "${step1_par_auto1}" == "par_automated_library_prep=\"yes\"" ]]; then
        module load r/$R_VERSION 
        Rscript ${PIPELINE_HOME}/scrna/scripts/step1/scrna_step1_auto.R $OUTPUT_DIR  $R_LIB_PATH 
        echo "Note: You are generating samples_info automatically"
  fi

  cp -r  $OUTPUT_DIR/samples_info/* $OUTPUT_DIR/step1 
  if  [ -f $JOB_OUTPUT_DIR/.tmp/sample.list ]; then rm $JOB_OUTPUT_DIR/.tmp/sample.list; touch $JOB_OUTPUT_DIR/.tmp/sample.list ; else touch $JOB_OUTPUT_DIR/.tmp/sample.list; fi 
  if  [ -f $JOB_OUTPUT_DIR/.tmp/sample_dir.list ]; then rm $JOB_OUTPUT_DIR/.tmp/sample_dir.list; touch $JOB_OUTPUT_DIR/.tmp/sample_dir.list ; else touch $JOB_OUTPUT_DIR/.tmp/sample_dir.list; fi  
  if [ -f $OUTPUT_DIR/job_info/.tmp/parameters.txt ]; then
    rm $OUTPUT_DIR/job_info/.tmp/parameters.txt
  fi
  
  grep "par_ref_dir_grch=" $OUTPUT_DIR/job_info/parameters/parameters.txt | sed 's/\"//g' | sed "s/'//g" | sed "s/[[:blank:]]//g" > $OUTPUT_DIR/job_info/.tmp/parameters.txt
  arr=(par_include_introns par_mempercode par_r1_length par_r2_length par_no_target_umi_filter par_expect_cells par_force_cells par_no_bam par_no_libraries) 
  for item in ${arr[@]}; do
      grep ${item}=  $OUTPUT_DIR/job_info/parameters/parameters.txt | sed 's/\"//g' | sed "s/'//g" | sed "s/[[:blank:]]//g" | sed 's/[A-Z]/\L&/g' >> $OUTPUT_DIR/job_info/.tmp/parameters.txt
  done

  search_dir=${OUTPUT_DIR}/step1
  for entry in "$search_dir"/*
  do
    echo "$entry" >> ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
    echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_info/.tmp/sample.list
  done
  
  source $OUTPUT_DIR/job_info/.tmp/parameters.txt
  bash ${PIPELINE_HOME}/scrna/scripts/step1/create_cellranger_scrna.sh $OUTPUT_DIR/job_info/.tmp/ scrna_cellranger.slurm.sh $OUTPUT_DIR/job_info/.tmp/parameters.txt

  while read item
  do
    cp ${PIPELINE_HOME}/scrna/scripts/step1/slurm.template $item
    sed -i $item/slurm.template -e 's/account0/'$ACCOUNT'/'
    cp $OUTPUT_DIR/job_info/.tmp/scrna_cellranger.slurm.sh $item
  done < ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list

  while read item
  do
      cd ${item}; bash  ${item}/scrna_cellranger.slurm.sh -r ouput_folder &
  done < ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list

  wait
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo -e  "-------------------------------------------" $VERSION >> $EXPECTED_DONE_FILES  
  echo -e  "--------Job submitted using pipeline-------" $VERSION >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "step 1  submited " >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step1/" >> $EXPECTED_DONE_FILES
fi 

# ===============================================
# ===============================================

if  [ ! -f $JOB_OUTPUT_DIR/.tmp/sample.list ]; then
    search_dir=${OUTPUT_DIR}/step1
    for entry in "$search_dir"/*
    do
      echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_info/.tmp/sample.list
    done
fi


if  [ ! -f $JOB_OUTPUT_DIR/.tmp/sample_dir.list ]; then
    search_dir=${OUTPUT_DIR}/step1
    for entry in "$search_dir"/*
    do
      echo "$entry" >> ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
    done
fi 

SAMPLE_SIZE=`wc -l < ${OUTPUT_DIR}/job_info/.tmp/sample.list`


# ===============================================
# STEP 2:
# ===============================================
#
STEP=step_2

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  2 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*4)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*4))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*80)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  if [  -d $OUTPUT_DIR/step2/objs2 ]; then 
    rm -rf  $OUTPUT_DIR/step2/objs2 ; mkdir -p $OUTPUT_DIR/step2/objs2 &
  else
    mkdir -p $OUTPUT_DIR/step2/objs2 &
  fi 
  if [  -d $OUTPUT_DIR/step2/figs2 ]; then  
    rm -rf  $OUTPUT_DIR/step2/figs2 ; mkdir -p $OUTPUT_DIR/step2/figs2  &
    else 
    mkdir -p $OUTPUT_DIR/step2/figs2 
  fi
  if [  -d $OUTPUT_DIR/step2/info2 ]; then    
    rm -rf  $OUTPUT_DIR/step2/info2 ; mkdir -p $OUTPUT_DIR/step2/info2 &  
    else 
    mkdir -p $OUTPUT_DIR/step2/info2   
  fi
fi 

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  2 ]]  &&  [[  ${MODE0[@]} =~ 1 ]] ; then
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 2 submitted following step 1"  >> $EXPECTED_DONE_FILES
  step_2="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_CELLRANGER \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step2/pipeline_step2.sh"
  step_2=$($step_2 | grep -oP "\d+")
  echo "STEP 2: $step_2 is submitted"
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  echo "[Q] STEP 2         : $step_2 " >> $EXPECTED_DONE_FILES
  DEPEND_step_2="--dependency=afterok:$step_2"; echo $DEPEND_step_2
  echo_general="STEP 2: Job Number:$step_2"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step2_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The result is under ${OUTPUT_DIR}/step2" >> $EXPECTED_DONE_FILES
elif [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  2  ]]  &&  [[  ${MODE0[@]} != 1 ]]; then
  # echo "just step 2 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 2 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_2="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step2/pipeline_step2.sh"
  step_2=$($step_2 | grep -oP "\d+")
  echo "[Q] STEP 2         : $step_2 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  DEPEND_step_2="--dependency=afterok:$step_2"
  echo_general="STEP 2: Job Number:$step_2"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step2_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/step2" >> $EXPECTED_DONE_FILES  
fi 

exit 0
