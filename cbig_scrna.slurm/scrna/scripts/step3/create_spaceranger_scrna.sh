#!/bin/bash

source $3
OUTFILE=$1/$2
if [ -f $OUTFILE ]; then
    rm $OUTFILE
fi

cat <<EOF > $OUTFILE
#!/bin/bash
unset RUN_NAME LIBRARY EXPECT_CELLS SLURM_TEMPLATE failure REF_DIR

Usage() {
        echo
        echo -e "Usage:\t\$0 [arguments]"
        echo -e "\tmandatory arguments:\n" \\
          "\t\t-r  (--run_name)         run name, serves as prefix for output files\n"
        echo -e "\toptional arguments:\n" \\
          "\t\t-l  (--library_csv)      /path/to/library.csv file       (default: ./library.csv)\n" \\
          "\t\t-t  (--slurm_template)   /path/to/slurm.template file    (default: ./slurm.template)\n" \\
          "\t\t-h  (--help)             run this help message and exit\n"
        echo
}

if ! options=\$(getopt --name \$(basename \$0) --alternative --unquoted --options hr:l:f:t: --longoptions run_name:,library_csv:,slurm_template:,help -- "\$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    exit 42
fi

set -- \$options
while [ \$# -gt 0 ]
do
    case \$1 in
    -h| --help) Usage; exit 0;;
    -r| --run_name) RUN_NAME="\$2"; shift ;;
    -l| --library_csv) LIBRARY="\$2"; shift ;;
    -t| --slurm_template) SLURM_TEMPLATE="\$2"; shift ;;
    (--) shift; break;;
    (-*) echo "\$0: error - unrecognized option \$1" 1>&2; failure=1;;
    (*) break;;
    esac
    shift
done

[[ -z \${RUN_NAME} ]] && echo "error: run name must be specified with -r or --run_name" && failure=1
[[ \$failure -eq 1 ]] && echo "ERRORS FOUND. Exiting" && exit 42

source \$OUTPUT_DIR/job_info/.tmp/parameters.txt
source \$OUTPUT_DIR/job_info/.tmp/temp_config.ini

LIBRARY=\${LIBRARY:-./library.csv}
SLURM_TEMPLATE=\${SLURM_TEMPLATE:-./slurm.template}

#if [[ \$MODULEUSE ]]; then module use \$MODULEUSE ; fi
#module load \${SPACERANGER}/\${SPACERANGER_VERSION}
export PATH=/home/fiorini9/projects/def-sfarhan/fiorini9/software/spaceranger/spaceranger-3.1.3:$PATH

pwd_dir=\$(pwd)
TEMPLOG=\$OUTPUT_DIR/job_info/logs/step_3_spaceranger\$(basename \${pwd_dir}).log
echo "SPACE RANGER is currently running on \$(basename \${pwd_dir}). Please leave it undisturbed until it finishes."

EOF

cat <<EOF >> $OUTFILE
spaceranger count \\
    --id=\${RUN_NAME} \\
    --fastqs=\${par_fastqs_spaceranger} \\
    --jobmode=\${SLURM_TEMPLATE} \\
    --create-bam=true \\
EOF

if [[ -n "${par_ref_dir_grch_spaceranger}" ]]; then
cat <<EOF >> $OUTFILE
    --transcriptome=\${par_ref_dir_grch_spaceranger} \\
EOF
fi

if [[ -n "${par_image_spaceranger}" ]]; then
cat <<EOF >> $OUTFILE
    --image=\${par_image_spaceranger} \\
EOF
fi

if [[ -n "${par_slide_spaceranger}" ]]; then
cat <<EOF >> $OUTFILE
    --slide=\${par_slide_spaceranger} \\
EOF
fi

if [[ -n "${par_slidefile_spaceranger}" ]]; then
cat <<EOF >> $OUTFILE
    --slidefile=\${par_slidefile_spaceranger} \\
EOF
fi

if [[ -n "${par_area_spaceranger}" ]]; then
cat <<EOF >> $OUTFILE
    --area=\${par_area_spaceranger} \\
EOF
fi

if [[ -n "${par_localcores_spaceranger}" ]]; then
cat <<EOF >> $OUTFILE
    --localcores=\${par_localcores_spaceranger} \\
EOF
fi

if [[ -n "${par_localmem_spaceranger}" ]]; then
cat <<EOF >> $OUTFILE
    --localmem=\${par_localmem_spaceranger} \\
EOF
fi


########## END 
cat <<EOF >> $OUTFILE
    2>&1|tee -a \${RUN_NAME}.\$(date +%Y%m%d_%H%M).log >  \${TEMPLOG}
EOF

cat <<EOF >> $OUTFILE
echo -e "The computation on  \$(basename \${pwd_dir}) is done. "
awk '/Pipestance/'  \${TEMPLOG}
EOF




###################################################
###################################################
###################################################
################################################### ++++++ OLD CODE BEYOND THIS POINT. 
## insert this at line 73 if needed. 
#if [[ -n "${par_create_bam_spaceranger}" ]]; then
#if  [[  ${par_create_bam_spaceranger}  =~  yes  ]]  ; then
#cat <<EOF >> $OUTFILE
#    --create-bam=true \\
#EOF
#fi
#fi


#if [[ -n "${par_fastqs_spaceranger}" ]]; then
#cat <<EOF >> $OUTFILE
#    --fastqs=\${par_fastqs_spaceranger} \\
#EOF
#fi


#EOF

#cat <<EOF >> $OUTFILE
#spaceranger count \\
#    --id=\${RUN_NAME} \\
#    --libraries=\${LIBRARY} \\
#    --jobmode=\${SLURM_TEMPLATE} \\
#    --create-bam=true \\
#EOF

