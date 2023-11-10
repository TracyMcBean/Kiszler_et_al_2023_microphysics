#!/bin/bash
set +e
# Script to run the wrapper and if wanted store the percentage of each 
# process tendency as well as create plots of the tendencies and percentage etc.


MAKE_PLOTS=true
STORE_PERCENTAGE=true
# Conda environment for plotting (can be exchanged with python environment)
my_conda_env=iconplotenv

BASE_DIR=$(pwd)
# Location of compiled microphysics wrapper
WRAPPER_DIR="/icon-2mom-mcrph-wrapper/run/"
# Where to find the simulation data needed as input
SIMULATION_DIR=""
# Where to store the output of the wrapper
DATA_DIR="/scratch/wrapper_output/"

DATE="20210501"
#while [ ${DATE:0:4} -le $(date -d "$DATE +1 year" +%Y) ]
# only loop over next 10 days
#while [ ${DATE:7:8} -le $(date -d "$DATE +2 day" +%d) ]
#for DATE in ${DATE_ARRAY[@]};
# Run over all days in 2021 and 2022 and 2023
while [ ${DATE:0:6} -le 202309 ]
do
    # reformat date from YYYY-MM-DD to YYYYMMDD
    #DATE=${DATE:0:4}${DATE:4:2}${DATE:8:2}  
    echo $DATE
    
    INFILE=${SIMULATION_DIR}/${DATE}_r600m_f2km/METEOGRAM_patch001_${DATE}_awipev.nc
    # hydrometeor masses from wrapper
    MASSFILE=${DATA_DIR}/wrapper_mass/wrapper_mass_${DATE}.nc
    # Tendencies from wrapper
    TENDFILE=${DATA_DIR}/wrapper_tend/wrapper_tend_${DATE}.nc

    # check if infile exists
    if [[ ! -f ${INFILE} ]]; then
        echo "File ${INFILE} does not exist."
    else # run the wrapper
        echo "File ${INFILE} exists"

        cd ${WRAPPER_DIR}
# this cannot be indented otherwise the EOF_NML isn't recognized
cat > namelist.nml <<-EOF_NML
    &settings_nml
    infile    = "${INFILE}"
    !infile    = "../testdata/METEOGRAM_patch001_20221004_awipev_3sec_3h.nc"
    outfile   = "${MASSFILE}"
    tendfile  = "${TENDFILE}"
    msg_level = 2
    inwp_gscp = 4
    dtime     = 3  ! Time step of microphysics is same as dtime in simulation
    /
    &mcrph_nml     ! This namelist is for changes applied to the microphysics
    do_cold_mcrph    = .true. ! If .false. only warm processes are run.
    do_iceriming     = .true.
    do_icenuc_homhet = .true.
    do_iceriming     = .true.
    /
    &invars_nml  
    input_3D  = .false.
    max_ncells = 50000
    invar_names = "P", "T", "W", "RHO", "QV", "QC", "QNC", "QR", "QNR", "QI", "QNI", "QS", "QNS", "QG", "QNG", "QH", "QNH", "NIACT", "height_2", "height"
    !invar_names = "w", "qv", "qc", "qi", "qr", "qs", "qg", "qh", "qnc", "qni", "qnr" ,"qns", "qng", "qnh", "ninact", "pres", "temp", "z_ifc", "z_mc", "rho" 
    /
EOF_NML

        if [[ -f ${MASSFILE} ]]; then
            echo "${MASSFILE} already exists."
            result=0
        else
            # Run wrapper
            ${WRAPPER_DIR}/run_2mom_wrapper  
            result=$?
            echo "result: $result"
        fi

        # check if wrapper failed
        if (( $result == 2 || $result == 1 )); then
        #if [ ${command_failed:-0} -eq 1 ]; then
            echo "Wrapper failed"
        else
            echo "Wrapper succeeded"
            #####------- Plotting --------#######
            cd ${BASE_DIR}
            # Load conda
            CONDA_PREFIX=$(conda info --base)
            source ${CONDA_PREFIX}/etc/profile.d/conda.sh
            conda activate ${my_conda_env}
            
            # This is currently not working for some reason the tendency name selection doesn't work anymore
            # if $STORE_PERCENTAGE; then
                # python store_tend_perc.py ${DATE} ${TENDFILE} ${MASSFILE} ${DATA_DIR}/wrapper_tend_perc/
            # fi
    # 
             if $MAKE_PLOTS; then
                 python ../plotting/plot_tendency.py ${DATE} ${DATE:0:4}${DATE:4:2} ${TENDFILE} ${MASSFILE}
             fi
    
            conda deactivate
        fi
    
    fi # finnish if file exists
    # next date
    DATE=$(date -d "$DATE +1 day" +%Y%m%d)

done
