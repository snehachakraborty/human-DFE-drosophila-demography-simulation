#!/bin/bash -l

# run parameters (h_rt: time, h_data: memory, highp: high priority which is good for long running jobs)
#$ -l h_rt=336:00:00,h_data=50G,highp

# output directory for stout (-o) and sterr (-e) logs
#$ -o /u/home/s/sneha_c/project-klohmuel/meixi_simulations/rescaled_1Mb_no_recom_map/logs/joblog.$JOB_ID.$TASK_ID
#$ -e /u/home/s/sneha_c/project-klohmuel/meixi_simulations/rescaled_1Mb_no_recom_map/logs/joblog.$JOB_ID.$TASK_ID

# name of the job that appears in the queue (see your job queue using qstat -u <username>)
#$ -N array-rescaled-dros-hum-sim-1-108

# account to which the job is charged
#$ -A klohmuel

# number of tasks in job array for parallelization
# this means the task IDs will be 1, 2, 3, ..., 108 and the :1 means that each task will have a step size of 1
#$ -t 1-108:1

# Notify when
#$ -m bea

slim_script="/u/home/s/sneha_c/project-klohmuel/meixi_simulations/rescaled_1Mb_no_recom_map/scripts/rescaled_drosophila_dem_human_dfe_ML.slim"



slim -d run=${SGE_TASK_ID} -d model=1 -d seed=19405 -d "out_path='/u/home/s/sneha_c/project-klohmuel/meixi_simulations/rescaled_1Mb_no_recom_map/out_files'" -d "sim_annots='/u/home/s/sneha_c/project-klohmuel/meixi_simulations/rescaled_1Mb_no_recom_map/sim_files/sim_${SGE_TASK_ID}_annots.csv'" -d "log_path='/u/home/s/sneha_c/project-klohmuel/meixi_simulations/rescaled_1Mb_no_recom_map/logs'" ${slim_script}


