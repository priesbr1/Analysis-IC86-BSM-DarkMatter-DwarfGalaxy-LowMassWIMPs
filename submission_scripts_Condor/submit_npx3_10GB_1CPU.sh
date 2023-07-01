#!/bin/sh

# A simple sh script for submitting one job to the npx3 cluster.
# Simply log into npx3.icecube.wisc.edu, run any command as you normally 
#  would but prepend "./submit-npx3.sh" and the command will be run on the
#  npx3 cluster.
#
# Eg #1:
# ./submit-npx3.sh root -l -b -q macro.C
# (NB: You must run in batch mode "-b" and quit root when macro is done "-q")
#
# Eg #2:
# ./submit-npx3.sh ./env-shell.sh python icetray_job.py
# (NB: You must execute env shell when submitting icetray jobs)
#
# This script will create directories to store your execution script, log files,
#  errors, and std output, so you need write permission in the local directory.

# This script creates a script to be executed and another script to submit it.
# The execution script must be available *at time of job execution!*, which may
#  not be until much later and so it's stored in a directory 'npx3-execs'.
# You may occasionally want to 'rm -rf npx3-*' directories if they get big.
# The submission script "2sub.sub" can be discarded immediately.

# This method of passing your job into a bash script as arguments may fail
#  if you have quotes or other special characters

#Quickest ever Condor tutorial:
#"condor_q" gives list of jobs running
#"condor_q $USER" gives list of your jobs
#"condor_rm $USER" removes your jobs 

  # Creating output directories
  mkdir -p npx3-execs npx3-logs npx3-out npx3-error

  # Creating execution script, do not delete until job has started!
  echo "#!/bin/bash" > npx3-execs/npx3-$$.sh
  echo "date" >> npx3-execs/npx3-$$.sh
  echo "hostname" >> npx3-execs/npx3-$$.sh
  echo "cd `pwd`" >> npx3-execs/npx3-$$.sh
  #Set up new tools
  echo "eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`" >> npx3-execs/npx3-$$.sh
  echo "$@" >> npx3-execs/npx3-$$.sh
  echo "date" >> npx3-execs/npx3-$$.sh

  chmod +x npx3-execs/npx3-$$.sh

  # Creating condor submission script (ClassAd)
  echo "Universe  = vanilla" > 2sub.sub
  echo "Executable = npx3-execs/npx3-$$.sh" >> 2sub.sub
  echo "Log = npx3-logs/npx3-$$.log" >> 2sub.sub
  echo "Output = npx3-out/npx3-$$.out" >> 2sub.sub
  echo "Error = npx3-error/npx3-$$.error" >> 2sub.sub
  echo "request_cpus = 1" >> 2sub.sub
  #echo "request_gpus = 1" >> 2sub.sub
  echo "+NATIVE_OS = true" >> 2sub.sub
  echo "Requirements = (OpSysMajorVer =?= 7)" >> 2sub.sub
  echo "request_disk = 100MB" >> 2sub.sub
  echo "Request_memory = 5 GB" >> 2sub.sub
  #echo "Notification = Complete" >> 2sub.sub 
  #echo "+AccountingGroup="1_week.mnisa"" >> 2sub.sub
  echo "Queue" >> 2sub.sub
 # echo "notify_user = mnisa@icecube.wisc.edu" >> 2sub.sub
  condor_submit 2sub.sub
