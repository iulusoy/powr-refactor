#!/bin/bash

# source ~/.powrconfig || exit

####################################################################
#
# Parameters: 
# $1 : Job-file
# $2 : dbx           ==> debugger running
#      to-nam        ==> Job to be started on machine <nam>
#
# This a greatly simplified version, compared to lars_submit.com, 
#   and does not support migration to machines outside our cluster.
# 
####################################################################

#Define jobfile
job=$POWR_WORK/wrjobs/$1
# Define path of job.log and job.cpr files
pathout=$POWR_WORK/output/
# Define logfile of this submit procedure
logfile=$pathout/submit.log
# Define proc.dir
proc=${procdir:-$POWR_WORK/proc.dir}
# Define list-of-machines
lom=${lom:-$POWR_WORK/wrjobs/list-of-machines}

echo ------------------------------------------------ >> $logfile
echo Start of submit.com on `date` >> $logfile
echo ------------------------------------------------ >> $logfile

if [ ! -x $job ]
 then
  echo File $job does not exist or has not execute permission
  echo File $job does not exist or has not execute permission >> $logfile
  echo submit ABORTED! 
  echo submit ABORTED! >> $logfile
  exit
fi

############################################
# $wrem = target machine; default: localhost
############################################
HOST=`hostname | cut -d '.' -f 1`
wrem=$HOST

if [ $# -ge 2 ]; then
# ---- Debugger mode ----------------------------
  if [ $2 = 'dbx' ]; then
    echo 'submit.com: $1 with dbx gestartet'
    $job dbx
    echo $1 with dbx gestartet >> $logfile
    exit

# ---- to-machine ------------------------------
  elif [ `echo $2 | cut -c1-3` = 'to-' ]; then
###################################################
# $wrem = target machine specified by to-<nam> 
# localhost is used instead, if specified machine:
#         - is not in list-of-machines
#         - does not reply on ping 
################################################### 
    wrem=`echo $2 | cut -c4-`
    if grep -q "^$wrem$" $lom; then
        nping=`ping -c 5 $wrem | grep loss | cut -f 4 -d' '`
        if [ $nping -lt 1 ]; then
          echo submit.com : Ping failed on machine $wrem
          echo Job instead started on localhost $HOST 
          echo submit.com : Ping failed on machine $wrem >> $logfile
          echo Job instead started on localhost $HOST >> $logfile
          wrem=$HOST
        fi
    else
      echo submit.com : Machine $wrem unknown
      echo              Job on actual machine started
      echo submit.com : Machine $wrem unknown >> $logfile
      echo Job started on local host $HOST >> $logfile
      wrem=$HOST
    fi

#------- Error exit 
  else
    echo ERROR: unknown option $2
    echo submit ABORTED
    echo ERROR: unknown option $2 >> $logfile
    echo submit ABORTED           >> $logfile
    exit
  fi
fi

############################################################################
# if submit already runs at the target host, the job is started;
# ... else, a new submit is started at the target host in a remote-shell (rsh)
############################################################################
if [ $wrem = $HOST ]; then
  $job  >> $pathout$1.log 2>> $pathout$1.cpr &
  echo $1 started on $HOST  
  echo
  echo $1 started on $HOST >> $logfile
  echo                     >> $logfile
  echo ------------------------------------------------------ >> $pathout$1.log
  echo ------------------------------------------------------ >> $pathout$1.cpr
  echo submit.com : Job started on `date` >> $pathout$1.log
  echo submit.com : Job started on `date` >> $pathout$1.cpr
  echo ------------------------------------------------------ >> $pathout$1.log
  echo ------------------------------------------------------ >> $pathout$1.cpr
else
  if [ $OSTYPE = 'osf1' ] ; then
   rsh $wrem $proc/submit.com $1 &
   echo $1 submitted "(rsh)" on `date` from machine $HOST to $wrem 
   echo $1 submitted "(rsh)" on `date` from machine $HOST to $wrem >> $logfile
  else
   ssh $wrem $proc/submit.com $1 &
   echo $1 submitted "(ssh)" on `date` from machine $HOST to $wrem 
   echo $1 submitted "(ssh)" on `date` from machine $HOST to $wrem >> $logfile
  fi
fi

# wait with the prompt
sleep 4
