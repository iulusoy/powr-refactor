#!/bin/bash

HOST=$(hostname -s)

if [ $# -eq 1 ] ; then
 if [ _$1 == _all ] ; then
     if [ _$POWR_INSTTYPE == _local ] ; then
         source ~/.powrconfig || exit
	 rhosts=$(cat ${POWR_WORK}/wrjobs/list-of-machines | grep "^[a-z]")
     elif [ _$POWR_INSTTYPE == _niscluster ] ; then
         source ~/.powrconfig || exit
         rhosts=$(cat ${POWR_HOSTS_FILE_NIS} | grep "^[a-z]")
     else  # Potsdam
	 rhosts=$(cat ~htodt/powr/wrjobs/list-of-machines | grep "^[a-z]")
     fi
 elif [ _${1:0:1} == _-  ] ; then
  if [ _$1 == _-help ] ; then
   echo "Valid arguments are hostnames, 'all', or none (localhost)"
  else
   echo "Unknown option: " $1
   echo "Valid option: -help"
   echo "Valid arguments are hostnames, 'all', or none (localhost)"
  fi
 else
  rhosts=$1
 fi
else
 rhosts=$HOST
fi

let nprogramtotal=0
let ncoretotal=0
for rhost in `echo $rhosts`
do
 let nprograms=0
 cputype=offline
 offline=false
 echo ========================================================================
 printf "HOST = %s " $rhost
 if [ _"$rhost" == _"$HOST" ] ; then
#     ps rhaxu OC > /tmp/$USER.psxtmp
     ps rh -eLo ruser,pid,pcpu,pmem,vsize,rss,tty,stat,start_time,time,args  > /tmp/$USER.psxtmp
     nprocs=$(grep "physical id" /proc/cpuinfo | sort -u | wc -l) # num. of phys. processors
     ncorpp=$(grep "core id" /proc/cpuinfo | sort -u | wc -l)     # num. of phys. cores per processor
     ncores=$(echo "$nprocs * $ncorpp" | bc -l)
     cputype=$(cat /proc/cpuinfo | grep -m 1 'model name' | cut -b 14- )
 else # only 1 (!) login per remote host
#     ssh $rhost "ps rhaxu OC; grep -m 1 'model name' /proc/cpuinfo | cut -b 14- ; grep -c MHz /proc/cpuinfo" > /tmp/$USER.psxtmp || continue
     ssh $rhost " ps rh -eLo ruser,pid,pcpu,pmem,vsize,rss,tty,stat,start_time,time,args ; grep -m 1 'model name' /proc/cpuinfo | cut -b 14- ; grep 'physical id' /proc/cpuinfo | sort -u | wc -l ; grep 'core id' /proc/cpuinfo | sort -u | wc -l" > /tmp/$USER.psxtmp || offline=true
     if [ "$offline"_ == "false_" ] ; then
      # necessary because of operation with let !!!
      nnlines=$(cat /tmp/$USER.psxtmp | wc -l)
      nprocs=$(tail -2  /tmp/$USER.psxtmp | head -1)
      ncorpp=$(tail -1  /tmp/$USER.psxtmp)
      ncores=$(echo "$nprocs * $ncorpp" | bc -l)
      cputype=$(tail -3 /tmp/$USER.psxtmp | head -1)
      let nlines=nnlines-3
      head -$nlines  /tmp/$USER.psxtmp >  /tmp/$USER.psxtmp2
      mv /tmp/$USER.psxtmp2  /tmp/$USER.psxtmp
     fi
 fi

# echo =================================================================
# echo "HOST = $rhost ( $cputype )"
# printf "HOST = %s " $rhost
printf "( %s )\n" "$cputype"
echo ========================================================================

 if [ "$offline"_ == "true_" ] ; then 
     echo " "
     continue
 fi

 if [ ! -e  /tmp/$USER.psxtmp ] ; then
  exit
 fi

 IFS='
'

 for line in $(cat /tmp/$USER.psxtmp) ; do
  cpuusage=$(echo $line | awk '{print $3}')
  cpuactive=$(echo "$cpuusage > 1.0" | bc -l)
  if [ "$cpuactive" == 1 ] ; then
    echo $line | awk '{printf "%-8s %8i %4s %4s %4s %6s %7s %s \n", 
                                 $1, $2, $3, $4, $8, $9, $10, $11}'
    program=$(echo $line | awk '{print $11}')
    for powrprogram in "wrstart
steal
adapter
coli
como
wrcont
formal
raytracing
montecarlo"
    do
     powrprogramactive=$(echo $line | grep "$powrprogram")   
     if [ _"$powrprogramactive" != _"" ] ; then
      let nprograms=nprograms+1
     fi 
    done
  fi  
 done
 efficiency=$(echo "100.* $nprograms / $ncores " | bc -l | awk '{printf "%4.1f",$1}')
 echo ----------------------------------------------------------------- 
 printf "$nprograms active PoWR programs, $ncores Cores available ($rhost) Efficiency = $efficiency %%  \n\n"
 let nprogramtotal=nprogramtotal+nprograms
 let ncoretotal=ncoretotal+ncores
done

if [ _$1 == _all ] ; then
 efficiency=$(echo "$nprogramtotal / $ncoretotal * 100." | bc -l | awk '{printf "%3.1f",$1}')
 printf "$nprogramtotal PoWR programs running on the cluster on $ncoretotal Cores = $efficiency %%  \n\n"
fi
