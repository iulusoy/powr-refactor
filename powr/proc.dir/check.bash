#!/bin/bash
#
# this script generates and executes the check.plot
# usage: source check.com $cmode $kn
# check modes:
# a : all -> plot for all jobnums
# h : all -> only ASCII table 
# l : last 1000 -> plot 
# c : all -> plot with colors offset for each 1000 jobs

unset usera
unalias cat
usemodel=no
cmode=$1 # check mode (see above) 
kn=$2    # chain number

# unalias cp
# -- SCAN ARGUMENTS --
if [ $# == 2 ] ; then
 # current user
 usera=$USER
 source ~/.powrconfig  # current user has bash
 POWR_WORKu=$POWR_WORK
 if [ _${POWR_WORKu} == _ ] ; then
     if [ -e  ~/.powrconfig ]  ; then
       alt_POWR_WORKu=`grep -m 1 POWR_WORK ~$usera/.powrconfig`
       # extract field left from = 
       testeq=`echo "$alt_POWR_WORKu" | grep '='`
       if [ _"$testeq" != _ ] ; then
          POWR_WORKu=`echo "$alt_POWR_WORKu" | cut -d= -f2`
       else
          # or 3rd field (csh style)
          POWR_WORKu=`echo "$alt_POWR_WORKu" | awk '{print $3}'`
       fi
       unset testeq
     else # no .powrconfig -> old default is work
       POWR_WORKu=~$usera/work/
     fi
 fi   
elif [ $# == 3 ] ; then
 source ~/.powrconfig   
 usera=$3
 # for other user: extract dir from .powrconfig (bash)
 eval powrconfig=~${usera}/.powrconfig 
 if [ -e "$powrconfig" ] ; then
   alt_POWR_WORKu=`grep -m 1 POWR_WORK $powrconfig`
   # other user has ~/ -> replace by ~user/
   testtild=`echo "$alt_POWR_WORKu" | grep '~/'`
   if [ _"$testtild" != _  ] ; then
     alt_POWR_WORKu=`echo "$alt_POWR_WORKu" | sed s+~/+~${usera}/+`
   fi
   unset testtild
   testeq=`echo "$alt_POWR_WORKu" | grep '='`
   if [ _"$testeq" != _ ] ; then
     # extract field left from = 
     POWR_WORKu=`echo "$alt_POWR_WORKu" | cut -d= -f2`
   else
     # or 3rd field (csh style)
     POWR_WORKu=`echo "$alt_POWR_WORKu" | awk '{print $3}'`
   fi
   unset testeq
 fi 
 if [ _${POWR_WORKu} == _ ] ; then
  # default directory for user without .powrconfig 
  POWR_WORKu=~$usera/work/
 fi
 eval POWR_WORKu=$POWR_WORKu
else # no further arguments: current dir   
 usemodel=yes
 if [ ! -e wruniq.out ] ; then
  if [ -e wruniq.out.gz ] ; then
   echo "Gunzipped wruniq.out.gz to /tmp/"
   cp -f wruniq.out.gz /tmp/${USER}_wruniq.out.gz
   gunzip -f /tmp/${USER}_wruniq.out.gz
   mv /tmp/wruniq.out /tmp/${USER}_wruniq.out
   wruniqfile=/tmp/${USER}_wruniq.out
   usemodel=special
  else
   echo "No argument given & Could not find input file wruniq.out: return"
   return # note that exit would also kill the terminal (because of source ...)
  fi
    else
      echo Current directory : wruniq.out
 fi 
fi
# -- CHANGE TO TMP DIRECTORY OR WRSTART --
if [ $usemodel == no ] ; then
 if [ _$POWR_INSTTYPE == _local ] 
 then 
     cd $POWR_WORKu/tmp_data/wruniq$kn 
     wruniqfile=$PWD/out
 elif [ _$POWR_INSTTYPE == _potsdam -o _$POWR_INSTTYPE == _niscluster ]
 then
     fwhere=`cat ${POWR_WORKu}/scratch/wruniq$kn/fwhere`
     cd /home/${fwhere}/tmp_data/${usera}/wruniq$kn
     wruniqfile=$PWD/out
 fi
elif [ $usemodel == yes ] ; then
 wruniqfile=$PWD/wruniq.out
fi

####### epsilon ###########

epsilon=`grep "^EPSILON" CARDS | cut -d= -f 2`
if [ _"$epsilon" == _"" ] ;  then  # epsilon was comment out in CARDS file
 epsilon=0.005
fi
logepsilon=$(perl -le "print log($epsilon)/log(10)")
logtauepsilon=$(grep " CORRLIMIT" CARDS | grep "^TAUMAX" | tail -1 | sed s/CORRLIMIT=/:/ | cut -d":" -f2 | cut -d " " -f1)
logtempepsilon1=$(grep "^NO TEMPERATURE CORRECTIONS WHILE COR" CARDS | tail -1 | sed s/.GT./:/ | cut -d: -f2 | awk '{print $1}')
logtempepsilon2=$(grep "^NO TEMPERATURE CORRECTIONS WHILE COR" CARDS | tail -1 | sed s/.GT./:/ | cut -d: -f2 | awk '{print $2}')
## {unset usera && exit}
# -- WHERE TO WRITE TMP FILES --
if [ -d "/tmp" ] ; then
 mkdir -p /tmp/${USER}
 tmppath="/tmp/${USER}"
else
 tmppath=$HOME
fi
jobfile=$tmppath/jobnum.dat
corrmfile=$tmppath/corrmax.dat
plotfile=$tmppath/corrmax.plot
inputfile=$tmppath/input.dat 
out=$tmppath/input.tmp
/bin/rm -f $out
# -- create tmp files
if [ $OSTYPE == "osf1" ] ; then
 grep CORM $wruniqfile | cut -b 92-97 > $corrmfile
 jobnrline=`cut -b 112-130 $wruniqfile | grep 'JOB NO.' | cut -d. -f2 | grep "[0-9]"`
 for line in $jobnrline ; do 
   # jobnrcheck=`expr match $line "[0-9]"`
   # if [ $jobnrcheck -gt 0 ] ; then
     echo $line >> $jobfile
   # fi
 done
else
 grep -a CORM $wruniqfile | cut -b 92-97 > $corrmfile
 jobnrline=`cut -b 117-130 $wruniqfile | grep -a "^JOB NO." | cut -d. -f2 | grep "[0-9]"`
 for line in $jobnrline ; do
     jobnrcheck=`expr $line : "[0-9]"`
   if  [ ${jobnrcheck:--99} -gt 0 ] ; then
       echo $line >> $jobfile       
   fi
 done
fi
nlines=`cat $jobfile | wc -l`
if [ _"$cmode" == _l ] ; then
 paste -d' ' $jobfile $corrmfile | tail -333 > $inputfile
else
 paste -d' ' $jobfile $corrmfile > $inputfile
fi
/bin/rm -f $jobfile $corrmfile 

# -- split the input file
if [ _"$cmode" == _c ] ; then
let oldnumber=999
let pltn=0
let skippl=0
let skipplc=0
OLDIFS=$IFS
IFS='
'
for line in `cat $inputfile` ; do
 n=`echo $line | cut -d " " -f1-1 | cut -d "." -f2-2`
 let number=$n
 let rnumber=$number
 let rnumber="$rnumber % 1000"
 if [ $pltn == 0 -a $skippl == 0 -a $number -gt $oldnumber ] ; then
   let skippl="$number / 1000"
   let skipplc=$skippl
   if [ $skipplc -gt 9 ] ; then
     let skipplc="$skipplc % 10" 
        if [ $skipplc -gt 9 ] ; then
          let skipplc="$skipplc % 10" 
        fi
   fi
 fi
 if [ $rnumber -lt  $oldnumber ] ; then
  if [ $pltn -gt 0 ] ; then
   echo 'FINISH' >> $out
  fi
  echo 'N=?' >> $out
  let pltn++
 fi
 echo "$line" >> $out
 let oldnumber=$rnumber
done
IFS=$OLDIFS
echo 'FINISH' >> $out
/bin/mv -f $out $inputfile
fi # split

if [ _"$cmode" == _a -o "_$cmode" == _c -o "_$cmode" == _l ] ; then
# -- CREATe PLOT FILE --
echo "PLOT : corrmax" > $plotfile
echo '\FONT=TIMES' >> $plotfile
# echo '\INSTRUCTION EXPORT' >> $plotfile
if  [ $nlines -lt 3 ] ; then
 sym1=1
else 
 sym1=0
fi
echo '\DEFINECOLOR=3 0.8 0.8 0.8' >> $plotfile
echo '\COLOR=3' >> $plotfile
echo '\LINUN XMIN -0.5 XMAX -0.5 0 0' >> $plotfile
echo '\LINUN XMIN -1.0 XMAX -1.0 0 0' >> $plotfile
echo '\LINUN XMIN -1.5 XMAX -1.5 0 0' >> $plotfile
echo '\LINUN XMIN -2.0 XMAX -2.0 0 0' >> $plotfile
echo '\LINUN XMIN -2.5 XMAX -2.5 0 0' >> $plotfile
echo '\RESETCOLOR=3' >> $plotfile
echo '\DEFINECOLOR=3 0.2 0.6 0.0' >> $plotfile
echo '\BGRLUN COLOR=-1'  >> $plotfile
echo '\COLOR=1'  >> $plotfile
echo '\LINUNLAB XMIN '$logepsilon' XMAX '$logepsilon' 0 0 0 0.2 EPSILON' >> $plotfile
if [ ${logtauepsilon:--99} != -99 ] ; then
 echo '\LINUNLAB XMIN '$logtauepsilon' XMAX '$logtauepsilon' 0 0 0 0.2 TAU-CORRLIMIT' >> $plotfile
fi
if [ ${logtempepsilon1:--99} != -99 ] ; then
 echo '\LINUNLAB XMIN '$logtempepsilon1' XMAX '$logtempepsilon1' 0 0 0 0.2 TEMP-CORRLIMIT START' >> $plotfile
fi
if [ ${logtempepsilon2:--99} != -99 ] ; then
 echo '\LINUNLAB XMIN '$logtempepsilon2 'XMAX '$logtempepsilon2' 0 0 0 0.2 TEMP-CORRLIMIT STOP' >> $plotfile
fi
echo '\BGRLUN OFF'  >> $plotfile
echo '\DEFINECOLOR=5 1.0 0.8 0' >> $plotfile
if [ "_$cmode" == _c ] ; then
 if [ $pltn -gt 1 ] ; then 
  let jobfirst="( $skippl ) * 1000"
  echo "\LUN XMAX YMAX R-0.5 U-.5 0.25 &2 IT1= $jobfirst" >> $plotfile
 fi
fi # split
 echo '\COLOR=1' >> $plotfile
echo "HEADER:" >> $plotfile
echo 'X-ACHSE:\CENTER\JOBNUM' >> $plotfile
echo 'Y-ACHSE:\CENTER\log CORRMAX' >> $plotfile
echo "         MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER" >> $plotfile
if [ _"$cmode" == _l ] ; then
 jobmin=`head -1 $inputfile | awk '{print  $1}'`
 jobmax=`tail -1 $inputfile | awk '{print  $1}'`
 echo "X:        20.CM        $jobmin        $jobmax       10.        100.       .0" >> $plotfile >> $plotfile
 unset jobmin jobmax
elif [  _"$cmode" == _c ] ; then    
 echo "X:        20.CM        .0        1000.       10.        100.       .0" >> $plotfile
else
 echo "X: AUTOX" >> $plotfile   
fi
echo "Y:        15.CM        -3.5       1.5        .1          .5        .0" >> $plotfile
if [ "_$cmode" == _c ] ; then
let rpl=1
#------------------------------------------------------------
while [ $rpl -le $pltn ] ; do
 let clx="1 + $rpl" 
 let roff="( $rpl + $skippl - 1 ) * 1000"
 let roffc="( $rpl + $skipplc - 1 ) * 1000"
 if [ $skipplc -gt 0 ] ; then
  let clx="$clx - ($skipplc - 1)"
  if [ $clx == 0 ] ; then
   clx=1
  fi
  echo "\LUN XMAX YMAX R-0.5 U-.5 0.25 &$clx +$roff" >> $plotfile
 fi
 if [ $clx -gt 9 ] ; then
  clx=1
 fi
 if  [ $clx == 0 ] ; then
  clx=1
 fi
 if  [ $rpl -gt 1 ] ; then
  # echo "\NEXTLUN &$clx --IT=$rpl" >> $plotfile
  echo "\NEXTLUN &$clx +$roff" >> $plotfile
 fi
 echo '*' >> $plotfile
 echo "N=? XYTABLE  COLOR=$clx" >> $plotfile
 echo "COMMAND SETNAME corrmax$rpl" >> $plotfile
 #if ( $rpl > 1 ) then
 echo "COMMAND X- $roff" >> $plotfile
 #endif
 echo "COMMAND INCLUDE $inputfile DATASET=$rpl" >> $plotfile
 echo '*' >> $plotfile
 let rpl++
done
#------------------------------------------------------------
echo "N=? XYTABLE COLOR=4 SYMBOL=$sym1" >> $plotfile
echo "COMMAND APPEND corrmax1" >> $plotfile
echo "FINISH" >> $plotfile
echo '*' >> $plotfile
else # no split
  #------------------------------------------------------------
echo "N=? XYTABLE COLOR=2" >> $plotfile
echo "COMMAND INCLUDE $inputfile" >> $plotfile
#------------------------------------------------------------
if [ "_$cmode" != _l ] ; then 
echo "N=? XYTABLE COLOR=1 SYMBOL=0" >> $plotfile
echo "1 0" >> $plotfile
echo "996 0" >> $plotfile
echo "FINISH" >> $plotfile
fi
#------------------------------------------------------------
fi # split  
echo "END" >> $plotfile

# -- CALL WRPLOT --
if [ ${PATH_WRPLOT:--99} != 99 ] ; then
 echo "found wrplot:" ${PATH_WRPLOT}
else
 echo "no wrplot"
 # export PATH_WRPLOT="/home/crater/htodt/linux-wrplot.dir"
fi 
$PATH_WRPLOT/proc.dir/wrplot.com $plotfile   || cat $inputfile
# graph -T X -y -3.5 1.5 $inputfile
#/bin/rm -f $inputfile $plotfile

elif [ "_$cmode" == _h ] ; then
 cat $inputfile
 set lasteps=`tail -1 $inputfile | awk '{print $2}'`
 set lasteps=$(perl -le "print 10**($lasteps)")
 ####### epsilon ###########
 set epsilon=`grep "^EPSILON" CARDS | cut -d= -f 2`
 set logepsilon=$(perl -le "print log($epsilon)/log(10)")
 echo "LOG EPSILON IS ${logepsilon}, last epsilon is $lasteps"

fi # end of representation branch

/bin/rm -rf $tmppath
unset usera
unset wruniqfile
unset usemodel
unset POWR_WORKu
