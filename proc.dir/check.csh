#!/bin/csh
# 
# this script generates and executes the check.plot
# usage: source check.com $cmode $kn [$user]
# check modes:
# a : all -> plot for all jobnums
# h : all -> only ASCII table 
# l : last 1000 -> plot 
# c : all -> plot with colors offset for each 1000 jobs

unset usera
set usemodel=no
set cmode=$1 # check mode (see above) 
set kn=$2    # chain number
# -- SCAN ARGUMENTS --
if ($#argv == 2) then
 # current user
 set usera=$USER
 source ~/.powrconfig_csh # current user has tcsh
 set POWR_WORKu=$POWR_WORK
 if (${?POWR_WORKu} == 0) then
     if (`filetest -e  ~/.powrconfig` == 1 ) then
       set alt_POWR_WORKu=`grep -m 1 POWR_WORK ~$usera/.powrconfig`
       # extract field left from = 
       if ( _`echo "$alt_POWR_WORKu" | grep '='` != _ ) then
          set POWR_WORKu=`echo "$alt_POWR_WORKu" | cut -d= -f2`
       else
          # or 3rd field (csh style)
          set POWR_WORKu=`echo "$alt_POWR_WORKu" | awk '{print $3}'`
       endif
     else
       set POWR_WORKu=~$usera/work/
     endif   
 endif
 # POWR_WORKu is still empty:
 if ( _$POWR_WORKu == _ ) then
   # default directory for user without .powrconfig 
   set POWR_WORKu=~$usera/work/
 endif 
else if ($#argv == 3) then
 source ~/.powrconfig_csh
 set usera=$3
 # for another user: extract dir from .powrconfig (bash)
 if (`filetest -e ~$usera/.powrconfig` == 1) then
   set alt_POWR_WORKu=`grep -m 1 POWR_WORK ~$usera/.powrconfig`
   # other user has ~/ -> replace by ~user/
   if ( _`echo "$alt_POWR_WORKu" |grep '~/'` != _ ) then
     set alt_POWR_WORKu=`echo "$alt_POWR_WORKu" | sed s+~/+~$usera/+`
   endif
   if ( _`echo "$alt_POWR_WORKu" | grep '='` != _ ) then
     # extract field left from = 
     set POWR_WORKu=`echo "$alt_POWR_WORKu" | cut -d= -f2`
   else
     # or 3rd field (csh style)
     set POWR_WORKu=`echo "$alt_POWR_WORKu" | awk '{print $3}'`
   endif
 endif
 if ( ${?POWR_WORKu} == 0 ) then
   # default directory for user without .powrconfig 
  set POWR_WORKu=~$usera/work/
 endif
else # no further arguments: current dir
 set usemodel=yes
 if ( ! -e wruniq.out ) then
  if ( -e wruniq.out.gz ) then
   echo "Gunzipped wruniq.out.gz to /tmp/"
   cp -f wruniq.out.gz /tmp/${USER}_wruniq.out.gz
   gunzip -f /tmp/${USER}_wruniq.out.gz
   mv /tmp/wruniq.out /tmp/${USER}_wruniq.out
   set wruniqfile=/tmp/${USER}_wruniq.out
   set usemodel=special
  else
   echo "No argument given & Could not find input file wruniq.out: exit"
   exit
  endif
 else
   echo Current directory : wruniq.out
 endif 
endif

# -- CHANGE TO TMP DIRECTORY OR WRSTART --
if ( $usemodel == no ) then
 if ( -d /home/`cat $POWR_WORKu/scratch/wruniq${kn}/fwhere`/tmp_data/${usera}/wruniq${kn} ) then
  cd /home/`cat $POWR_WORKu/scratch/wruniq${kn}/fwhere`/tmp_data/${usera}/wruniq${kn} 
 else if ( $POWR_INSTTYPE == local) then
  cd $POWR_TMP_DATA/wruniq${kn}
 else # allow check for tmp-hosts that are not NFS-mounted
  #   cd ~${usera}/powr/scratch/wrstart${kn}
  cd $POWR_WORKu/scratch/wruniq${kn}
  set rtmphost=`cat $POWR_WORKu/scratch/wruniq${kn}/fwhere`
  scp ${rtmphost}:/home/${rtmphost}/tmp_data/${usera}/wruniq${kn}/out . || echo '===== YOU CANNOT WRITE IN CURRENT DIR ===='
  scp ${rtmphost}:/home/${rtmphost}/tmp_data/${usera}/wruniq${kn}/steal.plot . 
  scp ${rtmphost}:/home/${rtmphost}/tmp_data/${usera}/wruniq${kn}/steal_last.plot .
  scp ${rtmphost}:/home/${rtmphost}/tmp_data/${usera}/wruniq${kn}/como.plot .
  /bin/cp -f $POWR_WORKu/wrdata${kn}/CARDS .  
 endif
 set wruniqfile=$PWD/out
else if ( $usemodel == yes ) then
 set wruniqfile=$PWD/wruniq.out
endif
####### epsilon ###########

set epsilon=`grep "^EPSILON" CARDS | cut -d= -f 2`
if ( _"$epsilon" == _"" ) then  # epsilon was comment out in CARDS file
 set epsilon=0.005
endif
set logepsilon=`perl -le "print log($epsilon)/log(10)"`
set logtauepsilon=`grep " CORRLIMIT" CARDS | grep "^TAUMAX" | tail -1 | sed s/CORRLIMIT=/:/ | cut -d":" -f2 | cut -d " " -f1 || echo -99`
set logtempepsilon1=`grep "^NO TEMPERATURE CORRECTIONS WHILE COR" CARDS | tail -1 | sed s/.GT./:/ | cut -d: -f2 | awk '{print $1}' || echo -99`
set logtempepsilon2=`grep "^NO TEMPERATURE CORRECTIONS WHILE COR" CARDS | tail -1 | sed s/.GT./:/ | cut -d: -f2 | awk '{print $2}'  || echo -99`

## {unset usera && exit}
# -- WHERE TO WRITE TMP FILES --
if ( -d "/tmp" ) then
 mkdir -p /tmp/${USER}
 set tmppath="/tmp/${USER}"
else
 set tmppath=$HOME
endif 
set jobfile=$tmppath/jobnum.dat
set corrmfile=$tmppath/corrmax.dat
set plotfile=$tmppath/corrmax.plot
set inputfile=$tmppath/input.dat 
set out=$tmppath/input.tmp
/bin/rm -f $out
# -- create tmp files
if ($OSTYPE == "osf1") then
 grep CORM $wruniqfile | cut -b 92-97 > $corrmfile
 set jobnrline=`cut -b 112-130 $wruniqfile | grep 'JOB NO.' | cut -d. -f2 | grep "[0-9]"`
 foreach line ($jobnrline)
   set jobnrcheck=`expr $line : "[0-9]"`
   if ($jobnrcheck > 0) then
     echo $line >> $jobfile
   endif
 end
else
 grep -a CORM $wruniqfile | cut -b 92-97 > $corrmfile
 set jobnrline=`cut -b 117-130 $wruniqfile | grep -a "^JOB NO." | cut -d. -f2 | grep "[0-9]"`
 foreach line ($jobnrline)
   set jobnrcheck=`expr $line : "[0-9]"`
   if ($jobnrcheck > 0) then
     echo $line >> $jobfile
   endif
 end
endif
set nlines=`cat $jobfile | wc -l`
if ( "_$cmode" == _l ) then
 paste -d' ' $jobfile $corrmfile | tail -333 >  $inputfile
else
 paste -d' ' $jobfile $corrmfile >  $inputfile
endif
/bin/rm -f $jobfile $corrmfile

# -- split the input file ----------
if ( "_$cmode" == _c ) then
set oldnumber=999
set pltn=0
set skippl=0
set skipplc=0
foreach line ("`cat $inputfile`")
 set n=`echo $line | cut -d " " -f1-1 | cut -d "." -f2-2`
 set number=$n
 set rnumber = $number
 @ rnumber = $rnumber % 1000
 if (($pltn == 0) && ($skippl == 0) && ($number > $oldnumber)) then
   @ skippl = $number / 1000
   set skipplc=$skippl
   if ($skipplc > 9) then
     @ skipplc = $skipplc % 10 
        if ($skipplc > 9) then
          @ skipplc = $skipplc % 10 
        endif
   endif
 endif
 if ($rnumber < $oldnumber) then
  if ($pltn > 0) then
   echo 'FINISH' >> $out
  endif
  echo 'N=?' >> $out
  @ pltn++
 endif
 echo "$line" >> $out
 set oldnumber=$rnumber
end
echo 'FINISH' >> $out
/bin/mv -f $out $inputfile
endif # split 
# ----------------------


if ( _"$cmode" == _a || "_$cmode" == _c || "_$cmode" == _l ) then
# -- CREATe PLOT FILE --
echo "PLOT : corrmax" > $plotfile
echo '\FONT=TIMES' >> $plotfile
# echo '\INSTRUCTION EXPORT' >> $plotfile
if ($nlines < 3) then
 set sym1=1
else 
 set sym1=0
endif
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
if ( $logtauepsilon != -99 ) then
 echo '\LINUNLAB XMIN '$logtauepsilon' XMAX '$logtauepsilon' 0 0 0 0.2 TAU-CORRLIMIT' >> $plotfile
endif
if ( $logtempepsilon1 != -99 ) then
 echo '\LINUNLAB XMIN '$logtempepsilon1' XMAX '$logtempepsilon1' 0 0 0 0.2 TEMP-CORRLIMIT START' >> $plotfile
endif
if ( $logtempepsilon2 != -99 ) then
 echo '\LINUNLAB XMIN '$logtempepsilon2' XMAX '$logtempepsilon2' 0 0 0 0.2 TEMP-CORRLIMIT STOP' >> $plotfile
endif
echo '\BGRLUN OFF'  >> $plotfile
echo '\DEFINECOLOR=5 1.0 0.8 0' >> $plotfile
if ( "_$cmode" == _c ) then
 if ( $pltn > 1 ) then 
  @ jobfirst = ( $skippl ) * 1000
  echo "\LUN XMAX YMAX R-0.5 U-.5 0.25 &2 IT1= $jobfirst" >> $plotfile
 endif
endif # split
echo '\COLOR=1' >> $plotfile
echo "HEADER:" >> $plotfile
echo 'X-ACHSE:\CENTER\JOBNUM' >> $plotfile
echo 'Y-ACHSE:\CENTER\log CORRMAX' >> $plotfile
echo "         MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER" >> $plotfile
if ( _"$cmode" == _l )  then
 set jobmin=`head -1 $inputfile | awk '{print  $1}'`
 set jobmax=`tail -1 $inputfile | awk '{print  $1}'`
 echo "X:        20.CM        $jobmin        $jobmax       10.        100.       .0" >> $plotfile 
 unset jobmin jobmax
else if ( _"$cmode" == _c ) then    
 echo "X:        20.CM        .0        1000.       10.        100.       .0" >> $plotfile
else
 echo "X: AUTOX" >> $plotfile   
endif
echo "Y:        15.CM        -3.5       1.5        .1          .5        .0" >> $plotfile
if ( "_$cmode" == _c ) then
set rpl=1
#------------------------------------------------------------
while ( $rpl <= $pltn )
@ clx = 1 + $rpl 
@ roff = ( $rpl + $skippl - 1 ) * 1000
@ roffc = ( $rpl + $skipplc - 1 ) * 1000
if ($skipplc > 0) then
  @ clx = $clx - ($skipplc - 1)
  if ($clx == 0) then
   set clx = 1
  endif
  echo "\LUN XMAX YMAX R-0.5 U-.5 0.25 &$clx +$roff" >> $plotfile
endif
if ($clx > 9) then
  set clx = 1
endif
if ($clx == 0) then
  set clx = 1
endif
if ($rpl > 1) then
# echo "\NEXTLUN &$clx --IT=$rpl" >> $plotfile
 echo "\NEXTLUN &$clx +$roff" >> $plotfile
endif
echo '*' >> $plotfile
echo "N=? XYTABLE  COLOR=$clx" >> $plotfile
echo "COMMAND SETNAME corrmax$rpl" >> $plotfile
#if ( $rpl > 1 ) then
echo "COMMAND X- $roff" >> $plotfile
#endif
echo "COMMAND INCLUDE $inputfile DATASET=$rpl" >> $plotfile
echo '*' >> $plotfile
@ rpl++
end # while loop for splitted data sets
-----------------------------------------------------------
echo "N=? XYTABLE COLOR=4 SYMBOL=$sym1" >> $plotfile
echo "COMMAND APPEND corrmax1" >> $plotfile
echo "FINISH" >> $plotfile
echo '*' >> $plotfile
else  # no split
#------------------------------------------------------------
echo "N=? XYTABLE COLOR=2" >> $plotfile
echo "COMMAND INCLUDE $inputfile" >> $plotfile
#------------------------------------------------------------
if ( "_$cmode" != _l ) then 
echo "N=? XYTABLE COLOR=1 SYMBOL=0" >> $plotfile
echo "1 0" >> $plotfile
echo "996 0" >> $plotfile
echo "FINISH" >> $plotfile
endif
#------------------------------------------------------------
endif # split
echo "END" >> $plotfile
# -- CALL WRPLOT --
if ( ${?PATH_WRPLOT} ) then
 echo "found wrplot:" ${PATH_WRPLOT}
else
 echo "no wrplot"
 setenv PATH_WRPLOT "/home/crater/htodt/linux-wrplot.dir"
 alias wrplot "$PATH_WRPLOT/proc.dir/wrplot.com"
endif 
wrplot $plotfile   || cat $inputfile

else if ( "_$cmode" == _h ) then
 cat $inputfile
 set lasteps=`tail -1 $inputfile | awk '{print $2}'`
 set lasteps=`perl -le "print 10**($lasteps)"`
 ####### epsilon ###########
 set epsilon=`grep "^EPSILON" CARDS | cut -d= -f 2`
 set logepsilon=`perl -le "print log($epsilon)/log(10)"`
 echo "LOG EPSILON IS ${logepsilon}, last epsilon is $lasteps"

endif # end of representation branch

# graph -T X -y -3.5 1.5 -X 'JOBNUM' -Y 'LOG CORRMAX' $inputfile
# -T X : X11 window
# -y ymin ymax : upper and lower limit of y-axis

#/bin/rm -f $inputfile $plotfile

/bin/rm -rf $tmppath
unset usera
unset wruniqfile
unset usemodel
unset POWR_WORKu
