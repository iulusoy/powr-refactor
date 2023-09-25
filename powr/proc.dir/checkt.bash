#!/bin/bash
source ~/.powrconfig
#
# this script generates and executes the checkt.plot
# usage: source checkt.com $kn [$user]
#
# -- SCAN ARGUMENTS --
unset usera
usemodel=no
if [ $# == 1 ] ; then
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
fi

# -- CHANGE TO TMP DIRECTORY OR WRSTART --
if [ $usemodel == no ] ; then
  if [ -d /home/`cat $POWR_WORKu/scratch/wruniq${kn}/fwhere`/tmp_data/${usera}/wruniq${kn} ]; then
    stealpath=/home/`cat $POWR_WORKu/scratch/wruniq${kn}/fwhere`/tmp_data/${usera}/wruniq${kn}
    cd $stealpath
  elif [ $POWR_INSTTYPE == local ]; then
    cd $POWR_WORKu/tmp_data/wruniq$1 
    stealpath=$POWR_WORKu/tmp_data/wruniq$1
  else
    cd ${POWR_WORKu}/scratch/wrstart${kn}
    stealpath=${POWR_WORKu}/scratch/wrstart${kn}
  fi
  stealfile=$stealpath/steal.plot
else
  stealfile=$PWD/wruniq.plot
fi
# -- WHERE TO WRITE TMP FILES --
if [ -d "/tmp" ] ; then
 mkdir -p /tmp/${USER}
 tmppath="/tmp/${USER}"
else 
 tmppath=$HOME
fi
plotfile=$tmppath/checkt.plot
scratchfile=$tmppath/scratch.checkt
# -- CREATe PLOT FILE --
echo "PLOT:  WR TEMPERATURE STRATIFICATION T(R) VERSUS LOG(R/R*-1)" > $plotfile
echo "KASDEF COLOR=8" >> $plotfile
echo "KASDEF LINUN XMIN 10. XMAX 10.  0. 0." >> $plotfile
echo "KASDEF LINUN XMIN 20. XMAX 20.  0. 0." >> $plotfile
echo "KASDEF COLOR=1" >> $plotfile
echo "KASDEF FONT=HELVET" >> $plotfile
if [ _$(grep ':UNLU' $stealfile) == _ ] ; then
 echo " HEADER :&2T(r)&1 only - NO TEMPERATURE CORRECTIONS APPLIED YET" >> $plotfile 
else
 echo " HEADER :T(r) &2(new)&1 compared to last iteration &4(old)&1 - 100 times enhanced" >> $plotfile  
fi
echo " X-ACHSE:\CENTER\Depth Index L" >> $plotfile
echo " Y-ACHSE:\CENTER\T / kK" >> $plotfile
echo "     MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER" >> $plotfile
echo " X:  AUTO 0.00000     -5.00000      1.50000     0.500000      1.00000      0.00000" >> $plotfile
echo " Y:  AUTO 0.00000      0.00000      100.000      10.0000      10.0000      0.00000" >> $plotfile
echo "" >> $plotfile
echo "N= ?   PLOTSYMBOL=0" >> $plotfile
echo "COMMAND WRITE FILE=$scratchfile" >> $plotfile
echo "COMMAND SKIP" >> $plotfile
echo 'COMMAND INCLUDE' $stealfile 'INCKEY=" PLOT   : WR TEMPERATURE STRATIFICATION"' >> $plotfile
echo "" >> $plotfile
echo "" >> $plotfile
echo "N=? XYTABLE SELECT=-1,2 COLOR=2" >> $plotfile
echo "COMMAND XY-SWAP" >> $plotfile
echo "COMMAND CF-YMAX dummy ND" >> $plotfile
echo "COMMAND XY-SWAP" >> $plotfile
echo 'COMMAND X* -1' >> $plotfile
echo 'COMMAND X+ $ND' >> $plotfile
echo "COMMAND INCLUDE $scratchfile" >> $plotfile
echo "" >> $plotfile
echo '*' "stealfile is $stealfile" >> $plotfile
if [ _$(grep ':UNLU' $stealfile) == _ ] ;  then
  echo '* UNLU not found in ' $stealfile >> $plotfile 
else
 echo "N=? COLOR=4" >> $plotfile
 echo '***  100times increased temperature correction' >> $plotfile
 echo 'COMMAND Y* -100' >> $plotfile
 echo "COMMAND ARI 2 + 1" >> $plotfile
 echo 'COMMAND INCLUDE' $stealfile 'INCKEY="PLOT   :UNLU" DATASET=1' >> $plotfile
fi
echo "" >> $plotfile
echo '*** to force YMIN to 0:' >> $plotfile
echo "N=1 XYTABLE SYMBOL=0" >> $plotfile
echo "0 0" >> $plotfile
echo "" >> $plotfile
# plot the TMIN line (not for OSF1)
if  [ $OSTYPE != "osf1" ] ; then
 grep "TMIN" -B 3 $stealfile | sed s/L0/R-0/ >> $plotfile
 echo "\COLOR = 1" >> $plotfile
fi
echo "END" >> $plotfile
# -- CALL WRPLOT --
if [ ${?PATH_WRPLOT} ] ; then
              :
 echo "Found wrplot."
else
 echo "No wrplot path set. I fix it for linux."
fi
$PATH_WRPLOT/proc.dir/wrplot.com $plotfile 
# /bin/rm -f $plotfile $scratchfile
# /bin/rm -rf $tmppath

unset usera
unset usemodel
