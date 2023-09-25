#!/bin/csh

# this script generates and executes the checkt.plot
# usage: source checkt.com $kn [$user]
#
# -- SCAN ARGUMENTS --
unset usera
set usemodel=no
set kn=$1
# -- SCAN ARGUMENTS --
if ($#argv == 1) then
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
else if ($#argv == 2) then
 set usera=$2
 source ~/.powrconfig_csh
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
else
 echo Current directory : wruniq.plot
 set usemodel=yes
endif

# -- CHANGE TO TMP DIRECTORY OR WRSTART --
if ( $usemodel == no ) then 
 if ( -d /home/`cat $POWR_WORKu/scratch/wruniq${kn}/fwhere`/tmp_data/${usera}/wruniq${kn} ) then
  cd /home/`cat $POWR_WORKu/scratch/wruniq${kn}/fwhere`/tmp_data/${usera}/wruniq${kn} 
 else if ( $POWR_INSTTYPE == local) then
  cd $POWR_TMP_DATA/wruniq${kn}
 else
  cd ${POWR_WORKu}/scratch/wrstart${kn}
  set stealpath=${POWR_WORKu}/scratch/wrstart${kn}
 endif
 set stealfile=$stealpath/steal.plot
else
  set stealfile=$PWD/wruniq.plot
endif


# -- WHERE TO WRITE TMP FILES --
if ( -d "/tmp" ) then
 mkdir -p /tmp/${USER}
 set tmppath="/tmp/${USER}"
else 
 set tmppath=$HOME
endif 
set plotfile=$tmppath/checkt.plot
set scratchfile=$tmppath/scratch.checkt
# -- CREATe PLOT FILE --
echo "PLOT:  WR TEMPERATURE STRATIFICATION T(R) VERSUS LOG(R/R*-1)" > $plotfile
echo "KASDEF COLOR=8" >> $plotfile
echo "KASDEF LINUN XMIN 10. XMAX 10.  0. 0." >> $plotfile
echo "KASDEF LINUN XMIN 20. XMAX 20.  0. 0." >> $plotfile
echo "KASDEF COLOR=1" >> $plotfile
echo "KASDEF FONT=HELVET" >> $plotfile
if ( _`grep ':UNLU' $stealfile` == _ ) then
 echo " HEADER :&2T(r)&1 only - NO TEMPERATURE CORRECTIONS APPLIED YET" >> $plotfile 
else
 echo " HEADER :T(r) &2(new)&1 compared to last iteration &4(old)&1 - 100 times enhanced" >> $plotfile  
endif
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
if ( _`grep ':UNLU' $stealfile` != _ ) then
 echo "N=? COLOR=4" >> $plotfile
 echo '***  100fach ueberhoehte Temperaturkorrektur' >> $plotfile
 echo 'COMMAND Y* -100' >> $plotfile
 echo "COMMAND ARI 2 + 1" >> $plotfile 
 echo 'COMMAND INCLUDE' $stealfile 'INCKEY="PLOT   :UNLU" DATASET=1' >> $plotfile
endif
echo "" >> $plotfile
echo '*** um YMIN auf 0 zu zwingen:' >> $plotfile
echo "N=1 XYTABLE SYMBOL=0" >> $plotfile
echo "0 0" >> $plotfile
echo "" >> $plotfile
# plot the TMIN line (not for OSF1)
if ($OSTYPE != "osf1") then
 grep "TMIN" -B 3 $stealfile | sed s/L0/R-0/ >> $plotfile
 echo "\COLOR = 1" >> $plotfile
endif
echo "END" >> $plotfile
# -- CALL WRPLOT --
if ( ${?PATH_WRPLOT} ) then
              :
 echo "Found wrplot."
else
 echo "No wrplot path set. I fix it for linux."
 setenv PATH_WRPLOT "/home/crater/htodt/linux-wrplot.dir"
 alias wrplot "$PATH_WRPLOT/proc.dir/wrplot.com"
endif 
wrplot $plotfile 
/bin/rm -f $plotfile $scratchfile
/bin/rm -rf $tmppath

unset usera
unset usemodel
