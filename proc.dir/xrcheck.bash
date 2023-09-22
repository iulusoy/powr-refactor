#!/bin/bash
source ~/.powrconfig
#
# this script generates and executes emergentflux.plot
# usage: source check.com $kn [$user] [+]
# optin + calls wrps
# -- SCAN ARGUMENTS --

unset usera
usemodel=no
postscript=false
currdate=`date +"%y%m%d%H%M%S"`
if [ $# == 1 ] ; then
 if [ $1 == + ] ; then
  postscript=true
  echo Current directory : wruniq.plot
  usemodel=yes
 else
  kn=$1
 fi 
 usera=$USER
elif [ $# == 2 ] ; then
 kn=$1
 if [ $2 == + ] ; then
  postscript=true
  usera=$USER
 else 
  usera=$2
 fi
elif [ $# == 3 ] ; then
 kn=$1
 usera=$2
 if [ $3 == + ] ; then
  postscript=true
 fi
elif [ $# == 0 ] ; then
 echo Current directory : wruniq.plot
 usemodel=yes
fi

# -- CHANGE TO TMP DIRECTORY OR WRSTART --
if [ $usemodel == no ] ; then
 if [ -d $POWR_TMP_DATA/wruniq${kn} ] ; then
  cd $POWR_TMP_DATA/wruniq${kn} 
 else
  cd $POWR_WORK/scratch/wrstart${kn}
  rtmphost=`cat $POWR_WORK/scratch/wruniq$1/fwhere`
  scp ${rtmphost}:$POWR_TMP_DATA/wruniq$1/out .
  scp ${rtmphost}:$POWR_TMP_DATA/wruniq$1/steal.plot .
  scp ${rtmphost}:$POWR_TMP_DATA/wruniq$1/como.plot .
  /bin/cp -f $POWR_WORK/powr/wrdata$1/CARDS .  
 fi
 wruniqpath=$PWD
 wruniqfile=/steal.plot
else # usemodel = yes 
 wruniqpath=$PWD
 wruniqfile=/wruniq.plot
fi
## {unset usera && exit}
# -- WHERE TO WRITE TMP FILES --
if [ -d "/tmp" ] ; then
 tmppath="/tmp"
else
 tmppath=$HOME
fi 
plotname=emergentflux$USER.plot
plotfile=$tmppath/$plotname
/bin/rm -f $plotfile
# -- CREATe PLOT FILE --
echo 'PAPERFORMAT A4Q' >> $plotfile
echo 'PLOT   :EMERGENT FLUX ' >> $plotfile    
echo '\FONT=HELVET' >> $plotfile
echo '\INBOX' >> $plotfile
echo '\OFS 2.5 2' >> $plotfile
# echo '\FONT=HELVET' >> $plotfile
echo '\INSTRUCTION EXPORT' >> $plotfile
echo '' >> $plotfile
echo '\VAR xrstart = 1.' >> $plotfile
# echo '\VAR xrend = 41.' >> $plotfile
echo '\VAR xrend = 50.' >> $plotfile
echo '\VAR xrend2 = 25.' >> $plotfile
echo '\VAR xrend3 = 124.' >> $plotfile
echo "\VAR MODEL   = $wruniqpath" >> $plotfile
echo '\EXPR FNAME = $MODEL' "// $wruniqfile" >> $plotfile
echo '' >> $plotfile
echo '' >> $plotfile
echo '\INSTRUCTION NO-EXPORT' >> $plotfile
echo '\LUN XMIN YMIN 0.5 0.5 0.2 flux from: $FNAME' >> $plotfile
echo '' >> $plotfile
echo '\LINUNLAB LOG3.099 -17 LOG62. -17 0 0 0 0.2 Einstein' >> $plotfile
echo '\LINUNLAB LOG4.1 -17.5 LOG160. -17.5 0 0 0 0.2 ROSAT' >> $plotfile
echo '\COLOR=4' >> $plotfile
echo '\PEN=10' >> $plotfile
echo '\LINUNLAB LOG$xrstart -15.5 LOG$xrend -15.5 0 0 0 0.4 Xray flux' >> $plotfile
echo '\PEN=5' >> $plotfile
echo '\LINUNLAB LOG$xrstart -16.5 LOG$xrend2 -16.5 0 0 0 0.3 Xray flux2' >> $plotfile
echo '\LINUNLAB LOG$xrstart -18 LOG$xrend3 -18 0 0 0 0.3 Xray flux3' >> $plotfile
echo '\COLOR=1' >> $plotfile
echo '\PEN=1' >> $plotfile
echo '' >> $plotfile
echo '' >> $plotfile
echo '\IDY 10' >> $plotfile
echo '\IDSIZE 0.2' >> $plotfile
echo '\IDLOG' >> $plotfile
# these wavelenths are apparently wrong
#echo '\IDENT        8.165 &EK-SHELL O 3' >> $plotfile
#echo '\IDENT       10.299 &EK-SHELL O 2' >> $plotfile
#echo '\IDENT       13.856 &EK-SHELL O 2' >> $plotfile
echo '\IDENT       19.101 &EK-SHELL O 5' >> $plotfile
# echo '\IDENT       23.408 &EBOUND-FREE O VI 22S' >> $plotfile 
echo '\IDY 12' >> $plotfile
echo '\IDENT      227.815 &4&EBOUND-FREE He II....1' >> $plotfile
echo '\IDENT      243.852 &EBOUND-FREE G  IV....2' >> $plotfile
echo '\IDENT      258.880 &EBOUND-FREE C 32S1S..1' >> $plotfile
echo '\IDENT      261.352 &EBOUND-FREE N III2P2.1' >> $plotfile
echo '\IDENT      299.596 &EBOUND-FREE C 32P3P..2' >> $plotfile
echo '\IDENT      352.725 &EBOUND-FREE O 22P4S..1' >> $plotfile
echo '\IDENT      504.209 &4&EBOUND-FREE HeI 1S1..1' >> $plotfile
echo '\IDENT      911.672 &4&EBOUND-FREE H I......1' >> $plotfile
echo '\IDENT     1424.035 &ETHOMSON    ELECTRON' >> $plotfile
echo '\IDENT     3646.688 &4&EBOUND-FREE H I......2' >> $plotfile
echo '' >> $plotfile
echo ' HEADER :&EEMERGENT FLUX OF MODEL: $MODEL' >> $plotfile
echo ' X-ACHSE:\CENTER\log (#l# / \A)' >> $plotfile
echo ' Y-ACHSE:\CENTER\log  F&T#l#&M / (erg cm&H-2&M s&H-1&M \A&H-1&M)  in 10pc' >> $plotfile
echo '     MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER' >> $plotfile
echo ' X:  26.CM      .00000      5.00000     0.100000     0.500000      0.00000' >> $plotfile
echo ' Y:  17.CM     -20.0000     -3.00000      1.00000      5.00000      0.00000' >> $plotfile
echo '' >> $plotfile
echo 'N=?   PLOTSYMBOL=  5' >> $plotfile
echo 'COMMAND SETNAME wrcontflux' >> $plotfile
echo 'COMMAND INCLUDE $FNAME INCKEY=" PLOT: EMERGENT" DATASET=1 ' >> $plotfile
echo '' >> $plotfile
echo 'N=? COLOR=2' >> $plotfile
echo 'COMMAND SETNAME coliflux' >> $plotfile
echo 'COMMAND INCLUDE $FNAME INCKEY=" PLOT: EMERGENT" DATASET=2 ' >> $plotfile
echo '' >> $plotfile
echo 'N=? COLOR=4 SYMBOL=0' >> $plotfile
echo 'COMMAND APPEND coliflux' >> $plotfile
echo 'COMMAND XDEX' >> $plotfile
echo 'COMMAND YDEX' >> $plotfile
echo 'COMMAND X-CUT $xrstart $xrend3' >> $plotfile
echo 'COMMAND CF-INT xrayflux3' >> $plotfile
echo 'COMMAND X-CUT $xrstart $xrend' >> $plotfile
echo 'COMMAND CF-INT xrayflux' >> $plotfile
echo 'COMMAND X-CUT $xrstart $xrend2' >> $plotfile
echo 'COMMAND CF-INT xrayflux2' >> $plotfile
echo 'FINISH' >> $plotfile
echo '' >> $plotfile
echo 'N=? COLOR=4 SYMBOL=0' >> $plotfile
echo 'COMMAND APPEND coliflux' >> $plotfile
echo 'COMMAND XDEX' >> $plotfile
echo 'COMMAND YDEX' >> $plotfile
echo 'COMMAND X-CUT $xrend 1.E20' >> $plotfile
echo 'COMMAND CF-INT bolflux' >> $plotfile
echo 'FINISH' >> $plotfile
echo '' >> $plotfile
echo '\CALC xoverbol = $xrayflux / $bolflux' >> $plotfile
echo '\CALC xoverbol2 = $xrayflux2 / $bolflux' >> $plotfile
echo '\CALC xoverbol3 = $xrayflux3 / $bolflux' >> $plotfile
echo '' >> $plotfile
echo '\CALC Lx = $xrayflux * 4 * 3.1416 * 10. * 10. * (3.09E18)*(3.09E18)' >> $plotfile
echo '\CALC Lx2 = $xrayflux2 * 4 * 3.1416 * 10. * 10. * (3.09E18)*(3.09E18)' >> $plotfile
echo '\CALC Lx3 = $xrayflux3 * 4 * 3.1416 * 10. * 10. * (3.09E18)*(3.09E18)' >> $plotfile
echo '\CALC logLx = log ( $Lx / 3.846E33 )' >> $plotfile
echo '' >> $plotfile
echo '\CALC Lbol = $bolflux * 4 * 3.1416 * 10. * 10. * (3.09E18)*(3.09E18)' >> $plotfile
echo '\CALC logLbol = log ( $Lbol / 3.846E33 )' >> $plotfile
echo '\CALC logxoverbol = LOG($xoverbol)' >> $plotfile
echo '\CALC logxoverbol2 = LOG($xoverbol2)' >> $plotfile
echo '\CALC logxoverbol3 = LOG($xoverbol3)' >> $plotfile
echo '' >> $plotfile
echo '\VAR-LIST' >> $plotfile
echo '\FORMAT (1PG10.2) $xoverbol' >> $plotfile
echo '\FORMAT (1PG10.3) $logxoverbol' >> $plotfile
echo '\FORMAT (1PG12.2) $bolflux' >> $plotfile
echo '\FORMAT (1PG12.2) $xrayflux' >> $plotfile
echo '\FORMAT (1PG12.2) $Lx' >> $plotfile
echo '\FORMAT (1PG12.2) $Lx2' >> $plotfile
echo '\FORMAT (1PG12.2) $Lx3' >> $plotfile
echo '\FORMAT (1PG12.2) $Lbol' >> $plotfile
echo '\FORMAT (1PG12.4) $logLbol' >> $plotfile
echo '\FORMAT (1PG12.4) $logLx' >> $plotfile
echo '' >> $plotfile
echo '\LUN XMIN YMAX L0.5 U-1 0.2 X-ray:  F&T#l#&Md#l# &4$xrstart &1- &4$xrend &1\A =  &2 $xrayflux &1erg\,s&H-1&M\,cm&H-2&M' >> $plotfile
echo '\NEXTLUN L&TX&M  = $Lx erg\,s&H-1&M  (log L&TX&M/L\S  = $logLx)' >> $plotfile
echo '\NEXTLUN L&TX2&M = $Lx2 erg\,s&H-1&M ( $xrstart - $xrend2 \A )' >> $plotfile
echo '\NEXTLUN bol: F&T#l#&Md#l# $xrend - ... \A =              &2 $bolflux &1erg\,s&H-1&M\,cm&H-2&M' >> $plotfile
echo '\NEXTLUN L&Tbol&M = $Lbol erg\,s&H-1&M (log L&Tbol&M/L\S  = $logLbol)' >> $plotfile
echo '\NEWSIZELUN 0.35' >> $plotfile
echo '\PEN = 7' >> $plotfile
echo '\NEXTLUN L&TX&M/L&Tbol&M = &2$xoverbol &1(log = &4$logxoverbol )' >> $plotfile
echo '\PEN = 3' >> $plotfile
echo '\LUN XMAX -15 R-1 0. 0.2 $xrstart -  $xrend2 \,\A : L&TX&M/L&Tbol&M = &2$xoverbol2 &1(log = &4$logxoverbol2 )' >> $plotfile
echo '\NEXTLUN $xrstart - $xrend3 \,\A : L&TX&M/L&Tbol&M = &2$xoverbol3 &1(log = &4$logxoverbol3 )' >> $plotfile
echo '\LUNINC XMAX YMAX R-0.5 U-1 0.15 CARDS "" XRAY' >> $plotfile
echo '' >> $plotfile
echo 'ENDE' >> $plotfile
if [ $postscript == true ] ; then
  cp CARDS /tmp
  cd /tmp
  $PATH_WRPLOT/proc.dir/wrps.com $plotfile #  || cat $inputfile
  cp ${plotfile:r}.ps $HOME/${plotname:r}_${currdate}.ps
  echo "====================================================================="
  echo "Moved ${plotfile:r}.ps to $HOME/${plotname:r}_${currdate}.ps"
  echo "====================================================================="
  /bin/rm -f $plotfile ${plotfile:r}.ps CARDS
  cd -
else
 $PATH_WRPLOT/proc.dir/wrplot.com $plotfile #  || cat $inputfile
fi
# graph -T X -y -3.5 1.5 $inputfile
# /bin/rm -f $plotfile

unset usera
unset wruniqfile
unset usemodel
