#!/bin/bash
# version 13.10.2020 by htodt
# copy Model data from wrdata$kn and output to recent directory
# also from different user
# usage: modsave.com 1 [-f] [-uRUSER]

# presets 
FORCE=false
RUSER=$USER
BGZIP=false # false = no gzip of model file

# at least the chain number kn must be given
if [ $# -lt 1 ]
then
 echo "Usage: modsave kn [-f] [-uRUSER]"
 exit
fi

# scan argument list, options start with - 
for arg ; do
  if [ `echo $arg | cut -c1` = - ] ; then
      if   [ `echo $arg | cut -c2` = u ] ; then
          RUSER=`echo $arg | cut -c3-`
      elif [ `echo $arg | cut -c2` = f ] ; then
          FORCE=true
      fi
  else # argument without -> is chain number
      NR=$arg
  fi
done

if [ _"$RUSER" == _"$USER" ] ; then # --------------------------
   source ~/.powrconfig 
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
else # other user give ---------------------------------------
 usera=$RUSER
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
fi # -------------------------------------------------------------
     
# do not write to non-empty directory
if [ `ls |wc -l` -ne 0 -a $FORCE == false ] ; then
  echo "Current directory $PWD is not empty: STOP ; try option -f"
  exit
fi

if [ ! -e $POWR_WORKu/wrdata$NR ] ; then
  echo "$POWR_WORKu/wrdata$NR does not exist: STOP" 
  exit
fi    

echo "Copying from ${POWR_WORKu}"
copt="-na"
cp ${copt}  $POWR_WORKu/output/formal$NR.out            formal.out
cp ${copt}  $POWR_WORKu/output/formal$NR.plot           formal.plot
cp ${copt}  $POWR_WORKu/output/wruniq$NR.out            wruniq.out
cp ${copt}  $POWR_WORKu/output/wruniq$NR.plot           wruniq.plot
cp ${copt}  $POWR_WORKu/output/wrstart$NR.out           wrstart.out
cp ${copt}  $POWR_WORKu/wrdata$NR/CARDS                 CARDS
cp ${copt}  $POWR_WORKu/wrdata$NR/FORMAL_CARDS          FORMAL_CARDS
if [ ! -e MODEL -a ! -e MODEL.gz ] ; then 
 cp ${copt} $POWR_WORKu/wrdata$NR/MODEL                 MODEL
fi 
cp ${copt}  $POWR_WORKu/wrdata$NR/NEWDATOM_INPUT        NEWDATOM_INPUT
cp ${copt}  $POWR_WORKu/wrdata$NR/NEWFORMAL_CARDS_INPUT NEWFORMAL_CARDS_INPUT
cp ${copt}  $POWR_WORKu/wrdata$NR/DATOM                 DATOM
if [ -e $POWR_WORKu/wrdata$NR/modinfo.kasdefs ] ; then
 cp ${copt} $POWR_WORKu/wrdata$NR/modinfo.kasdefs       modinfo.kasdefs
elif [ -e $POWR_WORKu/output/modinfo${NR}.kasdefs ] ; then
 cp ${copt} $POWR_WORKu/output/modinfo${NR}.kasdefs     modinfo.kasdefs        
fi
cp ${copt}  $POWR_WORKu/wrdata$NR/FGRID                 FGRID
cp ${copt}  $POWR_WORKu/wrdata$NR/FEDAT                 FEDAT
cp ${copt}  $POWR_WORKu/wrdata$NR/FEDAT_FORMAL          FEDAT_FORMAL

if [ _${RUSER} == _${USER} ] ; then 
    echo "SAVED" >> $POWR_WORKu/scratch/wruniq${NR}/name
fi

chopt="-f"
chp="444"
chmod ${chopt} ${chp}  formal.out
chmod ${chopt} ${chp}  formal.plot
chmod ${chopt} ${chp}  wruniq.out
chmod ${chopt} ${chp}  wruniq.plot
chmod ${chopt} ${chp}  wrstart.out
chmod ${chopt} ${chp}  CARDS
chmod ${chopt} ${chp}  FORMAL_CARDS
chmod ${chopt} ${chp}  MODEL
if [ $BGZIP == true ] ; then
    gzip MODEL && echo "Gzipped MODEL"
fi 
chmod ${chopt} ${chp}  DATOM
chmod ${chopt} ${chp}  NEWDATOM_INPUT
chmod ${chopt} ${chp}  NEWFORMAL_CARDS_INPUT
chmod ${chopt} ${chp}  modinfo.kasdefs
chmod ${chopt} ${chp}  FGRID

echo "Saved all data into current directory: $PWD"
unset POWR_WORKu copt chopt chp
