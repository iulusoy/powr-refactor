#!/bin/bash
# creates NEWFORMAL_CARDS_INPUT based on current FORMAL_CARDS 
# Detailed implementation that constructs the file from
# actual FORMAL_CARDS entries
# Andreas Sander, last update 15 Oct 2021

fc='FORMAL_CARDS'
nfc='NEWFORMAL_CARDS_INPUT'
rangedata='ranges.dummy'
cmddata='commands.dummy'
labdata='labs.dummy'

if [ -e $fc ] ; then
  if [ -e $nfc ] ; then
    mv $nfc NEWFORMAL_CARDS_INPUT_OLD
  fi  
  stdpath=`grep -A1 Standardpfad FORMAL_CARDS | tail -1 | cut -b3-`
  if [[ $stdpath != "" ]]; then
    echo 'STANDARDPATH '$stdpath >> $nfc
  else 
    stdpath=`grep STANDARDPATH FORMAL_CARDS | tail -1 | cut -b3-`
    echo $stdpath >> $nfc
  fi  
  
  # We allow an arbitrary number of spaces before the RANGE keyword, but no other characters
  grep -n '^[[:space:]]*RANGE' $fc > "$rangedata"
  curstart=1
  echo 'Creating NEWFORMAL_CARDS_INPUT with the following ranges: '
  rangeout=""
  
  while read rangeline 
  do 
    curend=`echo $rangeline | awk -F':' '{print $1}'`
    rangestr=`echo $rangeline | awk -F':' '{print $2}'`
    rangeblue=`echo $rangestr | awk '{print $2}'`
    rangered=`echo $rangestr | awk '{print $3}'`
    # get end of current BLEND block
    nextstart=`sed -n $curend,/-BLEND/= FORMAL_CARDS | tail -1`
    
    curend=$((curend-1))
    nextstart=$((nextstart+1))
    
    curlabs=''
    echo '' >> $nfc

    # check all lines between end of BLEND block and next RANGE command
    if [ -e "$cmddata" ] ; then
      rm -f $cmddata;
    fi
    sed -n ${curstart},${curend}p FORMAL_CARDS | grep -v '^*' > "$cmddata"
    while read cmdline 
    do    
      linewc=`echo $cmdline | wc -w`      
      if [ $linewc != 0 ]; then
        labcheck=`echo $cmdline | awk '{print ($1$2 == "STRINGCOMMENT")}'`
        if [ $labcheck != 0 ]; then
          # Range name keywords are stored
          newlab=`echo "$cmdline" | awk -F"* " '{print $2}'`
          curlabs=$curlabs""$newlab","
        else
          # Command lines are directly copied
          echo "$cmdline" >> $nfc                
        fi
      fi          
    done < "$cmddata"

    curlabs=`echo "$curlabs" | rev | cut -c 2- | rev`
    firstlab=`echo "$curlabs" | awk -F',' '{print $1}'`
    rangeout=$rangeout""$firstlab", "
    morelabs=`echo "$curlabs" | awk -F',' '{OFS=",";$1=""; print $0}'`
    morewc=`echo $morelabs | wc -w`      
    if [ $morewc != 0 ] ; then
      if [ -e "$labdata" ] ; then
        rm -f $labdata 
      fi
      echo "$morelabs" | awk -F',' '{for(i=1;i<=NF;i++) {print $i} }' > "$labdata"
      alias=" #ALIAS "
      while read thislab
      do
        labwc=`echo $thislab | wc -w`
        if [ $labwc != 0 ]; then
          alias=$alias"\"$thislab\" "
        fi
      done < "$labdata"
    else
      alias=""
    fi
    
    echo "RANGE "$rangeblue" "$rangered"   "$firstlab"   "$alias >> $nfc
    
    curstart=$nextstart
     
  done < "$rangedata"
  
  echo "" >> $nfc
  rangeout=`echo $rangeout | rev | cut -c 2- | rev`
  echo "$rangeout"
  
# remove dummy files
if [ -e $rangedata ] ; then
  rm -f $rangedata
fi
if [ -e $cmddata ] ; then
  rm -f $cmddata
fi
if [ -e $labdata ] ; then
  rm -f $labdata
fi
  
else
  echo "Error: FORMAL_CARDS not found"
fi

