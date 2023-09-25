#!/bin/bash

source ~/.powrconfig || exit

# clones the wrdata directory of an existing chain to another one
# usage: clonechain sourcechain targetchain [user] [to-host | from-host]

filelist1="CARDS DATOM FGRID FORMAL_CARDS MODEL backup NEWDATOM_INPUT NEWFORMAL_CARDS_INPUT"
filelist2="FEDAT FEDAT_FORMAL"
othost=$HOST
othost2=$HOST
if [ $# -lt 2 ] # number of arguments < 2 => print out help
then 
 echo "You must give source and target chain numbers"
 echo 'Usage: clonechain source target [user] [to-host | from-host]'
 exit
elif [ $# -eq 3 ] # number of arguments = 3 => analyze 3rd argument
then
 otuser=$3 
 if [ ${otuser:0:3} == "to-" ]
 then
   othost=${otuser:3}
   otuser=$USER
 elif [ ${otuser:0:5} == "from-" ]
 then
   othost2=${otuser:5}
   otuser=$USER
 fi
elif [ $# -eq 4 ] # number of arguments = 4 => 3rd = USER, 4th = to- | from- host
then
 otuser=$3
 othost=$4
 if [ ${othost:0:3} == "to-" ]
 then
   othost=${othost:3}
 elif [ ${othost:0:5} == "from-" ]
 then
   othost2=${othost:5}
   othost=$HOST
 fi
fi

schain=$1
if [ "$schain" -eq "$schain" -a "$schain" -gt 0 ] 2>/dev/null 
then
    true
else
    echo "1st argument is not an integer number > 0, i.e. not a chain number: $schain"
    exit 
fi

tchain=$2
if [ "$tchain" -eq "$tchain" -a "$tchain" -gt 0 ] 2>/dev/null
then
     true
else
    echo "2nd argument is not an integer number > 0, i.e. not a chain number: $tchain"
    exit 
fi

if [ $schain -eq $tchain -a $# -eq 2 ]
then
    echo "Source chain $schain and target chain $tchain are the same on the local host"
    exit
fi

echo ""
echo "Dear user $USER on host ${HOST}, are you sure that you want to"
printf " overwrite your wrdata${tchain} on $othost with ${otuser:-$USER}'s wrdata${schain} from $othost2 ? "
printf "[yes/no]"
printf "\n"
read reply
if [ _$reply == '_yes' ]
then
 echo 'ok, overwriting'

 # copy wrdata
 if [ _$otuser == '_' ]
 then
  otuser=${USER}
 fi
 if [ $HOST == $othost -a $HOST == $othost2 ] # all on same host : NFS copy
 then 
  for file in $filelist1
  do
      eval cp -f  ~${otuser}/powr/wrdata${schain}/${file}                 $POWR_WORK/wrdata${tchain}/${file}
  done
  for file in $filelist2
  do  
      eval cp -Pf ~${otuser}/powr/wrdata${schain}/${file}                 $POWR_WORK/wrdata${tchain}/${file}
  done
  echo "Cloned from Chain ~${otuser}: wrdata${schain} on " `date` >> $POWR_WORK/wrdata${tchain}/wrdata.hist
 elif [ $HOST == $othost2 ] # source is current host ; target is remote host
 then
  for file in $filelist1
  do
      eval scp -Cp  ~${otuser}/powr/wrdata${schain}/${file}                 ${othost}:work/wrdata${tchain}/${file}
  done
  for file in $filelist2
  do
      eval rsync -avz ~${otuser}/powr/wrdata${schain}/${file}              ${othost}:work/wrdata${tchain}/${file}
  done
  ssh $othost "echo Cloned from Chain ${otuser}: wrdata$schain on `date` >> ~/powr/wrdata${tchain}/wrdata.hist"
 elif [ $HOST == $othost ] # source is remote host ; target is current host
 then
  for file in $filelist1
  do 
     eval scp -Cp  ${othost2}:~${otuser}/powr/wrdata${schain}/${file}                 $POWR_WORK/wrdata${tchain}/${file}
  done
  for file in $filelist2
  do
     eval rsync -avz ${othost2}:~${otuser}/powr/wrdata${schain}/${file}              $POWR_WORK/wrdata${tchain}/${file}
  done
  echo "Cloned from Chain ${otuser}: wrdata$schain on" `date` >> ~/powr/wrdata${tchain}/wrdata.hist
 fi
else
 echo 'abort'
fi

exit
