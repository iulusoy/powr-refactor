#!/bin/bash

if [ -e ~/.powrconfig ] ; then
     source ~/.powrconfig
fi     
#
# Prozedur status.com
# usage: status.com -u user
#        status.com help
#        status.com [t]cpr wruniq$kn
#   ... and many more (see below)

# functions:

getchainnumbers(){
user_temp_kall=$1
kall=`eval ls $POWR_WORKu | grep -E "wrdata[0-9]+$" \
 | cut -c7- | sort -n | tr '\n' ' '`
}

HOSTNAME=`hostname | cut -d '.' -f 1`
us=$USER
n=$#
i=1
while [ $i -lt $n ]
do 
  act=`eval echo "\$"$i`
  i=`expr $i + 1`
  if [ $act = -u ] ; then
    us=`eval echo "\$"$i`
    n=`expr $n - 2`
  fi
done

if [ _$us != _$USER ] ; then
    eval powrconfig=~${us}/.powrconfig
    if [ -e "$powrconfig" ] ; then
     alt_POWR_WORKu=`grep -m 1 POWR_WORK $powrconfig`
     # other user has ~/ -> replace by ~user/
     testtild=`echo "$alt_POWR_WORKu" | grep '~/'`
     if [ _"$testtild" != _  ] ; then
       alt_POWR_WORKu=`echo "$alt_POWR_WORKu" | sed s+~/+~${us}/+`
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
      POWR_WORKu=~$us/work/
    fi
    eval POWR_WORKu=$POWR_WORKu
else
    eval POWR_WORKu=${POWR_WORK:-~/work/}
fi      

#
if [ $n -ge 1 ]
then
if [ $1 = 'help' ] 
then
echo No parameter : info
echo
echo First parameter : 
echo           "         info     :   Information ueber den Inhalt der status-files"
echo           "         spec     :   Information ueber den Inhalt der status-files"
echo           "                        fuer spezielle Ketten"
echo           "         log      :   Ausgabe der log-files"
echo           "         cpr      :   Ausgabe der cpr-files"
echo           " tlog o. ltail    :   Staendige Ausgabe des log files"
echo           " tcpr o. ctail    :   Staendige Ausgabe des cpr files"
echo           "         wall     :   Ausgabe des gesamten cpr-files | grep Wall"
echo           "         purge    :   Kuerzen der .out, .cpr und .log files"
echo           "         nextjob  :   Ausgabe der Files nextjob"
echo           "         break    :   Abspeichern der Job-files bei naechster "
echo           "                        Gelegenheit, break->fbreak"
echo           "         no_break :   Kein Abspeichern der Job-files bei naechster "
echo           "                        Gelegenheit, no_break->fbreak"
echo           "         cload    :   Laedt das aktuelle CARDS file auf alle Maschinen"
echo           "         run      :   Kein Abspeichern der Job-files bei naechster "
echo           "                        Gelegenheit, run->fbreak, Job kann nicht! "
echo           "                        automatisch schlafen gelegt werden"
echo           "         stop     :   Stoppen des Jobs bei naechster "
echo           "                        Gelegenheit, stop->fbreak"
echo           "         KILL     :   kills the job instantly"
echo           "         sleep    :   Pausieren des Jobs bei naechster "
echo           "                        Gelegenheit, sleep->fbreak"
echo           "         swait    :   Status auf wait geschaltet"
echo           "         wait     :   Job schaltet sich bei der naechsten Aktion auf wait"
echo           "    "
echo           "    fast options  :"
echo           "        fwait     :   "
echo           "        fstop     :   "
echo           "        fbreak    :   "
echo           "    "
echo           "          nice    :   schreibt das Argument nach queue und queue2"
echo           "                             stat nice <NR> wruniq<nr>"
echo           "                        oder stat nice wruniq<nr> (Setzt auf den Default 19)"
echo           "    "
echo           "         non      :   schreibt non_active nach status"
echo           "         off      :   Schreibt das Wort off in das file status"
echo           "      to-machine  :   Break und starten auf neuer Maschine"
echo           "         ws       :   echo der ws.master files"
echo           "         name     :   Gibt der Kette einen Namen"
echo           "        +name     :   Aendert den Namen der Kette"
echo           "         name+    :   Aendert den Namen der Kette"
echo           "       kmname     :   Gibt der Kette einen besonderenNamen"
echo           "      +kmname     :   Aendert den besonderen Namen der Kette   (Not supported so far!)"
echo           "       kmname+    :   Aendert den besonderen Namen der Kette   (Not supported so far!)"
exit
fi
fi
#
# Kein Zweiter Parameter :      Frage nach dem File
#
# Zweiter Parameter      :      Name der Kette
#
# Kein Dritter Parameter :      Ausgabe der letzten 30 Zeilen
#                               Kuerzen bis auf die 1000 letzten Zeilen
#
# Dritter Parameter      :      Zahl der Zeilen oder all (->10000)
#
#
#op=tmp_data
###op=tmp_2day
#
echo
if [ $n -eq 0 ]
then
  action=info
  what=kall
else
  action=$1
  what='none'
fi

eval h="~${us}"
getchainnumbers ${us}

if [ $n -ge 1 ]
then
  c1=`echo $1 | cut -c1`
  c2=`echo $1 | cut -c2`
  if [ $1 = 'carina' -o \
       $1 = 'orion'  -o \
       $1 = 'hydra'  -o \
       $1 = 'gemini' -o \
       $1 = 'cygnus' -o \
       $1 = 'vela'   -o \
       $1 = 'leo'    -o \
       $1 = 'dorado' -o \
       $1 = 'corona' -o \
       $1 = 'all' ]
  then
    action=info
    what=$1
  elif [ $c1 = 'k' -a $c2 != 'm' ]
  then
    kall=`echo $1 | cut -c2-4`
    if [ $n -ge 2 ]
    then
     what=' '
     i=$kall
     while [ $i -lt $2 ]
     do 
       i=`expr $i + 1`
       kall=$kall' '$i
     done
    fi
    what=kall
    action=info
  elif [ $1 = 'spec' ]
  then
    what=spec
    action=noinfo
  fi
fi
#
if [ $action = 'info' ]
then
  echo "STATUS : User = $us"
  echo '====================='
  echo -e "\twhere\tWRstart\tWRuniq\tnext\tFormal\tName"
  echo -------------------------------------------------------------------------------
  for name in $kall
  do
    st1=`cat -s $POWR_WORKu/scratch/wrstart$name/status | cut -c 1-7 | sed s/active/ACTIVE/`
    fw2=`cat -s $POWR_WORKu/scratch/wruniq$name/fwhere | cut -c1-7`

    st2=`cat -s $POWR_WORKu/scratch/wruniq$name/status | cut -c 1-7 | sed s/active/ACTIVE/ | sed s/sleep/SLEEP/ | sed s/non/Non/`
    na=' '

    if test -e  $POWR_WORKu/scratch/wruniq$name/name
     then
     na=`cat -s $POWR_WORKu/scratch/wruniq$name/name`
    fi

    if test -e $POWR_WORKu/scratch/wruniq$name/fbreak
      then
      fb2=`cat -s $POWR_WORKu/scratch/wruniq$name/fbreak | cut -c1-6`
    else
      fb2='      '
    fi

    st4=`cat -s $POWR_WORKu/scratch/formal$name/status | cut -c 1-7 | sed s/active/ACTIVE/`
    if [ $fw2 = $what -o $what = 'all' -o $what = 'kall' -o $what = 'k2' ]
    then
##        echo -e "Ket.$name\t"$st1"\t\t"$fw2"\t"$fw22"\t"$st2"\t"$fb2"\t\t"$st3"\t"$st4"\t\t"$na
        echo -e "Ket.$name\t"$fw2"\t"$st1"\t"$st2"\t"$fb2"\t"$st4"\t"$na
    fi
done
exit
fi
#
if [ $what = 'spec' ]
then
  echo
  echo -n 'Sthlp    : '
  cat -s $POWR_WORKu/scratch/sthlp/status | sed s/^active/ACTIVE/
  echo -n 'Sthlp_kette: '
  cat -s $POWR_WORKu/scratch/sthlp_kette/status | sed s/^active/ACTIVE/
  echo
  echo -n 'Formal_MP: ' 
  cat -s $POWR_WORKu/scratch/formal_mp/status | sed s/^active/ACTIVE/
  echo
  st1=`cat -s $POWR_WORKu/scratch/fohlp/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st2=`cat -s $POWR_WORKu/scratch/fohlp2/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st3=`cat -s $POWR_WORKu/scratch/fohlp3/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st4=`cat -s $POWR_WORKu/scratch/fohlp4/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st5=`cat -s $POWR_WORKu/scratch/fohlp5/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st6=`cat -s $POWR_WORKu/scratch/fohlp6/status | cut -c 1-7 | sed s/active/ACTIVE/`
  echo -e "Fohlp 1-6:\t"$st1"\t"$st2"\t"$st3"\t"$st4"\t"$st5"\t"$st6
  echo
  st1=`cat -s $POWR_WORKu/scratch/mmformal/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st2=`cat -s $POWR_WORKu/scratch/mmformal2/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st3=`cat -s $POWR_WORKu/scratch/mmformal3/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st4=`cat -s $POWR_WORKu/scratch/mmformal4/status | cut -c 1-7 | sed s/active/ACTIVE/`
  st5=`cat -s $POWR_WORKu/scratch/mmformal5/status | cut -c 1-7 | sed s/active/ACTIVE/`
  echo -e "MMFormal 1-5:\t"$st1"\t"$st2"\t"$st3"\t"$st4"\t"$st5
#  st1=`cat -s $POWR_WORKu/scratch/mmformal5/status | cut -c 1-7 | sed s/active/ACTIVE/`
#  st2=''
#  st3=''
#  st4=''
#  echo -e "MMformal 5  :\t"$st1"\t\t"$st2"\t\t"$st3"\t\t"$st4
  echo
  echo -n 'FormalRange: ' 
  cat -s $POWR_WORKu/scratch/formalrange/status | sed s/^active/ACTIVE/
  echo
  echo -n 'Formal_kette1: ' 
  cat -s $POWR_WORKu/scratch/formal_kette1/status | sed s/^active/ACTIVE/
  echo -n 'Formal_kette2: ' 
  cat -s $POWR_WORKu/scratch/formal_kette2/status | sed s/^active/ACTIVE/
  echo -n 'Formal_kette3: ' 
  cat -s $POWR_WORKu/scratch/formal_kette3/status | sed s/^active/ACTIVE/
  echo -n 'Formal_kette4: ' 
  cat -s $POWR_WORKu/scratch/formal_kette4/status | sed s/^active/ACTIVE/
  echo -n 'Formal_kette5: ' 
  cat -s $POWR_WORKu/scratch/formal_kette5/status | sed s/^active/ACTIVE/
  echo -n 'Formal_kette6: ' 
  cat -s $POWR_WORKu/scratch/formal_kette6/status | sed s/^active/ACTIVE/
  echo -n 'Formal_kette7: ' 
  cat -s $POWR_WORKu/scratch/formal_kette7/status | sed s/^active/ACTIVE/
  echo -n 'Formal_kette8: ' 
  cat -s $POWR_WORKu/scratch/formal_kette8/status | sed s/^active/ACTIVE/
  exit
fi
#
if [ $us != $USER ]
then
  echo User was $us
  echo ==============
fi

if [ $action = 'log' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which file do you want?'
    exit
  else
    log=$2
    if [ $n -eq 2 ]
    then
      vtail=30
    else
     if [ $3 = 'all' ]
     then
       vtail=10000
     else
       vtail=$3
     fi
    fi
  fi
  tail -$vtail $POWR_WORKu/output/$log.log
  exit
fi
#
if [ $action = 'ltail' -o $action = 'tlog' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which file do you want?'
    exit
  else
    log=$2
  fi
  tail -f $POWR_WORKu/output/$log.log
  exit
fi
#
if [ $action = 'cpr' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which file do you want?'
    exit
  else
    cpr=$2
    if [ $n -eq 2 ]
    then
      vtail=30
    else
     if [ $3 = 'all' ]
     then
       vtail=10000
     else
       vtail=$3
     fi
    fi
  fi
  if [ `echo $cpr | cut -c1-6` = 'wruniq' ]
  then
  tail -$vtail $POWR_WORKu/output/$cpr.cpr
  else
    tail -$vtail $POWR_WORKu/output/$cpr.cpr
  fi
  exit
fi
#
if [ $action = 'ctail' -o $action = 'tcpr' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which file do you want?'
    exit
  else
    cpr=$2
  fi
    tail -f $POWR_WORKu/output/$cpr.cpr
  exit
fi
#
if [ $action = 'wall' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which file do you want?'
    exit
  else
    cpr=$2
    if [ $n -eq 2 ]
    then
      vtail=30
    else
      if [ $3 = 'all' ]
      then
        vtail=100000
      else
        vtail=$3
      fi
    fi
  fi
  cat -s $POWR_WORKu/output/$cpr.cpr | grep Wall | tail -n $vtail
  exit
fi
#
#
if [ $action = 'out' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which file do you want?'
    exit
  else
    out=$2
    if [ $n -eq 2 ]
    then
      vtail=30
    else
     if [ $3 = 'all' ]
     then
       vtail=10000
     else
       vtail=$3
     fi
    fi
  fi
  cat -s $POWR_WORKu/output/$out.out | tail -n $vtail
  exit
fi

if [ $action = 'name' -o $action = '+name' -o $action = 'name+' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to name?'
    exit
  else
    k=$2
  fi
  if [ $n -eq 2 ]
  then
    text=' '
  else
    text=$3
  fi
  if [ $action = 'name' ]
  then
    echo $text > $POWR_WORKu/scratch/$k/name
  elif [ $action = 'name+' ]
  then
    echo $text >> $POWR_WORKu/scratch/$k/name
  elif [ $action = '+name' ]
  then
    cp $POWR_WORKu/scratch/$k/name $POWR_WORKu/scratch/$k/name.scratch
    echo $text > $POWR_WORKu/scratch/$k/name
    cat $POWR_WORKu/scratch/$k/name.scratch >> $POWR_WORKu/scratch/$k/name
    rm $POWR_WORKu/scratch/$k/name.scratch
  fi
  exit
fi

if [ $action = 'kmname' -o $action = '+kmname' -o $action = 'kmname+' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to name?'
    exit
  else
    k=$2
  fi
  if [ $n -eq 2 ]
  then
    text='_'
  else
    text=$3
  fi
  if [ $action = 'kmname' ]
  then
    echo $text > $POWR_WORKu/scratch/$k/kmname
  elif [ $action = 'kmname+' ]
  then
#    echo $text >> $POWR_WORuK/scratch/$k/kmname
    echo action kmname+ not supported so far
  elif [ $action = '+kmname' ]
  then
#    cp $POWR_WORKu/scratch/$k/name $POWR_WORKu/scratch/$k/kmname.scratch
#    echo $text > $POWR_WORKu/scratch/$k/kmname
#    cat $POWR_WORKu/scratch/$k/kmname.scratch >> $POWR_WORKu/scratch/$k/kmname
#    rm $POWR_WORKu/scratch/$k/kmname.scratch
    echo action +kmname not supported so far
  fi
  exit
fi

if [ $action = 'nice' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to name?'
    exit
  else
    k=$3
  fi
  if [ $n -eq 2 ]
  then
    text='19'
    k=$2
  else
    text=$2
  fi
  echo $text > $POWR_WORKu/scratch/$k/queue
  echo $text > $POWR_WORKu/scratch/$k/queue2
  exit
fi

if [ $action = 'KILL' ]
then
 jobname=`echo $2 | sed s/[0-9]//g`
 jobnumber=`echo $2 | sed s/[a-z]//g`
 scratchwhere=$jobname
 killprog="" # kill adapter, como, etc.
 if [ _$jobname = _formal ]
 then
    scratchwhere=wruniq
    killprog="formal" # kill specifically formal
 fi
 if [ _$jobname = _spawnformals ]
 then
    scratchwhere=wruniq
    killprog="spawnformals" # kill specifically formal
    echo 'stop' > $HOME/powr/scratch/wruniq${jobnumber}/fbreak
 fi
 scratchwhere=${scratchwhere}${jobnumber}
 killhost=`cat $HOME/powr/scratch/${scratchwhere}/fwhere`
 echo "KILL to $2 on $killhost $killprog"
 ssh $killhost "kill -KILL \`ps axu | grep $USER | grep ${jobnumber}_${killprog}  | grep -v 'grep' | head -1 | awk '{print \$2}'\`"
 exit
fi

if [ $action = 'break' -o \
     $action = 'no_break' -o \
     $action = 'srun' -o \
     $action = 'run' -o \
     $action = 'sleep' -o \
     $action = 'stop' -o \
     $action = 'wait' -o \
     $action = 'fwait' -o \
     $action = 'fstop' -o \
     $action = 'fbreak' -o \
     `echo $action | cut -c1-3` = 'to-' ]
then
  if [ $action = 'srun' ]
  then
    action=no_break
  fi
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to modify?'
    exit
  else
    k=$2
    kn=`echo $2 | cut -c7-10`
  fi

  echo $action > $POWR_WORKu/scratch/$k/fbreak

  exit
fi

if [ $action = 'cload' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to cload?'
    exit
  else
    k=$2
  fi
  kn=`echo $k | cut -c7-10`
  fw=`cat $POWR_WORKu/scratch/$k/fwhere`
  st=`cat $HOME/powr/scratch/$k/status`
  if [ $st != 'sleep' -a $st != 'active' ]
  then
    echo cload nur bei laufenden Ketten!
    echo stat cload was abortet
    exit
  fi
  cp $HOME/powr/wrdata$kn/CARDS /home/$fw/tmp_data/$us/$k/CARDS
  exit
fi

#Only for Jobs running on external Powrstations
if [ $action = 'resub' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to modify?'
    exit
  else
    k=$2
  fi
  echo 'Do resub only on external WS!'
  echo 'Command ignored'
  exit
fi

if [ $action = 'status' ]
then
  if [ $n -lt 3 ]
  then
    echo 'Bedienung von stat status incorrect'
    exit
  else
    k=$3
    action=$2
  fi
  echo $action > $POWR_WORKu/scratch/$k/status
  exit
fi

if [ $action = 'queue' ]
then
  if [ $n -lt 3 ]
  then
    echo 'Bedienung von stat status incorrect'
    exit
  else
    k=$3
    action=$2
  fi
  echo $action > $POWR_WORKu/scratch/$k/fwhere2
  exit
fi

if [ $action = 'clear' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which status do you want to clear?'
    exit
  else
    clear=$2
  fi
  echo 'clear' > $POWR_WORKu/scratch/$clear/status
  exit
fi
#
if [ $action = 'abort' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which status do you want to abort?'
    exit
  else
    abort=$2
  fi
  echo 'abort' > $POWR_WORKu/scratch/$abort/status
  exit
fi
#
#
if [ $action = 'osleep' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to osleep?'
    exit
  else
    abort=$2
  fi
  echo 'osleep' > $POWR_WORKu/scratch/$abort/fbreak
  exit
fi
#
if [ $action = 'swait' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to swait?'
    exit
  else
    abort=$2
  fi
  echo 'wait' > $POWR_WORKu/scratch/$abort/status
  exit
fi
if [ $action = 'non' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to non_active?'
    exit
  else
    abort=$2
  fi
  echo 'non_active' > $POWR_WORKu/scratch/$abort/status
  exit
fi
#
if [ $action = 'purge' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which file do you want?'
    exit
  else
    purge=$2
    if [ $n -eq 2 ]
    then
      vtail=1000
    else
       vtail=$3
    fi
  fi
  rm $POWR_WORKu/output/status.scratch
  mv $POWR_WORKu/output/$purge.out $POWR_WORKu/output/status.scratch
  tail -n $vtail $POWR_WORKu/output/status.scratch > $POWR_WORKu/output/$purge.out
#
  rm $POWR_WORKu/output/status.scratch
  mv $POWR_WORKu/output/$purge.log $POWR_WORKu/output/status.scratch
  tail -n $vtail $POWR_WORKu/output/status.scratch > $POWR_WORKu/output/$purge.log
#
  rm $POWR_WORKu/output/status.scratch
  mv $POWR_WORKu/output/$purge.cpr $POWR_WORKu/output/status.scratch
  tail -n $vtail $POWR_WORKu/output/status.scratch > $POWR_WORKu/output/$purge.cpr
  exit
fi
#
if [ $action = 'off' ]
then
  if [ $n -eq 1 ]
  then
    echo 'Which Job do you want to off?'
    exit
  else
    abort=$2
  fi
  echo 'off' > $POWR_WORKu/scratch/$abort/status
  exit
fi
#
echo Unknown action : $action
