#!/bin/bash
# Wrapper for displayatomicdata.bash
# Andreas Sander, last update 15 Oct 2021
# arguments: [-wrh | -htodt | -d <dir>] [<element>] 

source ~/.powrconfig || exit

script=$POWR_WORK/proc.dir/displayatomicdata.bash

if [ _"${POWR_INSTTYPE}" == _"potsdam" ] ; then
  pathhtodt="~htodt/datom/wrdata-archive/"
  pathwrh="~wrh/work/wrdata-archive/"
fi

# defaults
path=$pathwrh
filter=""
len=${#1}

if [[ "$1" == "-wrh" || "$1" == "-htodt" ]]; then
  if [ _"${POWR_INSTTYPE}" != _"potsdam" ]; then
    echo "Error: Options -wrh and -htodt are only possible inside the potsdam cluster."
    echo '       You can specify a custom repository with the -d option.'
    exit
  fi
  if [[ "$1" == "-wrh" ]]; then
    path=$pathwrh
  elif [[ "$1" == "-htodt" ]]; then
    path=$pathhtodt
  fi
  len=${#2}
  if [ $len -gt 0 ]; then
    filter=$2
  fi
elif [[ "$1" == '-d' ]]; then
  #show a different (custom) path?
  path=$2
  len=${#3}
  if [ $len -gt 0 ]; then
    filter=$3
  fi
elif [ $len -gt 0 ]; then
  filter=$1
fi

#expand tilde
eval script=$script
eval path=$path

if [ ${#filter} -gt 0 ]; then
  #ensure uppercase
  filter=`echo $filter | awk '{print toupper($0)}'`
  $script $path | grep -e "^$filter"_ -e '---'
  echo '### Results filtered for Element '$filter' ###'
else
  $script $path
fi
