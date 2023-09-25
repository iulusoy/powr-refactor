#!/bin/bash
# fully automatic update of levelcards after a new DATOM has
# been generated (meta-script)
# Andreas Sander, last update 15 Oct 2021

source ~/.powrconfig || exit

if [ -e CARDS ] ; then
  path=$POWR_WORK"/proc.dir"
  eval path=$path
  if [ ! -d $path ] ; then
    echo "Error: Script path $path could not be found!"
    exit
  fi

  levelcardsexe=$POWREXEPATH"/levelcards.exe.opt"
  eval levelcardsexe=$levelcardsexe
  if [ ! -e $levelcardsexe ]; then
    echo "Error: Program $levelcardsexe could not be found!"
    exit
  fi

  $path/rmlevels.bash
  $levelcardsexe
  $path/addlevels.bash
  echo 'AUTOLEVELS: CARDS have been updated'
else
  echo "Error: CARDS file could not be found"
fi
