#!/bin/bash
# creates NEWDATOM_INPUT based on current DATOM comment 
# Note: this is a simple implementation that works only for newer DATOM files
# Andreas Sander, last Update 15 Oct 2021

if [ -e DATOM ] ; then
  if [ -e NEWDATOM_INPUT ] ; then
    mv NEWDATOM_INPUT NEWDATOM_INPUT_OLD
  fi
  sed -n '/* NEWDATOM_INPUT was/,/^*=====/p' DATOM | sed 's/^\*//g' | sed -e 's/^ *//g' | tac | sed -e '1d' | tac | sed -e '1d' > NEWDATOM_INPUT
else
  echo "Error: DATOM not found"
fi
