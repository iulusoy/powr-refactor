#!/bin/bash
# Mini script for easier call of msinfo with INFO-D 

source ~/.powrconfig

var="$1"
modfile="$2"
if [ ${#var} -eq 0 ] ; then
  echo 'Syntax: msread <variable> [<ms-file>]'
  exit
fi

if [ ${#modfile} -eq 0 ] ; then 
# no file given: Select standard MODEL file
  if [ -e MODEL ] ; then
    modfile='MODEL'
  elif [ -e fort.3 ] ; then
    modfile='fort.3'
  else
    echo "msread error: Neither MODEL nor fort.3 could be found!"
    exit
  fi 
fi
echo ' * msread: Reading from mass-storage file '$modfile' ...'
format=auto

# Manual check for often used values that do not match auto format
charvars='MODHEAD GEFFKEY INCRIT VCRIT MODEL LEVEL'
for curvar in $charvars ; do
  if [[ $var == $curvar ]] ; then
    format='(A8)'
  fi
done

floatvars='JTOTL KTOTL NTOTL'
for curvar in $floatvars ; do
  if [[ $var == $curvar ]] ; then
    format='(G20.10)'
  fi
done

# msinfo call
$POWR_WORK/proc.dir/msinfo.com $modfile INFO-D "${var}" "${format}"
