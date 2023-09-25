#!/bin/bash
source ~/.powrconfig
IFS=''

help=false
if [ $# -eq 0 ]
then
  help=true
fi

if [ $# -eq 1 ]
then
  if [ $1 = 'help' ]
  then
    help=true
  fi
fi

if [ $help = 'true' ]
then
  echo 
  echo 
  echo 'How to use  msinfo :'
  echo '================================='
  echo ' '
  echo 'msinfo needs the name of the mass storgage file (fort.3, MODEL, FEDAT, ...)'
  echo 'and the command. The keyword "red" after the filname redirects output'
  echo 'to file /tmp/${USER}_msinfo.out'
  echo '                                                '
  echo "msinfo <file> [red] <command1> <command2> <...>"
  echo '                                                '
  echo 'The following commands are available :'
  echo "INFO      :  short info (default)"
  echo "INFO-L    :  long info"
  echo "INFO-D    :  prints the content of a variable"
  echo "               (additional parameters needed)"
  echo "  syntax  :  msinfo <file> INFO-D <name of the variable> <output format>"
  echo "             where output format is in Fortran style: Iw, Fw.d, Ew.dEe"
  echo 
  echo "example   : msinfo MODEL INFO red INFO-D 'ND' '(I6)'"
  echo 
  echo 
  exit
fi


model=$1
if [ ${model:0:1} == '/' ]
then
   mpath=""
elif [ ${model:0:1} == '~' ]
then
   mpath=""
else 
   mpath=$PWD
fi

cd /tmp || exit
ln -sf $mpath/$model /tmp/fort.33


if [ `echo $1 | cut -c1-5` = 'ASCII' ]
then
  cat $1 | grep 'NAME='
  exit
fi

rm -f msinfo.input
file=false
for arg
do
  arg2=`echo $arg | cut -c1-3`
  if [ $arg2 = 'red' ]
  then
    file=true
  else
    echo $arg >> msinfo.input
  fi
done

if [ $file = 'true' ]
then
  rm -f msinfo.out
  if [ ${HOSTTYPE%-linux} == 'x86_64' ] 
  then
      export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${POWR_WORK}/intellibs
      $POWR_WORK/exe.dir/msinfo.exe.opt > ${USER}_msinfo.out
  fi
  echo "Output written to /tmp/${USER}_msinfo.out"
else
  if [ ${HOSTTYPE%-linux} == 'x86_64' ] 
  then
      export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${POWR_WORK}/intellibs
      $POWR_WORK/exe.dir/msinfo.exe.opt
  fi
fi
rm -f fort.33
rm -f msinfo.input
