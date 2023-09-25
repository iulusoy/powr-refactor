#!/bin/bash
# displays number of levels etc. for all DATOM.EL_ION and FORMAL_CARDS files in 
# a given (1st argument) directory
# Helge. 11.03.2016

source ~/.powrconfig || exit

if [ $# -eq 0 ]
then
 adir=$POWR_ATOM_DATA
else
 adir=$1
fi

echo "------------- Content of $adir -----------------"

# this aaaaa substitution trick is needed for correct sorting order of S an SI
for datomionpath in `ls $adir/DATOM.*_* | sed s/_/aaaaa/ | sort -t.  -k2`
do
 datomionpath=`echo $datomionpath | sed s/aaaaa/_/`
 datomion=${datomionpath##*.}
 NLEVEL=`grep "^LEVEL" $datomionpath | wc -l`
 NLINE=`grep "^LINE" $datomionpath | wc -l`
 clines=$(echo "$NLEVEL * ( $NLEVEL - 1 ) / 2 " | bc -l)
 testlines=$(echo "$clines == $NLINE" | bc -l)
 if [ $testlines -ne 1 ] ; then
     echo "error : $datomion has $NLINE LINES but $clines trans from $NLEVEL"
     exit
 fi    
 NCONT=`grep "^CONTINUUM" $datomionpath | wc -l`
 NKSHELL=`grep "^K-SHELL" $datomionpath | wc -l`
 NDRTRANS=`grep "^DRTRANSIT" $datomionpath | wc -l`
 formalcard=$adir/FORMAL_CARDS.$datomion
 NMULT=`grep "^+MULTIPLET" $formalcard | wc -l `
 NFORMAL=`grep "^UPPERLEVEL" $formalcard | wc -l `
 printf "%-7s : %3d LEVEL, %4d LINE, %3d CONTIUUM, %1d K-SHELL, %3d DR-TRANSIT, %4d MULTIPLET/LINE \n"  ${datomion} $NLEVEL $NLINE $NCONT $NKSHELL $NDRTRANS $NFORMAL
done

echo "------------- Content of $adir -----------------"
