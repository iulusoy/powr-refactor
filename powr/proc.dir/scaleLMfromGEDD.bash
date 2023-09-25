#!/bin/bash
# Scale L for a given M (or M for a given L) by assuming constant GEDD
# from file modinfo.kasdefs
# arguments: [L|l|M|m] [<log(L/L_sun)>|<M/M_sun>]

file='modinfo.kasdefs'
gedd=`grep 'GEDD ' $file | awk -F'=' '{print $2}'`
mstar=`grep 'MSTAR ' $file | awk -F'=' '{print $2}'`
logl=`grep 'LOGL ' $file | awk -F'=' '{print $2}'`
lum=`perl -le "printf(\"%.2f\", ( 10**$logl ) );"`
vinf=`grep 'VFINAL ' $file | awk -F'=' '{print $2}'`
c='29979245800'

xlsun='3.85E33'
xmsunyr='6.303E25'

# Calculate Mass
# identifier can be uppercase oder lowercase
if [[ $1 == 'L' || $1 == 'l' ]]; then
  echo "For a luminosity of log L/Lsun = $2, the scaled mass for same Gamma_Edd is M/Msun = "
  perl -le "printf(\"%.2f\", ( 10**$2 / $lum * $mstar ) );"
elif [[ $1 == 'M' || $1 == 'm' ]]; then
  echo "For a mass of $2 Msun, the scaled luminosity for same Gamma_Edd is log L/Lsun = "
  perl -le "printf(\"%.2f\", ( log( $lum * $2 / $mstar ) / log(10)  ) );"
else 
  echo "Error: Unknown identfier or no argument given!"
fi

echo ""

