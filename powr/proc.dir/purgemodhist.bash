#!/bin/bash
# Truncate MODHIST file to current part
# Andreas Sander

in=MODHIST
out=MODHIST_new.tmp

if [ -e $out ]; then
  rm -f $out
fi

readarray file < $in
lines=$(( ${#file[@]} -1 ))

# read old MODHIST file backwards : jobnum should be decreasing!
lastjob=-1
for (( line=$lines, i=$lines;((line >= 0 && i > 0)); line--, i-- )) ; do
  aline=${file[$line]}
  jobnum=$(echo "$aline" | sed -e "s:^/::" -e "s/\..*//")
  if (($lastjob > 0)) ; then
    if (($jobnum > $lastjob)) ; then
      # stop if current job is larger than old job 
      # since this lines (and onwards) must be from an older calculation
      break
    else
      echo -ne "$aline" >> $out
    fi
  else 
    echo -ne "$aline" >> $out
  fi
  lastjob=$jobnum
done

# revert output and store into old MODHIST
sed -n '1!G;h;$p' $out > $in
