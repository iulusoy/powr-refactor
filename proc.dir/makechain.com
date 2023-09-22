#!/bin/bash

source ~/.powrconfig || exit

# this script generates a new chain for PoWR with number "kn" given via call
# created by H. Todt 05.03.2008
# usage: makechain -wrh -sforce 32

export workdir="$POWR_WORK"

if [ $# -lt 1 ] ; then
 echo "No argument given: stop -> print help" 
 echo 'Usage: makechain.bash [-wrh|-goetz] [-sforce|-force] $kn'
 echo 'Argument :'
 echo ' $kn              - chain number (integer >= 1)'
 echo 'Options  :'
 echo ' -wrh             - chain for wrh-PoWR version (default)'
 echo ' -hydro           - chain for hydro-PoWR version'
 echo ' -force           - overwrites wrdata$kn and scratch directories and scripts'
 echo ' -sforce          - overwrites only scripts' 
 echo ' -all             - perform given action on all existing(!) chains'
 echo ' -only$kn1-$kn2   - perform given action only on the given existing(!) chains'
 echo 'Examples :'
 echo ' makechain.bash -wrh -sforce 32'
 exit
fi

# test existence of directories:
if [ ! -d $POWR_WORK/wrjobs ]  ; then
 echo "Creating: $POWR_WORK/wrjobs"
 mkdir "$POWR_WORK/wrjobs"
fi

if [ ! -d $POWR_WORK/scratch ] ; then
 echo "Creating: $POWR_WORK/scratch"
 mkdir "$POWR_WORK/scratch"
fi

if [ ! -d $POWR_WORK/output ] ; then 
 echo "Creating: $POWR_WORK/output"
 mkdir "$POWR_WORK/output"
fi

# test existence of wrjobs/tmphosts, otherwise create it:
if [ ! -s $POWR_WORK/wrjobs/tmphosts ] ; then 
 echo "Could not find a VALID $POWR_WORK/wrjobs/tmphosts file: Creation enforced."
 homehost=`echo $HOME | cut -d/ -f3`
 echo "default $homehost" > $POWR_WORK/wrjobs/tmphosts
fi

# default for PoWR version:
powr='wrh'

forced='false'
fwrdata='false'
kn='void'

chainnumbers='one'
# possible values: 
# all
# $kn1-$kn2
firstchain='0'
lastchain='0'

let iforce=0
let ipowr=0
let ihydro=0

# scan arguments and options
for arg in $@  ; do 
 if [ `echo $arg | cut -c1` = - ] ; then # scan options
  if [ $arg == '-force' ] ; then         # overwrite existing chain
   forced='true'
   fwrdata='true'
   let "iforce = iforce + 1"
  elif [ $arg == '-sforce' ] ; then     # overwrite only existing scripts
   forced='true'
   fwrdata='false'
   let "iforce = iforce + 1"
  elif [ $arg == '-wrh' ] ; then        # use wrh PoWR
   powr='wrh'
   echo "Using wrh-PoWR version!"
   let "ipowr = ipowr + 1"
  elif [ $arg == '-hydro' ] ; then      # use hydro PoWR version
   powr='wrh'
   echo "Using wrh-PoWR version with hydrodynamic branch!"
   let "ipowr = ipowr + 1"
   let "ihydro = ihydro + 1"      
  elif [ $arg == '-all' ] ; then
    # update all chains
    echo "Renew all existing chains"
    chainnumbers='all'
 elif [ ${arg:0:5} == '-only' ] ; then 
   # update only chains in given range
   echo "Renew only existing chains in selected range:"
   chainnumbers=${arg:5}
   firstchain=${chainnumbers%-*}
   lastchain=${chainnumbers#*-}
   if [ ${firstchain} -gt ${lastchain} ] ; then
     tempchain=$firstchain
     firstchain=$lastchain
     lastchain=$tempchain
   fi
   chainnumbers='only'
 else  #unknown option
   echo "Unknown option: $arg stop"
   exit
  fi  
 else                             # scan for other arguments
  if [ $kn == 'void'  ] ;  then   # number uninitialized
   if [[ ! $arg -ge 1 ]] ; then   # argument is not a number 
    echo "One Argument must be a number >= 1, you gave: $arg" 
    exit
   else                           # kn is set to argument
    kn=$arg
   fi
  else                            # only one argument allowed
   echo "Sorry only one chain per explicit call"
   exit
  fi                       
 fi
done

# echo $chainnumbers $kn $firstchain $lastchain

# get all existing wrdata numbers:
if [ "${chainnumbers}" == 'all' ]  ; then
  thenumbers=`ls $POWR_WORK |grep wrdata | cut -c7- |grep "[0-9]" |sort -n`
 elif [ "${chainnumbers}" == 'only' ] ; then
   thenumbers=`ls $POWR_WORK |grep wrdata | cut -c7- |grep "[0-9]" |sort -n`
   for ifound in $thenumbers ; do
    if [ ${ifound} -ge ${firstchain} -a ${ifound} -le ${lastchain} ] ; then
      mynumbers="${mynumbers} $ifound"
    fi
   done 
   thenumbers=$mynumbers
 else
  thenumbers=$kn
fi


# check for inconsistent arguments
if [[ $iforce -gt 1 ]] ; then
 echo "More than one force option given: stop"
 exit
fi
if [[ $ipowr -gt 1 ]] ; then
 echo "More than one PoWR version option given: stop"
 exit
fi

if [ "${chainnumbers}" == 'one' -a "${kn}" == 'void' ] ; then
 echo "Your forgot the chain number: stop"
 exit
fi

echo "PoWR version is: " $powr

#==================================================================
#
# Loop over all entries in thenumbers - one or more chain numbers  
#
#==================================================================
for kn in $thenumbers ; do
 if [ -d $POWR_WORK/wrdata$kn -a $forced == 'false' ]  ;then
    echo "Chain $kn may already exist - therefore not generated. Try \"-[s]force\"."
    exit 1
 elif [ -d  $POWR_WORK/wrdata$kn -a $forced == 'true' ] ; then
    if [ $fwrdata == 'true' ] ; then
     echo "Creation of chain $kn forced"
    else
     echo "Update of script $kn forced"
    fi
 else
    echo "Chain $kn does not not exist yet, try to create."
    fwrdata='true'
 fi

#######################################################################
# Creation of Chain directory and necessary Data from $POWR_WORK/dummychain
#######################################################################

if [ $fwrdata == 'true' -o  ! -d $POWR_WORK/wrdata$kn ] ; then 
 mkdir -p $POWR_WORK/wrdata$kn
 cd $POWR_WORK/wrdata$kn || exit 
 cp $POWR_WORK/dummychain/DATOM                 $POWR_WORK/wrdata$kn/DATOM
 cp $POWR_WORK/dummychain/FGRID                 $POWR_WORK/wrdata$kn/FGRID
 cp $POWR_WORK/dummychain/CARDS_$powr           $POWR_WORK/wrdata$kn/CARDS
 cp $POWR_WORK/dummychain/FORMAL_CARDS          $POWR_WORK/wrdata$kn/FORMAL_CARDS
 cp $POWR_WORK/dummychain/NEWDATOM_INPUT        $POWR_WORK/wrdata$kn/NEWDATOM_INPUT
 cp $POWR_WORK/dummychain/NEWFORMAL_CARDS_INPUT $POWR_WORK/wrdata$kn/NEWFORMAL_CARDS_INPUT
 if [ -e $POWR_WORK/dummychain/MODEL ] ; then
  cp $POWR_WORK/dummychain/MODEL                $POWR_WORK/wrdata$kn/MODEL
 fi
 if [ -e ${POWR_ATOM_DATA}/G2_BIG_VD100_FeXVI-K2015-parity.mS ] ; then
  ln -sf ${POWR_ATOM_DATA}/G2_BIG_VD100_FeXVI-K2015-parity.mS $POWR_WORK/wrdata$kn/FEDAT 
 else
  echo "WARNING ${POWR_ATOM_DATA}/G2_BIG_VD100_FeXVI-K2015-parity.mS DOES NOT EXIST:"
  echo "FEDAT is not created, please add later."
 fi
 if [ -e ${POWR_ATOM_DATA}/G2_BIG_VD100_FeXVI-K2015-parity.mS ] ; then
  ln -sf ${POWR_ATOM_DATA}/G2_BIG_VD100_FeXVI-K2015-parity.mS $POWR_WORK/wrdata$kn/FEDAT_FORMAL
 else
  echo "WARNING ${POWR_ATOM_DATA}/G2_BIG_VD100_FeXVI-K2015-parity.mS DOES NOT EXIST:"
  echo "FEDAT_FORMAL is not created, please add later."
 fi  
 echo "wrdata$kn created"
else
 echo "No wrdata$kn created - just overwriting old scripts"
fi
 
##############################################################
# Creation of scratch directories and initial Data for "stat"
##############################################################

if [ $fwrdata == 'true' -o  ! -d $POWR_WORK/wrdata$kn ] ; then 

 # scratch wrstart
 mkdir -p            $POWR_WORK/scratch/wrstart$kn
 touch               $POWR_WORK/scratch/wrstart$kn/fwhere
 echo "nohost"     > $POWR_WORK/scratch/wrstart$kn/fwhere
 touch               $POWR_WORK/scratch/wrstart$kn/status
 echo "non_active" > $POWR_WORK/scratch/wrstart$kn/status

 # scratch wruniq
 mkdir -p            $POWR_WORK/scratch/wruniq$kn
 touch               $POWR_WORK/scratch/wruniq$kn/fwhere
 echo "nohost"     > $POWR_WORK/scratch/wruniq$kn/fwhere
 touch               $POWR_WORK/scratch/wruniq$kn/status
 echo "non_active" > $POWR_WORK/scratch/wruniq$kn/status
 touch               $POWR_WORK/scratch/wruniq$kn/fbreak
 echo "no_break"   > $POWR_WORK/scratch/wruniq$kn/fbreak
 touch               $POWR_WORK/scratch/wruniq$kn/queue
 echo "nohost"     > $POWR_WORK/scratch/wruniq$kn/queue

 # scratch miscellaneous
 mkdir -p            $POWR_WORK/scratch/newdatom$kn
 mkdir -p            $POWR_WORK/scratch/newformal_cards$kn
 mkdir -p            $POWR_WORK/scratch/njn$kn
 mkdir -p            $POWR_WORK/scratch/modify$kn
 mkdir -p            $POWR_WORK/scratch/steal$kn
 touch               $POWR_WORK/scratch/steal$kn/status
 echo "non_active" > $POWR_WORK/scratch/steal$kn/status

 # scratch formal
 mkdir -p            $POWR_WORK/scratch/formal$kn
 touch               $POWR_WORK/scratch/formal$kn/status
 echo "non_active" > $POWR_WORK/scratch/formal$kn/status
 echo "scratch subdirectories for chain $kn created"

else

 echo "No scratch$kn created - just overwriting old scripts"

fi

#############################################################
# Creation of wrjobs
#############################################################

cd $POWR_WORK/dummychain || exit

# newdatom ----------------------------------------------------

touch $POWR_WORK/wrjobs/newdatom$kn
echo '#!/bin/bash'                                    > $POWR_WORK/wrjobs/newdatom$kn
echo 'source ~/.powrconfig'                          >> $POWR_WORK/wrjobs/newdatom$kn
echo '. $POWR_WORK/wrjobs/newdatom_gen' "$kn" '$1'   >> $POWR_WORK/wrjobs/newdatom$kn
chmod 754                                               $POWR_WORK/wrjobs/newdatom$kn

# if [ ! -e $POWR_WORK/wrjobs/wrstart_${powr}_gen ] ; then
 cp  $POWR_WORK/dummychain/newdatom_gen                 $POWR_WORK/wrjobs/newdatom_gen
 echo "Copied newdatom_gen to wrjobs"
# fi

if [ ! -e $POWR_WORK/wrdata${kn}/NEWDATOM_INPUT ] ; then
  cp $POWR_WORK/dummychain/NEWDATOM_INPUT               $POWR_WORK/wrdata${kn}/NEWDATOM_INPUT
fi 

# newformal_cards ----------------------------------------------------

touch                                                        $POWR_WORK/wrjobs/newformal_cards$kn
echo '#!/bin/bash'                                         > $POWR_WORK/wrjobs/newformal_cards$kn
echo 'source ~/.powrconfig'                               >> $POWR_WORK/wrjobs/newformal_cards$kn
echo '. $POWR_WORK/wrjobs/newformal_cards_gen' "$kn" '$1' >> $POWR_WORK/wrjobs/newformal_cards$kn
chmod 754                                                    $POWR_WORK/wrjobs/newformal_cards$kn

# if [ ! -e $POWR_WORK/wrjobs/wrstart_${powr}_gen ] ; then
 cp  $POWR_WORK/dummychain/newformal_cards_gen               $POWR_WORK/wrjobs/newformal_cards_gen
 echo "Copied newformal_cards_gen to wrjobs"
# fi

if [ ! -e $POWR_WORK/wrdata${kn}/NEWFORMAL_CARDS_INPUT ] ; then
  cp $POWR_WORK/dummychain/NEWFORMAL_CARDS_INPUT             $POWR_WORK/wrdata${kn}/NEWFORMAL_CARDS_INPUT
fi 


# wrstart ----------------------------------------------------

touch                                                           $POWR_WORK/wrjobs/wrstart$kn
echo '#!/bin/bash'                                            > $POWR_WORK/wrjobs/wrstart$kn
echo 'source ~/.powrconfig'                                  >> $POWR_WORK/wrjobs/wrstart$kn
# echo 'echo $HOSTNAME > $POWR_WORK/scratch/'"wruniq$kn/fwhere"  >> $POWR_WORK/wrjobs/wrstart$kn
if [[ $ihydro -gt 0 ]]; then
 echo '. $POWR_WORK/wrjobs/wrstart_'${powr}_hydro "$kn" '$1' >> $POWR_WORK/wrjobs/wrstart$kn
 echo '# . $POWR_WORK/wrjobs/wrstart_'${powr}_gen "$kn" '$1' >> $POWR_WORK/wrjobs/wrstart$kn
else
  echo '. $POWR_WORK/wrjobs/wrstart_'${powr}_gen "$kn" '$1'  >> $POWR_WORK/wrjobs/wrstart$kn
fi
echo '# . $POWR_WORK/wrjobs/wrstart_'${powr}_vd20 "$kn" '$1' >> $POWR_WORK/wrjobs/wrstart$kn
echo '# . $POWR_WORK/wrjobs/wrstart_'${powr}_xxl "$kn" '$1'  >> $POWR_WORK/wrjobs/wrstart$kn
chmod 754 $POWR_WORK/wrjobs/wrstart$kn

# if [ ! -e $POWR_WORK/wrjobs/wrstart_${powr}_gen ] ; then
 cp  $POWR_WORK/dummychain/wrstart_${powr}_gen                  $POWR_WORK/wrjobs/wrstart_${powr}_gen
 cp  $POWR_WORK/dummychain/wrstart_${powr}_vd20                 $POWR_WORK/wrjobs/wrstart_${powr}_vd20
 cp  $POWR_WORK/dummychain/wrstart_${powr}_xxl                  $POWR_WORK/wrjobs/wrstart_${powr}_xxl
if [[ $ihydro -gt 0 ]]; then
 cp  $POWR_WORK/dummychain/wrstart_${powr}_hydro                $POWR_WORK/wrjobs/wrstart_${powr}_hydro
fi
 echo "Copied wrstart_${powr}_gen, wrstart_${powr}_vd20 wrstart_${powr}_xxl to wrjobs"
# fi

# wruniq ----------------------------------------------------------

touch                                                           $POWR_WORK/wrjobs/wruniq$kn
echo '#!/bin/bash'                                            > $POWR_WORK/wrjobs/wruniq$kn
echo 'source ~/.powrconfig'                                  >> $POWR_WORK/wrjobs/wruniq$kn
if [[ $ihydro -gt 0 ]]; then
  echo '. $POWR_WORK/wrjobs/wruniq_'${powr}_hydro "$kn" '$1' >> $POWR_WORK/wrjobs/wruniq$kn
  echo '# . $POWR_WORK/wrjobs/wruniq_'${powr}_gen "$kn" '$1' >> $POWR_WORK/wrjobs/wruniq$kn
else
  echo '. $POWR_WORK/wrjobs/wruniq_'${powr}_gen "$kn" '$1'   >> $POWR_WORK/wrjobs/wruniq$kn
fi
echo '# . $POWR_WORK/wrjobs/wruniq_'${powr}_vd20 "$kn" '$1'  >> $POWR_WORK/wrjobs/wruniq$kn
echo '# . $POWR_WORK/wrjobs/wruniq_'${powr}_xxl "$kn" '$1'   >> $POWR_WORK/wrjobs/wruniq$kn
chmod 754                                                       $POWR_WORK/wrjobs/wruniq$kn

# if [ ! -e $POWR_WORK/wrjobs/wruniq_${powr}_gen ] ; then
 cp  $POWR_WORK/dummychain/wruniq_${powr}_gen                   $POWR_WORK/wrjobs/wruniq_${powr}_gen
 cp  $POWR_WORK/dummychain/wruniq_${powr}_vd20                  $POWR_WORK/wrjobs/wruniq_${powr}_vd20
 cp  $POWR_WORK/dummychain/wruniq_${powr}_xxl                   $POWR_WORK/wrjobs/wruniq_${powr}_xxl
if [[ $ihydro -gt 0 ]]; then
  cp  $POWR_WORK/dummychain/wruniq_${powr}_hydro                $POWR_WORK/wrjobs/wruniq_${powr}_hydro
fi
 echo "Copied wruniq_${powr}_gen, wruniq_${powr}_vd20 and wruniq_${powr}_xxl to wrjobs"
# fi

# formal --------------------------------------------------------------

touch $POWR_WORK/wrjobs/formal$kn
echo '#!/bin/bash'                                            > $POWR_WORK/wrjobs/formal$kn
echo 'source ~/.powrconfig'                                  >> $POWR_WORK/wrjobs/formal$kn
if [[ $ihydro -gt 0 ]]; then
  echo '. $POWR_WORK/wrjobs/formal_'${powr}_hydro "$kn" '$1' >> $POWR_WORK/wrjobs/formal$kn
  echo '# . $POWR_WORK/wrjobs/formal_'${powr}_gen "$kn" '$1' >> $POWR_WORK/wrjobs/formal$kn
else
  echo '. $POWR_WORK/wrjobs/formal_'${powr}_gen "$kn" '$1'   >> $POWR_WORK/wrjobs/formal$kn
fi 
echo '# . $POWR_WORK/wrjobs/formal_'${powr}_vd20 "$kn" '$1'  >> $POWR_WORK/wrjobs/formal$kn
echo '# . $POWR_WORK/wrjobs/formal_'${powr}_xxl "$kn" '$1'   >> $POWR_WORK/wrjobs/formal$kn
chmod 754 $POWR_WORK/wrjobs/formal$kn

# if [ ! -e $POWR_WORK/wrjobs/wruniq_${powr}_gen ] ; then
 cp  $POWR_WORK/dummychain/formal_${powr}_gen                  $POWR_WORK/wrjobs/formal_${powr}_gen
 cp  $POWR_WORK/dummychain/formal_${powr}_vd20                 $POWR_WORK/wrjobs/formal_${powr}_vd20
 cp  $POWR_WORK/dummychain/formal_${powr}_xxl                  $POWR_WORK/wrjobs/formal_${powr}_xxl
 if [[ $ihydro -gt 0 ]]; then
   cp  $POWR_WORK/dummychain/formal_${powr}_hydro                  $POWR_WORK/wrjobs/formal_${powr}_hydro
 fi
 echo "Copied formal_${powr}_gen, formal_${powr}_vd20 and formal_${powr}_xxl to wrjobs"
# fi


# njn --------------------------------------------------------------

touch                                                          $POWR_WORK/wrjobs/njn$kn
echo '#!/bin/bash'                                           > $POWR_WORK/wrjobs/njn$kn
echo 'source ~/.powrconfig'                                 >> $POWR_WORK/wrjobs/njn$kn 
echo '. $POWR_WORK/wrjobs/njn_'${powr}_gen "$kn"            >> $POWR_WORK/wrjobs/njn$kn
chmod 754                                                      $POWR_WORK/wrjobs/njn$kn

# if [ ! -e $POWR_WORK/wrjobs/njn_${powr}_gen ] ; then
 cp  $POWR_WORK/dummychain/njn_${powr}_gen                     $POWR_WORK/wrjobs/njn_${powr}_gen
 echo "Copied njn_${powr}_gen to wrjobs"
# fi

# set_repeat ---------------------------------------------------------

touch $POWR_WORK/wrjobs/set_repeat$kn
echo '#!/bin/bash'                > $POWR_WORK/wrjobs/set_repeat$kn
echo 'source ~/.powrconfig'      >> $POWR_WORK/wrjobs/set_repeat$kn
echo "# Set_Repeat$kn"           >> $POWR_WORK/wrjobs/set_repeat$kn
echo "#"                         >> $POWR_WORK/wrjobs/set_repeat$kn
echo 'cd $POWR_WORK/'"wrdata$kn" >> $POWR_WORK/wrjobs/set_repeat$kn
# old version 
# echo "echo 'REPEAT' > next_job" >> $POWR_WORK/wrjobs/set_repeat$kn
# echo "echo 'REPEAT' > next_jobz" >> $POWR_WORK/wrjobs/set_repeat$kn
# new: 13.07.2012 (Andreas S. -> EDDIX in MODEL missing)
echo "echo 'WRCONT' > next_job"  >> $POWR_WORK/wrjobs/set_repeat$kn
echo "echo 'WRCONT' > next_jobz" >> $POWR_WORK/wrjobs/set_repeat$kn
chmod 754 $POWR_WORK/wrjobs/set_repeat$kn
echo "Created set_repeat${kn} in wrjobs"

# steal -- dummycat

touch $POWR_WORK/wrjobs/steal$kn
echo '#!/bin/bash'                       > $POWR_WORK/wrjobs/steal$kn
echo 'source ~/.powrconfig'             >> $POWR_WORK/wrjobs/steal$kn
echo '# Steal'                          >> $POWR_WORK/wrjobs/steal$kn
echo "kn=$kn"                           >> $POWR_WORK/wrjobs/steal$kn
cat $POWR_WORK/dummychain/steal_${powr} >>  $POWR_WORK/wrjobs/steal$kn
chmod 754 $POWR_WORK/wrjobs/steal$kn
echo "Created steal${kn} in wrjobs"

# steal_backup -- dummycat

touch $POWR_WORK/wrjobs/steal${kn}_backup
echo '#!/bin/bash'                              > $POWR_WORK/wrjobs/steal${kn}_backup
echo 'source ~/.powrconfig'                    >> $POWR_WORK/wrjobs/steal${kn}_backup
echo '# Steal'                                 >> $POWR_WORK/wrjobs/steal${kn}_backup
echo "kn=$kn"                                  >> $POWR_WORK/wrjobs/steal${kn}_backup
cat $POWR_WORK/dummychain/steal_backup_${powr} >> $POWR_WORK/wrjobs/steal${kn}_backup
chmod 754 $POWR_WORK/wrjobs/steal${kn}_backup
echo "Created steal${kn}_backup in wrjobs"

# modify -- dummycat

touch $POWR_WORK/wrjobs/modify$kn
echo '#!/bin/bash'                        > $POWR_WORK/wrjobs/modify$kn
echo 'source ~/.powrconfig'              >> $POWR_WORK/wrjobs/modify$kn
echo "# Modify$kn"                       >> $POWR_WORK/wrjobs/modify$kn
echo "kn=$kn"                            >> $POWR_WORK/wrjobs/modify$kn
cat $POWR_WORK/dummychain/modify_${powr} >> $POWR_WORK/wrjobs/modify$kn
chmod 754                                   $POWR_WORK/wrjobs/modify$kn
echo "Created modify${kn} in wrjobs"

# # formal -- dummycat
# 
# touch $POWR_WORK/wrjobs/formal$kn
# echo '#!/bin/bash' > $POWR_WORK/wrjobs/formal$kn
# echo "# Formal$kn" >> $POWR_WORK/wrjobs/formal$kn
# echo "kn=$kn" >> $POWR_WORK/wrjobs/formal$kn
# cat $POWR_WORK/dummychain/formal_${powr} >> $POWR_WORK/wrjobs/formal$kn
# chmod 754 $POWR_WORK/wrjobs/formal$kn

echo "wrjobs for wrdata${kn} created"
# echo "Please add $kn to status.com"

done
# End of loop over all thenumbers

exit
