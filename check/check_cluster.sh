#!/usr/bin/env bash
TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
DISPFREQ=$8
VERSION=$9
LPS=${10}
SETCUTOFF=${11}
QUEUE=${12}
PPN=${13}
CLIENTTMPDIR=${14}
EXCLUSIVE=${15}

# get current path
CURPATH=`pwd`

# directory for solutions
SOLDIR=$CURPATH/solutions
SOLUFILE=$CURPATH/testsets/${TSTNAME%%-}.solu
if test "${SETCUTOFF}" = 1 || test "${SETCUTOFF}" = true
then
    if test ! -e "${SOLUFILE}"
    then
	echo "Could not find solufile "${SOLUFILE}
	exit
    fi
fi

if test ! -e $CURPATH/results
then
    mkdir $CURPATH/results
fi

SETDIR=$CURPATH/../settings

# set settings name if necessary
if test "$SETNAME" = "default"
then
    SETTINGS=""
else
    SETTINGS=$SETDIR"/"$SETNAME".set"
fi

# check if the settings file exists
if test $SETNAME != "default"
then
    if test ! -e $SETTINGS
    then
	echo skipping test due to the non-existence of the settings file $SETTINGS
	exit
    fi
fi

# check if queue has been defined
if test "$QUEUE" = ""
then
    echo Skipping test since the queue name has not been defined.
    exit
fi

# check if number of nodes has been defined
if test "$PPN" = ""
then
    echo Skipping test since the number of nodes has not been defined.
    exit
fi

# check if client tmp-dir has been defined
if test "$CLIENTTMPDIR" = ""
then
    echo Skipping test since the path for the tmp-dir on the clients has not been defined.
    exit
fi

# check if the slurm blades should be used exclusively
if test "${EXCLUSIVE}" = "false"
then
    EXCLUSIVE=""
else
    EXCLUSIVE=" --exclusive"
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
# NOTE: the jobs should have a hard running time of more than 5 minutes; if not so, these
#       jobs get automatically assigned in the "exrpess" queue; this queue has only 4 CPUs
#       available 
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`

EVALFILE=$CURPATH/results/check.$QUEUE.$TSTNAME.$BINID.$SETNAME.eval
echo > $EVALFILE

# format time in seconds into format  dd-hh:mm:ss
formattime() {
    local T="${1}"
    local D=$((T / 60 / 60 / 24))
    local H=$((T / 60 / 60 % 24))
    local M=$((T / 60 % 60))
    local S=$((T % 60))
    local F="%02d"
    printf -v PRETTYTIME "${F}-${F}:${F}:${F}" "${D}" "${H}" "${M}" "${S}"
    echo "${PRETTYTIME}"
    }

ACCOUNT=""
if test "${CLUSTERQUEUE}" != "moskito" && test "${CLUSTERQUEUE}" != "prio"
then
    ACCOUNT="--account=dopt"

    # we add 10% to the hard memory limit and additional 100MB to the hard memory limit
    HARDMEMLIMIT=$(((MEMLIMIT + 100) + (MEMLIMIT / 10)))

    # we add 100% to the hard time limit and additional 600 seconds in case of small time limits
    HARDTIMELIMIT=$(((TIMELIMIT + 600) + TIMELIMIT))

    # format is (d-)HH:MM:SS
    HARDTIMELIMIT=$(formattime "${HARDTIMELIMIT}")
    # echo "${HARDTIMELIMIT}"
fi

echo
echo "testset:        " $TSTNAME
echo "queue:          " $QUEUE
echo "account:        " $ACCOUNT
echo "memlimit:       " $MEMLIMIT
echo "hard memlimit:  " $HARDMEMLIMIT
echo "timelimit:      " $TIMELIMIT
echo "hard timelimit: " $HARDTIMELIMIT
echo

# counter to define file names for a test set uniquely 
COUNT=1

for i in `cat testsets/$TSTNAME.test` DONE
do
  if test "$i" = "DONE"
  then
      break
  fi

  FILENAME=$i
  INSTANCE=$CURPATH/$FILENAME

  # check if problem instance exists 
  if test -e $INSTANCE
  then

      SHORTFILENAME=`basename $FILENAME .gz`
      SHORTFILENAME=`basename $SHORTFILENAME .mps`
      SHORTFILENAME=`basename $SHORTFILENAME .lp`
      SHORTFILENAME=`basename $SHORTFILENAME .opb`
      SHORTFILENAME=`basename $SHORTFILENAME .col`
      SHORTFILENAME=`basename $SHORTFILENAME .clq`
      SHORTFILENAME=`basename $SHORTFILENAME .txt`
      SHORTFILENAME=`basename $SHORTFILENAME .dimacs`
      SHORTFILENAME=`basename $SHORTFILENAME .wclq`

      BASENAME=$USER.$QUEUE.$TSTNAME.$COUNT"_"$SHORTFILENAME.$BINID.$SETNAME

      SETFILE=$CURPATH/results/$BASENAME.set

      echo $CURPATH/results/$BASENAME >> $EVALFILE

      # get objective value from solution file
      OBJVAL=""
      if test "${SETCUTOFF}" = 1 || test "${SETCUTOFF}" = true
      then
	  # get the objective value from the solution file: grep for the instance name and only use entries with an optimal or best known value;
	  # if there are multiple entries for this instance in the solution file, sort them by objective value and take the objective value
	  # written in the last line, i.e., the largest value;
	  # as a double-check, we do the same again, but reverse the sorting to get the smallest value
	  OBJECTIVEVAL=$(grep " ${SHORTFILENAME} " "${SOLUFILE}" | grep -e =opt= -e =best= | sort -k 3 -g | tail -n 1 | awk '{print $3}')

	  if test "${OBJECTIVEVAL}" = ""
	  then
	      echo "Could not find value for "${SHORTFILENAME}" in solufile "${SOLUFILE}
	      exit
	  fi
	  OBJVAL="-c "${OBJECTIVEVAL}
      fi

      if test "${CLUSTERQUEUE}" != "moskito" && test "${CLUSTERQUEUE}" != "prio"
      then
	 export CURPATH
	 export BINNAME
	 export INSTANCE
	 export BASENAME
	 export SOLNAME
	 export OBJVAL
	 export SETTINGS
	 export TIMELIMIT
	 export MEMLIMIT
	 export NODELIMIT
	 export DISPFREQ
	 export CLIENTTMPDIR

#	 echo sbatch --job-name="bacs_"$SHORTFILENAME --mem="${HARDMEMLIMIT}" --partition="${QUEUE}" "${ACCOUNT}" --time="${HARDTIMELIMIT}" --cpu-freq=highm1 ${EXCLUSIVE} --output=/dev/null runcluster.sh
         sbatch      --job-name="bacs_"$SHORTFILENAME --mem="${HARDMEMLIMIT}" --partition="${QUEUE}" "${ACCOUNT}" --time="${HARDTIMELIMIT}" --cpu-freq=highm1 ${EXCLUSIVE} --output=/dev/null runcluster.sh
      else
	 qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N $BINID":"$SHORTFILENAME -v CURPATH=$CURPATH,BINNAME=$BINNAME,INSTANCE=$INSTANCE,BASENAME=$BASENAME,SOLNAME=$SOLNAME,SETTINGS=$SETTINGS,TIMELIMIT=$TIMELIMIT,MEMLIMIT=$MEMLIMIT,NODELIMIT=$NODELIMIT,DISPFREQ=$DISPFREQ,CLIENTTMPDIR=$CLIENTTMPDIR -q $QUEUE -o /dev/null -e /dev/null runcluster.sh
      fi
      COUNT=`expr $COUNT + 1`
  else
      echo "input file "$INSTANCE" not found!"
  fi
done
