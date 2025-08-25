#!/bin/bash
TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
DISPFREQ=$8
CONTINUE=$9
LOCK={$10}
VERSION=${11}
LPS=${12}
USEBESTSOL=${13}
SETCUTOFF=${14}

SETDIR=../settings

if test ! -e results
then
    mkdir results
fi
if test ! -e locks
then
    mkdir locks
fi

LOCKFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.lock
RUNFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.run.$BINID
DONEFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.done

OUTFILE=results/check.$TSTNAME.$BINID.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINID.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINID.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINID.$SETNAME.tex
SETFILE=results/check.$TSTNAME.$BINID.$SETNAME.set

SOLDIR=solutions

SOLUFILE=testsets/${TSTNAME%%-*}.solu

if test "$LOCK" = "true"
then
    if test -e $DONEFILE
    then
	echo skipping test due to existing done file $DONEFILE
	exit
    fi
    if test -e $LOCKFILE
    then
	if test -e $RUNFILE
        then
	    echo continuing aborted run with run file $RUNFILE
	else
	    echo skipping test due to existing lock file $LOCKFILE
	    exit
	fi
    fi
    date > $LOCKFILE
    date > $RUNFILE
fi

if test ! -e $OUTFILE
then
    CONTINUE=false
fi

if test "$CONTINUE" = "true"
then
    MVORCP=cp
else
    MVORCP=mv
fi

DATEINT=`date +"%s"`
if test -e $OUTFILE
then
    $MVORCP $OUTFILE $OUTFILE.old-$DATEINT
fi
if test -e $ERRFILE
then
    $MVORCP $ERRFILE $ERRFILE.old-$DATEINT
fi

if test "$CONTINUE" = "true"
then
    LASTPROB=`./getlastprob.awk $OUTFILE`
    echo Continuing benchmark. Last solved instance: $LASTPROB
    echo "" >> $OUTFILE
    echo "----- Continuing from here. Last solved: $LASTPROB -----" >> $OUTFILE
    echo "" >> $OUTFILE
else
    LASTPROB=""
fi

uname -a >>$OUTFILE
uname -a >>$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

# set settings name if necessary
if test "$SETNAME" = "default"
then
    SETTINGS=""
else
    SETTINGS="-s "$SETDIR"/"$SETNAME".set"

    if test ! -e $SETTINGS
    then
	echo skipping test due to the non-existence of the settings file $SETTINGS
	exit
    fi
fi

# process test instances
STOREFILENAME=""
for i in `cat testsets/$TSTNAME.test` DONE
do
    if test "$i" = "DONE"
    then
	date > $DONEFILE
	break
    fi

    FILENAME=$i

    if test "$LASTPROB" = ""
    then
	LASTPROB=""
	if test -f $FILENAME
	then
	    echo @01 $FILENAME ===========
	    echo @01 $FILENAME ===========                >> $ERRFILE

	    if test "$USEBESTSOL" = "true" 
	    then
		SOLBNAME=`basename $FILENAME .gz`
		SOLBNAME=`basename $SOLBNAME .col`
		SOLBNAME=`basename $SOLBNAME .clq`
		SOLNAME=$SOLDIR"/"$SOLBNAME.sol
	    else
		SOLNAME=""
	    fi

	    # get objective value from solution file
	    OBJVAL=""
	    if test "${SETCUTOFF}" = 1 || test "${SETCUTOFF}" = true
	    then
		if test -e "${SOLUFILE}"
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
		else
		    echo "${SOLUFILE} not found"
		fi
	    fi

	    echo -----------------------------
	    date
	    date >>$ERRFILE
	    echo -----------------------------
	    date +"@03 %s"
	    bash -c " ulimit -f 200000; ../$BINNAME $FILENAME $SETTINGS -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $OBJVAL" 2>>$ERRFILE

	    retcode=$PIPESTATUS
	    if test "${retcode}" != 0
	    then
	       echo "${FILENAME} returned with error code ${retcode}." >>"${ERRFILE}"
	    fi

	    date +"@04 %s"
	    echo -----------------------------
	    date
	    date >>$ERRFILE
	    echo -----------------------------
	    echo
	    echo =ready=
	else
	    echo @02 FILE NOT FOUND: $i ===========
	    echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
	fi
    else
	echo skipping $FILENAME
	if test "$LASTPROB" = "$FILENAME"
	then
	    LASTPROB=""
        fi
    fi
done | tee -a $OUTFILE

date >>$OUTFILE
date >>$ERRFILE

if test -e $DONEFILE
then
    ./evalcheck.sh $OUTFILE

    if test "$LOCK" = "true"
    then
	rm -f $RUNFILE
    fi
fi
