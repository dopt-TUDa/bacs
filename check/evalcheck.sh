#!/bin/bash

export LANG=C

AWKARGS=""
FILES=""
for i in $@
do
    if test ! -e $i
    then
	AWKARGS="$AWKARGS $i"
    else
	FILES="$FILES $i"
    fi
done

for i in $FILES
do
    NAME=`basename $i .out`
    DIR=`dirname $i`
    OUTFILE=$DIR/$NAME.out
    RESFILE=$DIR/$NAME.res
    TEXFILE=$DIR/$NAME.tex
    PAVFILE=$DIR/$NAME.pav
    ERRFILE=$DIR/$NAME.err

    TSTNAME=`echo $NAME | sed 's/check.\([a-zA-Z0-9_]*\).*/\1/g'`

    if test -f testsets/$TSTNAME.test
    then
	TESTFILE=testsets/$TSTNAME.test
    else
	TESTFILE=""
    fi

    if test -f testsets/$TSTNAME.solu
    then
	SOLUFILE=testsets/$TSTNAME.solu
    else if test -f testsets/all.solu
    then
	SOLUFILE=testsets/all.solu
    else
        SOLUFILE=""
    fi
    fi

    gawk -f check.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" -v "ERRFILE=${ERRFILE}" $AWKARGS $TESTFILE $SOLUFILE $OUTFILE | tee $RESFILE
done
