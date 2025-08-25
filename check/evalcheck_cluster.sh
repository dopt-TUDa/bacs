#!/usr/bin/env bash

export LANG=C

REMOVE=0
AWKARGS=""
FILES=""

for i in $@
do
  if test ! -e $i
  then
      if test "$i" = "-r"
      then
	  REMOVE=1
      else
	  AWKARGS="$AWKARGS $i"
      fi
  else	
      FILES="$FILES $i"
  fi
done

for FILE in $FILES
do
 
  DIR=`dirname $FILE`
  EVALFILE=`basename $FILE .eval`
  EVALFILE=`basename $EVALFILE .out`

  OUTFILE=$DIR/$EVALFILE.out 
  ERRFILE=$DIR/$EVALFILE.err
  SETFILE=$DIR/$EVALFILE.set
  RESFILE=$DIR/$EVALFILE.res
  TEXFILE=$DIR/$EVALFILE.tex
  PAVFILE=$DIR/$EVALFILE.pav
  
  # check if the eval file exists; if this is the case construct the overall solution files
  if test -e $DIR/$EVALFILE.eval
  then
      echo > $OUTFILE
      echo > $ERRFILE
      echo create overall output and error file for $EVALFILE

      for i in `cat $DIR/$EVALFILE.eval` DONE
	do
	if test "$i" = "DONE"
	then
	    break
	fi
	
	FILE=$i.out
	if test -e $FILE
	then
	    cat $FILE >> $OUTFILE
	    if test "$REMOVE" = "1"
	    then
		rm -f $FILE
	    fi
	else
	    echo Missing $i
	fi
      
	FILE=$i.err
	if test -e $FILE
	then
	    cat $FILE >> $ERRFILE
	    if test "$REMOVE" = "1"
	    then
		rm -f $FILE
	    fi
	fi
      
	FILE=$i.set
	if test -e $FILE
	then
	    cp $FILE $SETFILE
	    if test "$REMOVE" = "1"
	    then
		rm -f $FILE
	    fi
	fi
      
	FILE=$i.tmp
	if test -e $FILE
        then
	    if test "$REMOVE" = "1"
	    then
		rm -f $FILE
	    fi
	fi
      done
      
      if test "$REMOVE" = "1"
      then
	  rm -f $DIR/$EVALFILE.eval
      fi
  fi

  # check if the out file exists
  if test -e $DIR/$EVALFILE.out
  then

      echo create results for $EVALFILE

      # detect used queue
      QUEUE=`echo $EVALFILE | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g'`

      # detect test set
      if test "$QUEUE" = ""
      then
          TSTNAME=$QUEUE
      else
          TSTNAME=`echo $EVALFILE | sed 's/check.'$QUEUE'.\([a-zA-Z0-9_-]*\).*/\1/g'`
      fi

      # detect test used solver
      SOLVER=`echo $EVALFILE | sed 's/check.\([a-zA-Z0-9_-]*\).\([a-zA-Z0-9_-]*\).\([a-zA-Z0-9_]*\).*/\3/g'`
      
      echo "Queue   " $QUEUE
      echo "Testset " $TSTNAME
      echo "Solver  " $SOLVER

      if test -f testsets/$TSTNAME.test
      then
	  TESTFILE=testsets/$TSTNAME.test
      else
	  TESTFILE=""
      fi

      if test -f testsets/$TSTNAME.solu
      then
	  SOLUFILE=testsets/$TSTNAME.solu
      else 
	  if test -f testsets/all.solu
	  then
	      SOLUFILE=testsets/all.solu
	  else
	      SOLUFILE=""
	  fi
      fi

      if test  "$SOLVER" = "cplex"
      then
	  awk -f check_cplex.awk -v "TEXFILE=$TEXFILE" $AWKARGS $SOLUFILE $OUTFILE | tee $RESFILE
      else
	  awk -f check.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" $AWKARGS $TESTFILE $SOLUFILE $OUTFILE | tee $RESFILE
      fi
  fi
done
