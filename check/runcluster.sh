#!/usr/bin/env bash

# check if tmp-path exists
if test ! -d $CLIENTTMPDIR
then
    exit
fi

OUTFILE=$CLIENTTMPDIR/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/$BASENAME.err

uname -a                            > $OUTFILE
uname -a                            > $ERRFILE
if test $(uname) == Linux                   # -b does not work with top on macOS
then
    top -b -n 1 | head -n 15        >> $OUTFILE
fi
echo
echo @01 $INSTANCE ===========      >> $OUTFILE 
echo @01 $INSTANCE ===========      >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date +"@03 %s"                      >> $OUTFILE

if test "$SETTINGS" = ""
then
    echo $CURPATH/../$BINNAME $INSTANCE -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $OBJVAL >> $ERRFILE
    $CURPATH/../$BINNAME $INSTANCE -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $OBJVAL >> $OUTFILE 2>>$ERRFILE
else
    echo $CURPATH/../$BINNAME $INSTANCE -s $SETTINGS -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $OBJVAL >> $ERRFILE
    $CURPATH/../$BINNAME $INSTANCE -s $SETTINGS -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $OBJVAL >> $OUTFILE 2>>$ERRFILE
fi

retcode=${PIPESTATUS[0]}
if test "${retcode}" != 0
then
    echo "${INSTANCE} returned with error code ${retcode}." >>"${ERRFILE}"
fi

date +"@04 %s"                      >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $ERRFILE
echo                                >> $OUTFILE
echo =ready=                        >> $OUTFILE

mv $OUTFILE $CURPATH/results/$BASENAME.out
mv $ERRFILE $CURPATH/results/$BASENAME.err

#chmod g+r $ERRFILE
#chmod g+r $CURPATH/results/$BASENAME.out
#chmod g+r $CURPATH/results/$BASENAME.set
