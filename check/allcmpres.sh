#!/usr/bin/env bash

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

TESTSETS=""
# for i in `ls -1 --color=none $FILES | sed 's!\(.*\)check\.\([^ .]*\)\.\([^ ]*\)\.res!\2!g' | sort -u`
for i in `ls -1 --color=none $FILES | sed 's!\(.*\)check\.\([^ .]*\)\.\([^ ]*\)\.res!\2!g'`
do
    TESTSETS="$TESTSETS $i"
done

export LC_NUMERIC=C

for i in $TESTSETS
do
    echo
    echo ====vvvv==== $i ====vvvv====
    awk -f cmpres.awk $AWKARGS texcmpfile="cmpres.$i.tex" `ls -1 --color=none $FILES | grep "$i\..*\.res"`
    echo ====^^^^==== $i ====^^^^====
done
