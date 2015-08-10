#!/bin/bash

PARENTDIR=$1
OBJECT=$2
DATE=$3
CALS=$4
RAWS=$5
MDAR=$6
MFLA=$7
DSUB=$8
SUBS=$9

PARENTDIR=/mnt/gfsproject/naiad/njones/moddragonflydata/$OBJECT
CALFILESDIR=$PARENTDIR/$CALS/$DATE
RAWLIGHTSDIR=$PARENTDIR/$RAWS/$DATE
MASTERDARKDIR=$PARENTDIR/$MDAR/$DATE
MASTERFLATDIR=$PARENTDIR/$MFLA/$DATE
DARKSUBDIR=$PARENTDIR/$DSUB/$DATE
FULLCALDIR=$PARENTDIR/$SUBS/$DATE

if [ ! -d "$CALFILESDIR"]; then
mkdir $CALFILESDIR
fi
if [ ! -d "$MASTERDARKDIR"]; then
mkdir $MASTERDARKDIR
fi
if [ ! -d "$MASTERFLATDIR"]; then
mkdir $MASTERFLATDIR
fi
if [ ! -d "$DARKSUBDIR"]; then
mkdir $DARKSUBDIR
fi
if [ ! -d "$FULLCALDIR"]; then
mkdir $FULLCALDIR
fi
if [ ! -d "$RAWLIGHTSDIR"]; then
mkdir $RAWLIGHTSDIR
fi

echo Create master darks
create_masterdarks -v -o $MASTERDARKDIR $CALFILESDIR
echo Subtract darks from raw light
subtract_masterdarks -v -m $MASTERDARKDIR -o $DARKSUBDIR $RAWLIGHTSDIR
echo Subtract darks from flats
subtract_masterdarks -v -m $MASTERDARKDIR -o $DARKSUBDIR $CALFILESDIR
echo Create master flats
create_masterflats -v -o $MASTERFLATDIR $DARKSUBDIR
echo Divide by flats
divideby_masterflats -v -m $MASTERFLATDIR -o $FULLCALDIR $DARKSUBDIR