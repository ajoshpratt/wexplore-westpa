#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

TPR=CH4.tpr
TOP=CH4.top
NDX=CH4.ndx
OUT=_dist_$$.xvg

function cleanup() {
    rm -f $OUT
}

trap cleanup EXIT

# Get progress coordinate
rm mdout.mdp
rm $WEST_STRUCT_DATA_REF.tpr
rm \#*
grompp -f md-continue.mdp -c $WEST_STRUCT_DATA_REF.gro -o $WEST_STRUCT_DATA_REF.tpr -p $TOP -t $WEST_STRUCT_DATA_REF.trr -n $NDX -maxwarn 2
rm mdout.mdp
echo "6 7" | g_dist -f $WEST_STRUCT_DATA_REF.gro -s $WEST_STRUCT_DATA_REF.tpr -n $NDX -o $OUT -xvg none || exit 1
#awk '{print $2*10;}' < $OUT > $WEST_PCOORD_RETURN || exit 1
awk '{print $2*10;}' < $OUT > $$.PCOORD.1 || exit 1
echo "6 7" | g_dist -f $WEST_STRUCT_DATA_REF.gro -s $WEST_STRUCT_DATA_REF.tpr -n $NDX -o $OUT -xvg none || exit 1
awk '{print $2*10;}' < $OUT > $$.PCOORD.2 || exit 1
paste $$.PCOORD.1 $$.PCOORD.2 > $WEST_PCOORD_RETURN || exit 1

if [ -n "$SEG_DEBUG" ] ; then
    head -v $WEST_PCOORD_RETURN
fi

rm $$.PCOORD.*
