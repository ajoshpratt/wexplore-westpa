#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -pv $WEST_CURRENT_SEG_DATA_REF || exit 1
cd $WEST_CURRENT_SEG_DATA_REF || exit 1

TOP=CH4.top
NDX=CH4.ndx

if [[ "$USE_LOCAL_SCRATCH" == "1" ]] ; then
    # make scratch directory
    WORKDIR=/$SCRATCH/$USER/$WEST_JOBID/$WEST_CURRENT_SEG_DATA_REF
    $SWROOT/bin/mkdir -p $WORKDIR || exit 1
    cd $WORKDIR || exit 1
    STAGEIN="$SWROOT/bin/rsync -a"
else
    STAGEIN="$SWROOT/bin/ln -sv"
fi

function cleanup() {
    # Clean up
    if [[ "$USE_LOCAL_SCRATCH" == "1" ]] ; then
        $SWROOT/bin/cp seg.{xtc,trr,edr,tpr,gro,log} whole.xtc nojump.xtc $WEST_CURRENT_SEG_DATA_REF || (sleep 10 && $SWROOT/bin/cp seg.{xtc,trr,edr,tpr,gro,log} whole.xtc nojump.xtc $WEST_CURRENT_SEG_DATA_REF) || exit 1
        cd $WEST_CURRENT_SEG_DATA_REF
        $SWROOT/bin/rm -Rf $WORKDIR
    else
        $SWROOT/bin/rm -f *.xvg *.itp *.mdp *.ndx *.top state.cpt
    fi
}

trap cleanup EXIT

# Set up the run
ln -sv $WEST_SIM_ROOT/{$TOP,$NDX,*.mdp,*.itp} .

# if [ "$WEST_PARENT_SEG_ID" -lt "0" ]; then
#     # this is a start/restart
#     ln -sv $WEST_PARENT_SEG_DATA_REF/unbound.gro .
#     grompp -f md-genvel.mdp -c unbound.gro -p $TOP -o seg.tpr -n $NDX || exit 1
# else
#     ln -sv $WEST_PARENT_SEG_DATA_REF/seg.gro ./parent.gro
#     ln -sv $WEST_PARENT_SEG_DATA_REF/seg.trr ./parent.trr
#     grompp -f md-continue.mdp -c parent.gro -t parent.trr -p $TOP -o seg.tpr -n $NDX || exit 1
# fi

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        ln -sv $WEST_PARENT_DATA_REF/seg.gro ./parent.gro
        ln -sv $WEST_PARENT_DATA_REF/seg.trr ./parent.trr
        ln -sv $WEST_PARENT_DATA_REF/seg.edr ./parent.edr
        grompp -f md-continue.mdp -c parent.gro -t parent.trr -p $TOP -o seg.tpr -n $NDX -maxwarn 2 -e parent.edr || exit 1
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory; $WEST_PARENT_DATA_REF contains the reference to the
        # appropriate basis state or generated initial state
        ln -sv $WEST_PARENT_DATA_REF.gro ./initial.gro
        ln -sv $WEST_PARENT_DATA_REF.trr ./initial.trr
        ln -sv $WEST_PARENT_DATA_REF.edr ./initial.edr
        grompp -f md-continue.mdp -c initial.gro -t initial.trr -p $TOP -o seg.tpr -n $NDX -maxwarn 2 -e initial.edr || exit 1
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac
    

# Propagate segment
mdrun -s seg.tpr -o seg.trr -x seg.xtc -c seg.gro \
      -e seg.edr -g seg.log -nt 1 -ntomp 1 -ntomp_pme 1 \
      || exit 1

# Get progress coordinate
echo -e "2 \n" | trjconv    -f seg.xtc  -s seg.tpr  -n $NDX -o whole.xtc -pbc nojump || exit 1
echo "6 7" | g_dist -f whole.xtc -s seg.tpr -n $NDX -xvg none || exit 1
awk '{print $2*10;}' < dist.xvg > $$.PCOORD.1 || exit 1
echo "6 7" | g_dist -f whole.xtc -s seg.tpr -n $NDX -xvg none || exit 1
awk '{print $2*10;}' < dist.xvg > $$.PCOORD.2 || exit 1
paste $$.PCOORD.1 $$.PCOORD.2 > $WEST_PCOORD_RETURN || exit 1
#rm $$.PCOORD.*
