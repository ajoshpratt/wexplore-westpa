#!/bin/bash -l
source env.sh
ps aux | grep w_run | grep -v grep
pkill -9 -f w_run

SFX=.d$$
mv /tmp/traj_segs{,$SFX}
mv /tmp/seg_logs{,$SFX}
mv istates{,$SFX}
rm -Rf /tmp/traj_segs$SFX /tmp/seg_logs$SFX istates$SFX & disown %1
rm -f system.h5 west.h5 seg_logs.tar
mkdir /tmp/seg_logs /tmp/traj_segs istates
ln -sv /tmp/traj_segs .
ln -sv /tmp/seg_logs .

BSTATE_ARGS="--bstates-from BASIS_STATES"
TSTATE_ARGS="--tstate unbound,10.0,1.0"

$WEST_ROOT/bin/w_init $BSTATE_ARGS $TSTATE_ARGS --segs-per-state 1 --serial
