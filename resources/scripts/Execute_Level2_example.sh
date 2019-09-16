#!/bin/sh
#
# Make sure to set the bath to your I3_BUILD correctly:
cd /home/mduvoort/icetray/icerec/build
./env-shell.sh python IcePacker/resources/scripts/Level2_example.py -g /data/exp/IceCube/2008/filtered/level2/burnSample/PFFilt_PhysicsTrig_PhysicsFiltering_Run00110750_Level2_GCD.i3.gz -i /data/exp/IceCube/2008/filtered/level2/burnSample/PFFilt_PhysicsTrig_PhysicsFiltering_Run00110750_Level2_Part00000000.i3.gz,/data/exp/IceCube/2008/filtered/level2/burnSample/PFFilt_PhysicsTrig_PhysicsFiltering_Run00110750_Level2_Part00000001.i3.gz -o PFFilt_PhysicsTrig_PhysicsFiltering_Run00110750_Level2.IceEvent.root -n 100
