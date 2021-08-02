#!/bin/bash
 
#Build script for interThermalPhaseChangeFoam
wclean libso incompressibleTwoPhaseThermalMixture
wclean interThermalPhaseChangeFoam
wclean ./Libraries/DynamicKistlerContactAngle
wmake libso incompressibleTwoPhaseThermalMixture
wmake interThermalPhaseChangeFoam
cd Libraries
./Allwmake.sh
