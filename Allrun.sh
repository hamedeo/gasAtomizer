#!/bin/bash

# List of directories to execute commands in
dir1=~/OpenFOAM/deo-v2212/run/secondary_breakup/secondaryBreakupBenzKHCase2
dir2=~/OpenFOAM/deo-v2212/run/secondary_breakup/secondaryBreakupBenzTabCase2
dir3=~/OpenFOAM/deo-v2212/run/secondary_breakup/secondaryBreakupBenzKHCase3
dir4=~/OpenFOAM/deo-v2212/run/secondary_breakup/secondaryBreakupBenzTabCase3
dir5=~/OpenFOAM/deo-v2212/run/secondary_breakup/secondaryBreakupBenzKHCase4
dir6=~/OpenFOAM/deo-v2212/run/secondary_breakup/secondaryBreakupBenzTabCase4
directories=($dir1 $dir2 $dir3 $dir4 $dir5 $dir6)

# Loop through each directory and execute commands
for dir in "${directories[@]}"
do 
    echo "Executing commands in $dir"
    cd "$dir" 
    sleep 5
    echo "Running foamCleanTutorials"   
    foamCleanTutorials
    echo "Running blockMesh "
    blockMesh 
    echo "Running decomposePar"
    decomposePar
    echo "Running mpirun -np 2 sprayFoam -parallel"
    mpirun -np 2 sprayFoam -parallel
    sleep 5
    echo "Running foamToVTK"
    foamToVTK
    sleep 5
    echo "Running foam.foam"
    touch foam.foam
    echo "Foam done"
done


echo "Foam done"
sleep 10

# Shutdown the system with sudo privileges and delay
# Prompt user for confirmation
read -p "Do you want to shutdown the system? (y/n) " confirm

# Wait for user input for 15 minutes
read -t 900 -n 1 -s -r -p "Press any key to cancel shutdown or wait 1 minute to shutdown..." response

# Check if user confirmed shutdown or no input was given within 15 minutes
if [[ "$confirm" == "y" ]] || [[ -z "$response" ]]
then
    echo "Shutting down the system in 30 seconds..."
    sleep 30
    sudo shutdown now
else
    echo "Shutdown cancelled."
fi
