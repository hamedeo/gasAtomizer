#!/bin/bash


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Prevent screen from blanking
xset s off
xset -dpms
xset s noblank
# ************************************************************************* #

wait_time=1
n=0

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Prompt the user to confirm shutdown
read -p "Do you want to shut down the system after execution? (y/n) " confirm

if [[ $confirm == "y" ]]; then
    # Prompt the user to enter the wait time
    read -p "Enter the time to wait before shutting down (in minutes): " wait_time
    echo "System will shut down in $wait_time minutes after execution."
fi
# ************************************************************************* #


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Case Execution  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
start_time=$(date +%s.%N)

# Original directory
orDir=$(pwd)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Directory List  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# List of directories to execute commands in
directories=$(find ./ -maxdepth 1 \( -name "secondaryBreakupBenz_*" \) -not \( -path "./secondaryBreakupBenz_baseCase" -o -path "./secondaryBreakupBenz_Test" -o -path "./secondaryBreakupBenz_KH_Case2" \))
echo "$directories"
sleep 30
# Prompt the user to confirm shutdown
read -p "Do you confirm the directory list? (y/n) " dirlistconfirm

if [[ $dirlistconfirm != "y" ]]; then
    echo "Please check the code again"
    sleep 2
    exit 1
else
    echo "Great"
    sleep 1
fi

# directories=($dir1 $dir2 $dir3 $dir4 $dir5 $dir6)
# ************************************************************************* #


# Loop through each directory and execute commands
for dir in $directories
do  
    n=$((n+1))
    echo "Executing commands in $dir"
    cd "$dir" 
    sleep 5
    echo "clean the case"   
    foamCleanTutorials
    
    # Start time of each case
    start=$(date +%s.%N)                

    echo "Running blockMesh && decomposePar "
    blockMesh && decomposePar    

    echo "Running mpirun -np 2 sprayFoam -parallel"
    mpirun -np 2 sprayFoam -parallel

    sleep 5
    echo "Translate foam to VTK in each processor file"
    cd ./processor0
    foamToVTK
    cd ../processor1
    foamToVTK
    cd ..

    echo "Create foam.foam"
    touch foam.foam

    echo "Foam $dir done"

    #  End time of each case 
    end=$(date +%s.%N)                  
    
    # Beep superMario as a case complition
    beep -n -f 146 -l 20 -D 0 -n -f 369 -l 20 -D 0 -n -f 659 -l 20 -D 47 -n -f 146 -l 20 -D 0 -n -f 369 -l 20 -D 0 -n -f 659 -l 20 -D 189 -n -f 146 -l 20 -D 0 -n -f 369 -l 20 -D 0 -n -f 659 -l 20 -D 190 -n -f 146 -l 20 -D 0 -n -f 369 -l 20 -D 0 -n -f 523 -l 20 -D 47 -n -f 146 -l 20 -D 0 -n -f 369 -l 20 -D 0 -n -f 659 -l 20 -D 190 -n -f 391 -l 20 -D 0 -n -f 493 -l 20 -D 0 -n -f 783 -l 20 -D 475 -n -f 195 -l 20 -D 0 -n -f 391 -l 20 -n -f 195 -l 20 -D 0 -n -f 391 -l 20 -D 473
    
    elapsed=$(echo "scale=3; ($end - $start) / 1" | bc)   # Calculate the elapsed time in seconds
    echo "case$n execution time: $elapsed seconds"
    # Back to original directory
    cd $orDir
done

end_time=$(date +%s.%N)
total_time=$(echo "$end_time - $start_time" | bc)
echo "Whole execution time $total_time seconds"

echo "All Foam done"

# Play Victory as BEEP
beep -l 375 -f 392 -n -l 750 -f 523 -n -l 463 -f 392 -n -l 187 -f 440 -n -l 750 -f 494 -n -l 375 -f 330 -n -l 375 -f 330 -n -l 750 -f 440 -n -l 463 -f 392 -n -l 187 -f 349 -n -l 750 -f 392 -n -l 463 -f 262 -n -l 187 -f 262 -n -l 750 -f 294 -n -l 463 -f 294 -n -l 187 -f 330 -n -l 750 -f 349 -n -l 463 -f 349 -n -l 187 -f 392 -n -l 750 -f 440 -n -l 375 -f 494 -n -l 375 -f 523 -n -l 1125 -f 587
sleep 10
# ************************************************************************* #

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Turn off screen blanking prevention
xset s on
xset +dpms
xset s blank
# ************************************************************************* #


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Shutdown the system
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Shutdown the system with sudo privileges and delay
if [[ "$confirm" != "y" ]] then
    # Wait for user input for $wait_time minutes
    read -t $(($wait_time * 60)) -n 2 -s -r -p "Press any key to cancel shutdown or wait $wait_time minute to shutdown..." response

    # Check if user confirmed shutdown or no input was given within $wait_time minutes
    if [[ $response != "" ]] then
        echo "Shutdown cancelled."
    else
        echo "Shutting down the system in 30 seconds..."
        sleep 30
        sudo shutdown now
    fi

else
    echo "Enjoy results"
fi
# ************************************************************************* #