#!/bin/bash


# Original directory
orDir=$(pwd)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Directory List  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# List of directories to execute commands in
directories=$(find ./ -maxdepth 1 \( -name "secondaryBreakupBenz*" \) -not \( -path "./secondaryBreakupBenz_BaseCase" -o -path "./secondaryBreakupBenz_Test" \))
echo "$directories"
# Prompt the user to confirm shutdown
read -t 10 -n 2 -p "Do you confirm the directory list? (y/n) " dirlistconfirm

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



# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Command to execute  * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

# Create output directory
mkdir "./output_dir"

# Process each directory
for dir in $directories
do  
    mkdir "output_dir/$dir"
    mkdir "output_dir/$dir/lagrangian"
    echo "Executing commands in $dir"
    cp -r "$dir/processor0/VTK/lagrangian" "output_dir/$dir"
    cp -r "$dir/processor1/VTK/lagrangian" "output_dir/$dir"

    cd $orDir
done
# ************************************************************************* #


echo "Done"























