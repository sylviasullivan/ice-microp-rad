# This version assumes that you start with files that are not yet
# horizontally regridded.

# Build a command to remap the model level output to the pressure levels
# given in PMEAN_*-*.txt

# Directory and inputs to build the file names.
basedir='/scratch/b/b380873/tropic_run5/'
echo 'File prefix, e.g. WINDTH_3D, RAD_3D'
read fileprefix
echo 'First timestep, e.g. 1, 25, 52'
read lowtimestep
echo 'Last timestep'
read hightimestep

for timestep in $(seq $lowtimestep $hightimestep); do
    echo $timestep

    # How many leading zeros for the timestep in the filename?
    if [ ${#timestep} -eq 1 ]; then
       timestepprefix='000'
    elif [ ${#timestep} -eq 2 ]; then
       timestepprefix='00'
    elif [ ${#timestep} -eq 3 ]; then
       timestepprefix='0'
    else
       echo 'Incorrect time step specified.'
    fi

    # Read the pressure levels from the txt file. Add them line by line to $pressures.
    # c counts how many pressures are specified.
    #c=0
    #pressures=
    #while read -r line; do pressures=$pressures$line','; c=$((c+1)); done < PMEAN_48-72.txt
    #echo $c' pressures specified'

    # Remove the last comma from $pressures.
    #pressures="${pressures%?}"

    # Simpler route. Pressures must [=] Pa.
    pressures='90000,85000,75000,50000,25000'

    # Assemble the filenames and command.
    part1='cdo -ap2pl,'
    part2='-selvar,w,air_pressure,pres_sfc'
    inputfile=$basedir$fileprefix'_icon_tropic_'$timestepprefix$timestep'.nc'
    outputfile=$basedir'W_3D_icon_tropic_'$timestepprefix$timestep'_PL.nc'

    # Rename pres as air_pressure in the $inputfile.
    cdo chname,pres,air_pressure $inputfile $basedir'temp.nc'
    cmd=$part1$pressures' '$part2' '$basedir'temp.nc '$outputfile

    # Evaluate the command as you would from the command line.
    eval $cmd
done
