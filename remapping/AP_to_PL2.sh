# This version assumes that you start with files that have already
# been vertically regridded.

# Build a command to remap the model level output to the pressure levels
# given in PMEAN_*-*.txt

# Directory and inputs to build the file names.
basedir='/scratch/b/b380873/1V2M1A1R/'
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
    c=0
    pressures=
    while read -r line; do pressures=$pressures$line','; c=$((c+1)); done < PMEAN_48-72.txt
    echo $c' pressures specified'

    # Remove the last comma from $pressures.
    pressures="${pressures%?}"

    # Simpler route. Pressures must [=] Pa.
    #pressures='90000,85000,75000,50000,25000'
    #pressures='5000,7000,10000,15000,20000,25000,30000,40000,50000,60000,70000,80000,85000,90000,92500,95000,97500,100000'

    # Assemble the filenames and command.
    part1='cdo -ap2pl,'
#    part2='-selvar,w,air_pressure,pres_sfc'
    inputfile=$basedir$fileprefix'_icon_tropic_'$timestepprefix$timestep'.nc'
    outputfile=$basedir$fileprefix'_icon_tropic_'$timestepprefix$timestep'_PL2.nc'

    # Rename pres as air_pressure in the $inputfile.
    cdo chname,pres,air_pressure $inputfile $basedir'temp.nc'
    cmd=$part1$pressures' '$part2' '$basedir'temp.nc '$outputfile

    # Evaluate the command as you would from the command line.
    echo $cmd
    eval $cmd
done
