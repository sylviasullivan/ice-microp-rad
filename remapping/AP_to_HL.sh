# This version assumes that you start with files that have already
# been vertically regridded.

# Build a command to remap the model level output to the heights
# given in HEIGHTS_*.txt

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

    # Read the heights from the txt file. Add them line by line to $heights.
    # c counts how many heights are specified.
    c=0
    heights=
    #while read -r line; do heights=$heights$line','; c=$((c+1)); done < HEIGHTS_fig2.txt
    #while read -r line; do heights=$heights$line','; c=$((c+1)); done < HEIGHTS_CloudSat.txt
    while read -r line; do heights=$heights$line','; c=$((c+1)); done < HEIGHTS_full_range.txt
    echo $c' heights specified'

    # Remove the last comma from $heights.
    heights="${heights%?}"

    # Simpler route. Heights must [=] m.
    #heights='1000,2000,3000,...'

    # Assemble the filenames and command.
    part1='cdo -ap2hl,'
#    part2='-selvar,w,air_pressure,pres_sfc'
    inputfile=$basedir$fileprefix'_icon_tropic_'$timestepprefix$timestep'_remapdis_0.025.nc'
    outputfile=$basedir$fileprefix'_icon_tropic_'$timestepprefix$timestep'_remapdis_0.025_HL.nc'

    # Rename pres as air_pressure in the $inputfile.
    echo 'Change pres to air_pressure'
    cdo chname,pres,air_pressure $inputfile $basedir'temp.nc'
    cmd=$part1$heights' '$part2' '$basedir'temp.nc '$outputfile

    # Evaluate the command as you would from the command line.
    #echo $cmd
    eval $cmd
done
