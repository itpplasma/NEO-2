basename_directories="n2_vshift"

matlabscript="path('/afs/itp.tugraz.at/user/buchholz/Programs/neo-2/OctaveScripts', path); try, export_2spec_Matyas, catch, exit(1), end, exit(0);"
octavescript="export_2spec_Matyas;"
target=$(pwd)

dirs=$(ls -d ${basename_directories}*)

# for each folder wih a velocity shift ...
for d in ${dirs} ; do
    echo "$d"
    cd $d

    # ... merge the neo2_multispecies output files ...
    ./h5merge_multispec.x
    eval 'matlab -nodesktop -nosplash -r "$matlabscript"'
    #~ eval 'octave --no-gui --eval "$octavescript"'

    if (( $? == 0 )) ; then
        echo "$d:    OK"
    else
        echo "$d:    failed"
    fi

    cd $target

done

get_NTV_torque_int.sh
