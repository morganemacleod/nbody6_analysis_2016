cd ..
for d in king_w*/ ; do
    echo "$d";
    cd $d
    pwd
    mkdir hr
    python ../nbody_analysis/split_hr.py
    cd ..
    echo "... done!"
done
cd nbody_analysis
