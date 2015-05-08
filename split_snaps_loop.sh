cd ..
for d in king_w*/ ; do
    echo "$d";
    cd $d
    pwd
    mkdir snaps
    python ../nbody_analysis/split_snaps.py
    cd ..
    echo "... done!"
done
cd nbody_analysis
