cd ../100K/
for d in king_w7_A/ ; do
    echo "$d";
    cd $d
    pwd
    python ../../nbody_analysis/read_snaps.py
    cd ..
    echo "... done!"
done
cd ../nbody_analysis