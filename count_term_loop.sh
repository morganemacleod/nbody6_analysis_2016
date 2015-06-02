cd ../100K/
for d in king_w7*/ ; do
    echo "$d";
    cd $d
    pwd
    python ../../nbody_analysis/count_termination.py > termination_stats.txt
    cd ..
    echo "... done!"
done
cd ../nbody_analysis