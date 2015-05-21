cd ../100K/
for d in king_w7bh*/ ; do
    echo "$d";
    cd $d
    pwd
    sh ../../nbody_analysis/group_data.sh
    cd ..
    echo "... done!"
done
cd nbody_analysis