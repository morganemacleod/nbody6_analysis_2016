# This script is to be run in the simulation output directory. It gets the 
#  output of the output<n> folders and groups it together into single, larger, files.
#
#
#
#echo "ADJUST ... "
#grep ADJUST output*/output.dat -h > adjust.dat
#
#echo "bh..."
#grep found output*/output.dat -h | cut -c 16- > bh.dat
#
#echo "fort.40 ..."
#grep "" output*/fort.40 > fort.40 -h
#
#echo "fort.50 ..."
#grep "" output*/fort.50 -h | cut -c 10- > fort.50
#
#echo "fort.57 ..."
#grep "" output*/fort.57 > fort.57 -h
#
#echo "fort.58 ..."
#grep "" output*/fort.58 > fort.58 -h 
#
#echo "fort.45 ..."
#grep "" output*/fort.45 > fort.45 -h
#
#echo "fort.82 ..."
#grep "" output*/fort.82 > fort.82 -h
#
#echo "fort.83 ..."
#grep "" output*/fort.83 > fort.83 -h
#
#echo "snapdata ..."
#grep "" output*/snapdata.dat > snapdata.dat -h
#
#echo "mttime ..."
#grep "Time:" output*/output.dat -h | cut -c 8- > mttime.dat
#
#echo "core/hm ..."
#grep "hmr" output*/output.dat -h | cut -c 32-  > core_hm_radii.dat
#
#echo "inspirals ..."
#grep INSPIRAL output*/output.dat -h | cut -c 19- > inspiral.dat
#
#echo "scales ..."
#grep PHYSICAL output1/output.dat > scales.dat
# 
echo "t= ..."
grep "T =" output*/output.dat > time_summary.dat