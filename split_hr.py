#
# SPLIT the HR diagnostic of singles into individual files
#
# need to: 
# mkdir hr
#
#

out = 0
count = 1
timeout = open("hr/times.dat","w")
with open("fort.83") as f: 
    for line in f:
        #print line[0:9]
        if line[0:9] == " ## BEGIN":
            if out: out.close()
            out = open("hr/hr_single_%05i.dat" % count, "w")
            timeout.write( line.split()[3] )
            timeout.write("\n")
            count = count+1
        elif line[0:9] != " ## END  ": 
            out.write(line)


#            timeout.write( line.split()[1] )
#            timeout.write("\n")
