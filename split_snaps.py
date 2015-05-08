out = 0
count = 1
with open("snapdata.dat") as f: 
    for line in f:
        if line[0] == "t":
            if out: out.close()
            out = open("snaps/snap_%05i.dat" % count, "w")
            count = count+1
            print "writing ... snaps/snap_%05i.dat" % count
        else: out.write(line)

print "about to write times..."
timeout = open("snaps/times.dat","w")
with open("snapdata.dat") as f: 
    for line in f:
        if line[0] == "t":
            timeout.write( line.split()[1] )
            timeout.write("\n")
print "done writing times"
