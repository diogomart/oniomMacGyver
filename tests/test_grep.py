import c_grep
out = c_grep.c_grep("big_scan_oniom.log", ['Z-mat', "Optimized P", "Converged?", "scan point","SCF Done"])
#print(out)
#print(type(out))

for line in out:
    print(line)

