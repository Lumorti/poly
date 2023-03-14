#!/bin/bash
iters=100000
cores=8
#./mub -f -c ${cores} -d 2 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d2n3.log
#./mub -f -c ${cores} -d 2 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d2n4.log
#./mub -f -c ${cores} -d 3 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d3n3.log
#./mub -f -c ${cores} -d 3 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d3n4.log
#./mub -f -c ${cores} -d 3 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d3n5.log
#./mub -f -c ${cores} -d 4 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n2.log
#./mub -f -c ${cores} -d 4 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n3.log
#./mub -f -c ${cores} -d 4 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n4.log
#./mub -f -c ${cores} -d 4 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n5.log
#./mub -f -c ${cores} -d 4 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n6.log
#./mub -f -c ${cores} -d 5 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n2.log
#./mub -f -c ${cores} -d 5 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n3.log
#./mub -f -c ${cores} -d 5 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n4.log
#./mub -f -c ${cores} -d 5 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n5.log
#./mub -f -c ${cores} -d 5 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n6.log
#./mub -f -c ${cores} -d 5 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n7.log
#./mub -f -c ${cores} -d 6 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d6n2.log
#./mub -f -c ${cores} -d 6 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d6n3.log
#./mub -f -c ${cores} -d 6 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d6n4.log
#./mub -f -c ${cores} -d 10 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d10n2.log
#./mub -f -c ${cores} -d 10 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d10n3.log
#./mub -f -c ${cores} -d 10 -n 4 -v 2 -i ${iters} -a 0.3 -b 1e-13| tee data/d10n4.log
#./mub -f -c ${cores} -d 12 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d12n2.log
#./mub -f -c ${cores} -d 12 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d12n3.log
#./mub -f -c ${cores} -d 12 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d12n4.log
#./mub -f -c ${cores} -d 14 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d14n2.log
#./mub -f -c ${cores} -d 14 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d14n3.log
#./mub -f -c ${cores} -d 14 -n 4 -v 2 -i ${iters} -a 0.5 -b 1e-13| tee data/d14n4.log
#./mub -f -c ${cores} -d 15 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d15n2.log
#./mub -f -c ${cores} -d 15 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d15n3.log
#./mub -f -c ${cores} -d 15 -n 4 -v 2 -i ${iters} -a 0.3 -b 1e-13| tee data/d15n4.log
./mub -f -c ${cores} -d 16 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d16n2.log
./mub -f -c ${cores} -d 16 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d16n3.log
./mub -f -c ${cores} -d 16 -n 4 -v 2 -i ${iters} -a 0.3 -b 1e-13| tee data/d16n4.log
./mub -f -c ${cores} -d 18 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d18n2.log
./mub -f -c ${cores} -d 18 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-13| tee data/d18n3.log
./mub -f -c ${cores} -d 18 -n 4 -v 2 -i ${iters} -a 0.3 -b 1e-13| tee data/d18n4.log
#grep -H "finished in" data/d2n3.log
#grep -H "finished in" data/d2n3.log
#grep -H "finished in" data/d2n4.log
#grep -H "finished in" data/d3n3.log
#grep -H "finished in" data/d3n4.log
#grep -H "finished in" data/d3n5.log
#grep -H "finished in" data/d4n2.log
#grep -H "finished in" data/d4n3.log
#grep -H "finished in" data/d4n4.log
#grep -H "finished in" data/d4n5.log
#grep -H "finished in" data/d4n6.log
#grep -H "finished in" data/d5n2.log
#grep -H "finished in" data/d5n3.log
#grep -H "finished in" data/d5n4.log
#grep -H "finished in" data/d5n5.log
#grep -H "finished in" data/d5n6.log
#grep -H "finished in" data/d5n7.log
#grep -H "finished in" data/d6n2.log
#grep -H "finished in" data/d6n3.log
#grep -H "finished in" data/d6n4.log
#grep -H "finished in" data/d10n2.log
#grep -H "finished in" data/d10n3.log
#grep -H "finished in" data/d10n4.log
#grep -H "finished in" data/d12n2.log
#grep -H "finished in" data/d12n3.log
#grep -H "finished in" data/d12n4.log
#grep -H "finished in" data/d14n2.log
#grep -H "finished in" data/d14n3.log
#grep -H "finished in" data/d14n4.log
#grep -H "finished in" data/d15n2.log
#grep -H "finished in" data/d15n3.log
#grep -H "finished in" data/d15n4.log
grep -H "finished in" data/d16n2.log
grep -H "finished in" data/d16n3.log
grep -H "finished in" data/d16n4.log
grep -H "finished in" data/d18n2.log
grep -H "finished in" data/d18n3.log
grep -H "finished in" data/d18n4.log
