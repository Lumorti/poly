#!/bin/bash
alpha=0.9
iters=100000
#./mub -f -d 2 -n 3 -v 2 -i ${iters} -a ${alpha} | tee data/d2n3.log
./mub -f -d 2 -n 4 -v 2 -i ${iters} -a ${alpha} | tee data/d2n4.log
#./mub -f -d 3 -n 3 -v 2 -i ${iters} -a ${alpha} | tee data/d3n3.log
#./mub -f -d 3 -n 4 -v 2 -i ${iters} -a ${alpha} | tee data/d3n4.log
./mub -f -d 3 -n 5 -v 2 -i ${iters} -a ${alpha} | tee data/d3n5.log
#./mub -f -d 4 -n 2 -v 2 -i ${iters} -a ${alpha} | tee data/d4n2.log
#./mub -f -d 4 -n 3 -v 2 -i ${iters} -a ${alpha} | tee data/d4n3.log
#./mub -f -d 4 -n 4 -v 2 -i ${iters} -a ${alpha} | tee data/d4n4.log
#./mub -f -d 4 -n 5 -v 2 -i ${iters} -a ${alpha} | tee data/d4n5.log
./mub -f -d 4 -n 6 -v 2 -i ${iters} -a ${alpha} | tee data/d4n6.log
#./mub -f -d 5 -n 2 -v 2 -i ${iters} -a ${alpha} | tee data/d5n2.log
#./mub -f -d 5 -n 3 -v 2 -i ${iters} -a ${alpha} | tee data/d5n3.log
#./mub -f -d 5 -n 4 -v 2 -i ${iters} -a ${alpha} | tee data/d5n4.log
#./mub -f -d 5 -n 5 -v 2 -i ${iters} -a ${alpha} | tee data/d5n5.log
#./mub -f -d 5 -n 6 -v 2 -i ${iters} -a ${alpha} | tee data/d5n6.log
./mub -f -d 5 -n 7 -v 2 -i ${iters} -a ${alpha} | tee data/d5n7.log
#./mub -f -d 6 -n 2 -v 2 -i ${iters} -a ${alpha} | tee data/d6n2.log
#./mub -f -d 6 -n 3 -v 2 -i ${iters} -a ${alpha} | tee data/d6n3.log
./mub -f -d 6 -n 4 -v 2 -i ${iters} -a ${alpha} | tee data/d6n4.log
grep -H "finished in" data/d2n3.log
grep -H "finished in" data/d2n4.log
grep -H "finished in" data/d3n3.log
grep -H "finished in" data/d3n4.log
grep -H "finished in" data/d3n5.log
grep -H "finished in" data/d4n2.log
grep -H "finished in" data/d4n3.log
grep -H "finished in" data/d4n4.log
grep -H "finished in" data/d4n5.log
grep -H "finished in" data/d4n6.log
grep -H "finished in" data/d5n2.log
grep -H "finished in" data/d5n3.log
grep -H "finished in" data/d5n4.log
grep -H "finished in" data/d5n5.log
grep -H "finished in" data/d5n6.log
grep -H "finished in" data/d5n7.log
grep -H "finished in" data/d6n2.log
grep -H "finished in" data/d6n3.log
grep -H "finished in" data/d6n4.log