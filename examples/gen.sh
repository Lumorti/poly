#!/bin/bash

# Variables
cores=4
iters=100000

#./mub -f -c ${cores} -d 2 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d2n3.log
#./mub -f -c ${cores} -d 2 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d2n4.log

#./mub -f -c ${cores} -d 3 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d3n2.log
#./mub -f -c ${cores} -d 3 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d3n3.log
#./mub -f -c ${cores} -d 3 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d3n4.log
#./mub -f -c ${cores} -d 3 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d3n5.log

#./mub -f -c ${cores} -d 4 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n2.log
#./mub -f -c ${cores} -d 4 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n3.log
#./mub -f -c ${cores} -d 4 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n4.log
#./mub -f -c ${cores} -d 4 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d4n5.log

#./mub -f -c ${cores} -d 5 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n2.log
#./mub -f -c ${cores} -d 5 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n3.log
#./mub -f -c ${cores} -d 5 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n4.log
#./mub -f -c ${cores} -d 5 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n5.log
#./mub -f -c ${cores} -d 5 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n6.log
#./mub -f -c ${cores} -d 5 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d5n7.log

#./mub -f -c ${cores} -d 6 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d6n2.log
#./mub -f -c ${cores} -d 6 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d6n3.log
#./mub -f -c ${cores} -d 6 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d6n4.log

#./mub -f -c ${cores} -d 7 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d7n2.log
#./mub -f -c ${cores} -d 7 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d7n3.log
#./mub -f -c ${cores} -d 7 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d7n4.log
#./mub -f -c ${cores} -d 7 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d7n5.log
#./mub -f -c ${cores} -d 7 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d7n6.log
#./mub -f -c ${cores} -d 7 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d7n7.log
#./mub -f -c ${cores} -d 7 -n 8 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d7n8.log
#./mub -f -c ${cores} -d 7 -n 9 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d7n9.log

#./mub -f -c ${cores} -d 8 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n2.log
#./mub -f -c ${cores} -d 8 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n3.log
#./mub -f -c ${cores} -d 8 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n4.log
#./mub -f -c ${cores} -d 8 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n5.log
#./mub -f -c ${cores} -d 8 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n6.log
#./mub -f -c ${cores} -d 8 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n7.log
#./mub -f -c ${cores} -d 8 -n 8 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n8.log
#./mub -f -c ${cores} -2 -d 8 -n 9 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n9.log
#./mub -f -c ${cores} -2 -d 8 -n 10 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d8n10.log

#./mub -f -c ${cores} -d 9 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n2.log
#./mub -f -c ${cores} -d 9 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n3.log
#./mub -f -c ${cores} -d 9 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n4.log
#./mub -f -c ${cores} -d 9 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n5.log
#./mub -f -c ${cores} -d 9 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n6.log
#./mub -f -c ${cores} -d 9 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n7.log
#./mub -f -c ${cores} -d 9 -n 8 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n8.log
#./mub -f -2 -c ${cores} -d 9 -n 9 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n9.log
#./mub -f -2 -c ${cores} -d 9 -n 10 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n10.log
#./mub -f -2 -c ${cores} -d 9 -n 11 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d9n11.log

#./mub -f -c ${cores} -d 10 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d10n2.log
#./mub -f -c ${cores} -d 10 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d10n3.log
#./mub -f -c ${cores} -d 10 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d10n4.log

#./mub -f -c ${cores} -d 11 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d11n2.log
#./mub -f -c ${cores} -d 11 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d11n3.log
#./mub -f -c ${cores} -d 11 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d11n4.log
#./mub -f -c ${cores} -d 11 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d11n5.log
#./mub -f -c ${cores} -d 11 -n 6 -v 2 -i ${iters} -a 0.7 -b 1e-15 | tee data/d11n6.log
#./mub -f -c ${cores} -d 11 -n 7 -v 2 -i ${iters} -a 0.7 -b 1e-15 | tee data/d11n7.log
#./mub -f -c ${cores} -d 11 -n 8 -v 2 -i ${iters} -a 0.7 -b 1e-15 | tee data/d11n8.log
#./mub -2 -f -c ${cores} -d 11 -n 9 -v 2 -i ${iters} -a 0.7 -b 1e-15 | tee data/d11n9.log
#./mub -2 -f -c ${cores} -d 11 -n 10 -v 2 -i ${iters} -a 0.7 -b 1e-15 | tee data/d11n10.log
#./mub -2 -f -c ${cores} -d 11 -n 11 -v 2 -i ${iters} -a 0.7 -b 1e-15 | tee data/d11n11.log
#./mub -2 -f -c ${cores} -d 11 -n 12 -v 2 -i ${iters} -a 0.6 -b 1e-15 | tee data/d11n12.log
./mub -2 -f -c ${cores} -d 11 -n 13 -v 2 -i ${iters} -a 0.6 -b 1e-15 | tee data/d11n13.log

#./mub -f -c ${cores} -d 12 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d12n2.log
#./mub -f -c ${cores} -d 12 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d12n3.log
#./mub -f -c ${cores} -d 12 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d12n4.log

#./mub -f -c ${cores} -d 13 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n2.log
#./mub -f -c ${cores} -d 13 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n3.log
#./mub -f -c ${cores} -d 13 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n4.log
#./mub -f -c ${cores} -d 13 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n5.log
#./mub -f -c ${cores} -d 13 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n6.log
#./mub -f -c ${cores} -d 13 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n7.log
#./mub -f -c ${cores} -d 13 -n 8 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n8.log
#./mub -f -c ${cores} -d 13 -n 9 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n9.log
#./mub -f -c ${cores} -d 13 -n 10 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n10.log
#./mub -f -c ${cores} -d 13 -n 11 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n11.log
#./mub -f -c ${cores} -d 13 -n 12 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n12.log
#./mub -f -c ${cores} -d 13 -n 13 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n13.log
#./mub -f -c ${cores} -d 13 -n 14 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n14.log
#./mub -f -c ${cores} -d 13 -n 15 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d13n15.log

#./mub -f -c ${cores} -d 14 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n2.log
#./mub -f -c ${cores} -d 14 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n3.log
#./mub -f -c ${cores} -d 14 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n4.log
#./mub -f -c ${cores} -d 14 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n5.log
#./mub -f -c ${cores} -d 14 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n6.log
#./mub -f -c ${cores} -d 14 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n7.log
#./mub -f -c ${cores} -d 14 -n 8 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n8.log
#./mub -f -c ${cores} -d 14 -n 9 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n9.log
#./mub -f -c ${cores} -d 14 -n 10 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n10.log
#./mub -f -c ${cores} -d 14 -n 11 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n11.log
#./mub -f -c ${cores} -d 14 -n 12 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n12.log
#./mub -f -c ${cores} -d 14 -n 13 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n13.log
#./mub -f -c ${cores} -d 14 -n 14 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n14.log
#./mub -f -c ${cores} -d 14 -n 15 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n15.log
#./mub -f -c ${cores} -d 14 -n 16 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d14n16.log

#./mub -f -c ${cores} -d 15 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d15n2.log
#./mub -f -c ${cores} -d 15 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d15n3.log
#./mub -f -c ${cores} -d 15 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d15n4.log

#./mub -f -c ${cores} -d 16 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n2.log
#./mub -f -c ${cores} -d 16 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n3.log
#./mub -f -c ${cores} -d 16 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n4.log
#./mub -f -c ${cores} -d 16 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n5.log
#./mub -f -c ${cores} -d 16 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n6.log
#./mub -f -c ${cores} -d 16 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n7.log
#./mub -f -c ${cores} -d 16 -n 8 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n8.log
#./mub -f -c ${cores} -d 16 -n 9 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n9.log
#./mub -f -c ${cores} -d 16 -n 10 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n10.log
#./mub -f -c ${cores} -d 16 -n 11 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n11.log
#./mub -f -c ${cores} -d 16 -n 12 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n12.log
#./mub -f -c ${cores} -d 16 -n 13 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n13.log
#./mub -f -c ${cores} -d 16 -n 14 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n14.log
#./mub -f -c ${cores} -d 16 -n 15 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n15.log
#./mub -f -c ${cores} -d 16 -n 16 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n16.log
#./mub -f -c ${cores} -d 16 -n 17 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n17.log
#./mub -f -c ${cores} -d 16 -n 18 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d16n18.log

#./mub -f -c ${cores} -d 17 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n2.log
#./mub -f -c ${cores} -d 17 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n3.log
#./mub -f -c ${cores} -d 17 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n4.log
#./mub -f -c ${cores} -d 17 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n5.log
#./mub -f -c ${cores} -d 17 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n6.log
#./mub -f -c ${cores} -d 17 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n7.log
#./mub -f -c ${cores} -d 17 -n 8 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n8.log
#./mub -f -c ${cores} -d 17 -n 9 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n9.log
#./mub -f -c ${cores} -d 17 -n 10 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n10.log
#./mub -f -c ${cores} -d 17 -n 11 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n11.log
#./mub -f -c ${cores} -d 17 -n 12 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n12.log
#./mub -f -c ${cores} -d 17 -n 13 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n13.log
#./mub -f -c ${cores} -d 17 -n 14 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n14.log
#./mub -f -c ${cores} -d 17 -n 15 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n15.log
#./mub -f -c ${cores} -d 17 -n 16 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n16.log
#./mub -f -c ${cores} -d 17 -n 17 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n17.log
#./mub -f -c ${cores} -d 17 -n 18 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n18.log
#./mub -f -c ${cores} -d 17 -n 19 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d17n19.log

#./mub -f -c ${cores} -d 18 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d18n2.log
#./mub -f -c ${cores} -d 18 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d18n3.log
#./mub -f -c ${cores} -d 18 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d18n4.log

#./mub -f -c ${cores} -d 19 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n2.log
#./mub -f -c ${cores} -d 19 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n3.log
#./mub -f -c ${cores} -d 19 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n4.log
#./mub -f -c ${cores} -d 19 -n 5 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n5.log
#./mub -f -c ${cores} -d 19 -n 6 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n6.log
#./mub -f -c ${cores} -d 19 -n 7 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n7.log
#./mub -f -c ${cores} -d 19 -n 8 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n8.log
#./mub -f -c ${cores} -d 19 -n 9 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n9.log
#./mub -f -c ${cores} -d 19 -n 10 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n10.log
#./mub -f -c ${cores} -d 19 -n 11 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n11.log
#./mub -f -c ${cores} -d 19 -n 12 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n12.log
#./mub -f -c ${cores} -d 19 -n 13 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n13.log
#./mub -f -c ${cores} -d 19 -n 14 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n14.log
#./mub -f -c ${cores} -d 19 -n 15 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n15.log
#./mub -f -c ${cores} -d 19 -n 16 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n16.log
#./mub -f -c ${cores} -d 19 -n 17 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n17.log
#./mub -f -c ${cores} -d 19 -n 18 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n18.log
#./mub -f -c ${cores} -d 19 -n 19 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n19.log
#./mub -f -c ${cores} -d 19 -n 20 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n20.log
#./mub -f -c ${cores} -d 19 -n 21 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d19n21.log

#./mub -f -c ${cores} -d 20 -n 2 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d20n2.log
#./mub -f -c ${cores} -d 20 -n 3 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d20n3.log
#./mub -f -c ${cores} -d 20 -n 4 -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee data/d20n4.log

# Read all log files in the data directory
for file in $(ls data/d11*.log | sort); do
    if [ -f "$file" ]; then
        grep -H "iteration" "$file"
    fi
done

