#!/bin/bash

# Variables
cores=4
iters=100000

#./mub -f -c ${cores} -d 3 -N 3,3,3,3   -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N3333.log
#./mub -f -c ${cores} -d 3 -N 1,1,1,1,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N11111.log
#./mub -f -c ${cores} -d 3 -N 2,1,1,1,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N21111.log
#./mub -f -c ${cores} -d 3 -N 2,2,1,1,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N22111.log
#./mub -f -c ${cores} -d 3 -N 3,1,1,1,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N31111.log
#./mub -f -c ${cores} -d 3 -N 2,2,2,1,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N22211.log
#./mub -f -c ${cores} -d 3 -N 3,2,1,1,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N32111.log
#./mub -f -c ${cores} -d 3 -N 2,2,2,2,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N22221.log
#./mub -f -c ${cores} -d 3 -N 3,2,2,1,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N32211.log
#./mub -f -c ${cores} -d 3 -N 3,3,1,1,1 -v 2 -i ${iters} -a 0.95 -b 1e-10 | tee data/d3N33111.log

#./mub -f -c ${cores} -d 4 -N 4,4,4,4,4   -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N44444.log
#./mub -f -c ${cores} -d 4 -N 1,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N111111.log
#./mub -f -c ${cores} -d 4 -N 2,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N211111.log
#./mub -f -c ${cores} -d 4 -N 2,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N221111.log
#./mub -f -c ${cores} -d 4 -N 3,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N311111.log
#./mub -f -c ${cores} -d 4 -N 2,2,2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N222111.log
#./mub -f -c ${cores} -d 4 -N 3,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N321111.log
#./mub -f -c ${cores} -d 4 -N 4,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N411111.log
#./mub -f -c ${cores} -d 4 -N 2,2,2,2,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N222211.log
#./mub -f -c ${cores} -d 4 -N 3,2,2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N322111.log
#./mub -f -c ${cores} -d 4 -N 4,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N421111.log
#./mub -f -c ${cores} -d 4 -N 3,3,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N331111.log
#./mub -f -c ${cores} -d 4 -N 2,2,2,2,2,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N222221.log
#./mub -f -c ${cores} -d 4 -N 3,2,2,2,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N322211.log
#./mub -f -c ${cores} -d 4 -N 3,3,2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N332111.log
#./mub -f -c ${cores} -d 4 -N 4,2,2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N422111.log
#./mub -f -c ${cores} -d 4 -N 4,3,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d4N431111.log

#./mub -f -c ${cores} -d 5 -N 5,5,5,5,5,5   -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N555555.log
#./mub -f -c ${cores} -d 5 -N 1,1,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N1111111.log
#./mub -f -c ${cores} -d 5 -N 2,1,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N2111111.log
#./mub -f -c ${cores} -d 5 -N 2,2,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N2211111.log
#./mub -f -c ${cores} -d 5 -N 3,1,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3111111.log
#./mub -f -c ${cores} -d 5 -N 2,2,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N2221111.log
#./mub -f -c ${cores} -d 5 -N 3,2,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3211111.log
#./mub -f -c ${cores} -d 5 -N 4,1,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N4111111.log
#./mub -f -c ${cores} -d 5 -N 2,2,2,2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N2222111.log
#./mub -f -c ${cores} -d 5 -N 3,2,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3221111.log
#./mub -f -c ${cores} -d 5 -N 3,3,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3311111.log
#./mub -f -c ${cores} -d 5 -N 4,2,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N4211111.log
#./mub -f -c ${cores} -d 5 -N 5,1,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N5111111.log
#./mub -f -c ${cores} -d 5 -N 2,2,2,2,2,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N2222211.log
#./mub -f -c ${cores} -d 5 -N 3,2,2,2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3222111.log
#./mub -f -c ${cores} -d 5 -N 3,3,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3321111.log
#./mub -f -c ${cores} -d 5 -N 4,2,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N4221111.log
#./mub -f -c ${cores} -d 5 -N 5,2,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N5211111.log
#./mub -f -c ${cores} -d 5 -N 2,2,2,2,2,2,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N2222221.log
#./mub -f -c ${cores} -d 5 -N 3,2,2,2,2,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3222211.log
#./mub -f -c ${cores} -d 5 -N 4,2,2,2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N4222111.log
#./mub -f -c ${cores} -d 5 -N 3,3,2,2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3322111.log
#./mub -f -c ${cores} -d 5 -N 3,3,3,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N3331111.log
#./mub -f -c ${cores} -d 5 -N 4,3,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N4321111.log
#./mub -f -c ${cores} -d 5 -N 5,2,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N5221111.log
#./mub -f -c ${cores} -d 5 -N 4,4,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N4411111.log
#./mub -f -c ${cores} -d 5 -N 5,3,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d5N5311111.log

#./mub -f -c ${cores} -d 6 -N 6,6,6   -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N666.log
#./mub -f -c ${cores} -d 6 -N 1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N1111.log
#./mub -f -c ${cores} -d 6 -N 2,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N2111.log
#./mub -f -c ${cores} -d 6 -N 2,2,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N2211.log
#./mub -f -c ${cores} -d 6 -N 3,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N3111.log
#./mub -f -c ${cores} -d 6 -N 2,2,2,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N2221.log
#./mub -f -c ${cores} -d 6 -N 3,2,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N3211.log
#./mub -f -c ${cores} -d 6 -N 4,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N4111.log
#./mub -f -c ${cores} -d 6 -N 2,2,2,2 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N2222.log
#./mub -f -c ${cores} -d 6 -N 3,2,2,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N3221.log
#./mub -f -c ${cores} -d 6 -N 3,3,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N3311.log
#./mub -f -c ${cores} -d 6 -N 4,2,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N4211.log
#./mub -f -c ${cores} -d 6 -N 3,2,2,2 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N3222.log
#./mub -f -c ${cores} -d 6 -N 6,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6111.log
#./mub -f -c ${cores} -d 6 -N 3,3,2,2 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N3322.log
#./mub -f -c ${cores} -d 6 -N 6,2,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6211.log
#./mub -f -c ${cores} -d 6 -N 3,3,3,2 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N3332.log
#./mub -f -c ${cores} -d 6 -N 6,3,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6311.log
#./mub -f -c ${cores} -d 6 -N 3,3,3,3 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N3333.log
#./mub -f -c ${cores} -d 6 -N 6,4,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6411.log
#./mub -f -c ${cores} -d 6 -N 4,3,3,3 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N4333.log
#./mub -f -c ${cores} -d 6 -N 6,5,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6511.log
#./mub -f -c ${cores} -d 6 -N 4,4,3,3 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N4433.log
#./mub -f -c ${cores} -d 6 -N 6,6,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6611.log
#./mub -f -c ${cores} -d 6 -N 4,4,4,3 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N4443.log
#./mub -f -c ${cores} -d 6 -N 6,6,2,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6621.log
#./mub -f -c ${cores} -d 6 -N 4,4,4,4 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N4444.log
#./mub -f -c ${cores} -d 6 -N 6,6,3,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6631.log
#./mub -f -c ${cores} -d 6 -N 6,5,2,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6521.log
#./mub -f -c ${cores} -d 6 -N 6,5,3,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6531.log
#./mub -f -c ${cores} -d 6 -N 6,5,4,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N6541.log
#./mub -f -c ${cores} -d 6 -N 6,2,2,1,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N62211111.log
#./mub -f -c ${cores} -d 6 -N 6,2,2,2,1,1,1,1 -v 2 -i ${iters} -a 0.9 -b 1e-10 | tee data/d6N62221111.log

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
#./mub -2 -f -c ${cores} -d 11 -n 13 -v 2 -i ${iters} -a 0.6 -b 1e-15 | tee data/d11n13.log

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
for file in $(ls data/d6N*.log | sort); do
    if [ -f "$file" ]; then
        grep -H "vars" "$file"
        grep -H "iteration" "$file"
    fi
done

