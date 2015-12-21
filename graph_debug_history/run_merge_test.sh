#!/bin/bash

set -e

#./metagraph -v -k 3 -O test_out test.fa
#./metagraph -v -k 3 -O test4_out test4.fa
#./metagraph -v -k 3 -O test5_out test5.fa
#./metagraph -v -k 3 -O test6_out test6.fa
#./metagraph -v -k 3 -O test7_out test7.fa
#./metagraph -v -k 3 -O test8_out test8.fa
#./metagraph -v -k 3 -O test9_out test9.fa
#./metagraph -v -k 3 -O test10_out test10.fa
#./metagraph -v -k 3 -O test11_out test11.fa
#./metagraph -v -k 3 -O test12_out test12.fa
#./metagraph -v -k 3 -O test13_out test13.fa
#./metagraph -v -k 3 -O test14_out test14.fa
#./metagraph -v -k 3 -O test15_out test15.fa
#./metagraph -v -k 3 -O test16_out test16.fa

#./metagraph -v -m test_out,test4_out,test5_out -O test_merge_out BBB
#./metagraph -v -k 3 -O test_merge2_out test.fa test4.fa test5.fa

#./metagraph -v -m test_out,test4_out,test5_out,test6_out,test7_out,test8_out,test9_out,test10_out,test11_out,test12_out,test13_out,test14_out,test15_out,test16_out -O test_merge_out BBB
./metagraph -v -m test_out,test6_out -O test_merge_out BBB
md5sum test_merge_out.W.dbg
md5sum test_merge_out.l.dbg

#./metagraph -v -m test6_out,test_out -O test_merge4_out BBB
#md5sum test_merge4_out.W.dbg
#md5sum test_merge4_out.l.dbg

#./metagraph -v -k 3 -O test_merge2_out test.fa test4.fa test5.fa test6.fa test7.fa test8.fa test9.fa test10.fa test11.fa test12.fa test13.fa test14.fa test15.fa test16.fa
./metagraph -v -k 3 -O test_merge2_out test.fa test6.fa
md5sum test_merge2_out.W.dbg
md5sum test_merge2_out.l.dbg

./metagraph -v -c test_merge_out,test_merge2_out BBB

#./metagraph -v -k 3 -O test_merge3_out test6.fa test.fa
#md5sum test_merge3_out.W.dbg
#md5sum test_merge3_out.l.dbg
