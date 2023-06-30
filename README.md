# variable-encoding-framework
We propose a novel variable-length encoding framework, which can be applied to AMQ data structures to improve the space efficiency and reduce the false positive rate.

## Dependencies
* `OpenSSL`
  
## Build and run
```sh
mkdir build
cd build
cmake ..
make
./test/1_throughput
```

## Evaluation
|Algorithm| Description|
|:----:|----|
|Cuckoo Filter (CF)|B. Fan, D. G. Andersen, M. Kaminsky, and M. Mitzenmacher, “Cuckoo filter: Practically better than bloom,” in Proceedings of ACM International Conference on Emerging Networking Experiments and Technologies. ACM, 2014, pp. 75–88. Implementation: https://github.com/efficient/cuckoofilter|
|Configurable Bucket Cuckoo Filter (CBCF)|P. Reviriego, J. Mart´ınez, D. Larrabeiti, and S. Pontarelli, “Cuckoo filters and bloom filters: Comparison and application to packet classification,” IEEE Transactions on Network and Service Management, vol. 17, no. 4, pp. 2690–2701, 2020. Implementation: https://github.com/jmladron/cbcf|
|Tagged Cuckoo Filter (TCF)|K. Huang and T. Yang, “Tagged cuckoo filters,” in Proceedings of International Conference on Computer Communications and Networks. IEEE, 2021, pp. 1–10.|
|Quotient Filter (QF)|M. A. Bender, M. Farach-Colton, R. Johnson, R. Kraner, B. C. Kuszmaul, D. Medjedovic, P. Montes, P. Shetty, R. P. Spillane, and E. Zadok, “Don’t thrash: How to cache your hash on flash,” Proceedings of the VLDB Endowment, vol. 5, no. 11, pp. 1627–1637, 2012. Implementation: https://github.com/vedantk/quotient-filter |
|Vector Quotient Filter (VQF)|P. Pandey, A. Conway, J. Durie, M. A. Bender, M. Farach-Colton, and R. Johnson, “Vector quotient filters: Overcoming the time/space trade-off in filter design,” in Proceedings of International Conference on Management of Data. ACM, 2021, pp. 1386–1399. Implementation: https://github.com/splatlab/vqf|
|Counting Bloom Filter (CBF)|L. Fan, P. Cao, J. M. Almeida, and A. Z. Broder, “Summary cache: A scalable wide-area web cache sharing protocol,” IEEE/ACM Transactions on Networking, vol. 8, no. 3, pp. 281–293, 2000.|
|Variable-increment Counting Bloom Filter (VCBF)|O. Rottenstreich, Y. Kanizo, and I. Keslassy, “The variable-increment counting bloom filter,” IEEE/ACM Transactions on Networking, vol. 22, no. 4, pp. 1092–1105, 2014.|
|Tandem Counting Bloom Filter (TCBF)|P. Reviriego and O. Rottenstreich, “The tandem counting bloom filter - it takes two counters to tango,” IEEE/ACM Transactions on Networking, vol. 27, no. 6, pp. 2252–2265, 2019.|


## To generate YCSB workloads
```sh
cd YCSB
wget https://github.com/brianfrankcooper/YCSB/releases/download/0.17.0/ycsb-0.17.0.tar.gz
tar -xvf ycsb-0.17.0.tar.gz
mv ycsb-0.17.0 YCSB
#Then run workload generator
mkdir workloads
./generate_all_workloads.sh
```
