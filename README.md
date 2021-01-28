CoCuckoo
========

CoCuckoo is a concurrent cuckoo hashing scheme for cloud storage systems. This repository is the implementation of CoCuckoo. For more details about CoCuckoo, please refer to [our paper](https://csyhua.github.io/csyhua/hua-USENIX-ATC2019.pdf) in USENIX ATC 2019.

> Yuanyuan Sun, Yu Hua, Zhangyu Chen, Yuncheng Guo, "Mitigating Asymmetric Read and Write Costs in Cuckoo Hashing for Storage Systems", Proceedings of the USENIX Annual Technical Conference (USENIX ATC), 2019.


Usage
-----
1. Include the header file (i.e., "CoCuckoo.h")

2. The prototype of `CoCuckoo` template is:

```c++
template <size_t NKEY, size_t NVAL, size_t POWER>
class CoCuckoo;
```

The required parameters for `CoCuckoo` include **NKEY** (number of bytes in key), **NVAL** (number of bytes in value), and **POWER** (hash table size = 2^POWER).
Typical examples of using CoCuckoo can be found in the [./test](test/) folder.
