# SPPRC-PARALLEL
Here is an attempt to accelerate SPPRC algorithm in BOOST library.

The scenario is to find the shortest path from the starting point to the endpoint in the graph. Since there is more than one label, there may be lots of non-dominated skylines. The SPPRC algorithm in the BOOST library gives me a framework, but it is not fast enough to meet my requirements. So I tried different ways to accelerate it. Ultimately, I could get a 10x acceleration compared to the original one.

## Things I tried
### Multithreading (1000%)
Parallel computing at vertex level is the most intuitive and effective, as checking the dominance on different vertices is relatively independent. But the lock must be carefully designed. Each vertex has a lock, and there is a possibility of deadlock when adding derived labels to other vertices. The solution is to design a buffer for each verbex, which has its own independent lock. When other vertices add labels to a vertex, the labels go to the buffer instead of the non buffer being checked for dominance. Then the buffer will be flushed before checking dominance next time.

### SIMD (5-10%)
Checking dominance requires multiple comparisons. One can use AVX instructions to accelerate this process. But the improvement is relatively small. 

### Custom queue with two locks of head and tail (no improvement)
I observe that 10% of time is spent waiting for the queue lock. Often, the queue is very long, and I wonder if it can become faster using fine-grain lock, which decouples producers and consumers competition for the queue lock. With the help of a virtual head node, I implemented this queue using a linked list. Unfortuantely it is much slower than std::queue on operations like pop, which results in no final speed improvement.

### Lockfree concurrent queue (no improvement)
I tried [concurrent queue](https://github.com/cameron314/concurrentqueue) to replace the std::queue with mutex, but it became slower.

## Things I didn't try
### Store labels in vector
To support fast deletion and less memory usage, pointers to labels are stored in linked lists. However, it is possible to use a vector container to store labels data directly and use a delete marker to indicate deletion. This will take the advantage of the locality principle and greatly improve the traversal speed (approx 10x). The biggest challenge is how to represent labels and ensure that labels are valid in the task queue without taking up too much memory. 

## Test
I tested the algorithm on [real traffic data](https://github.com/bstabler/TransportationNetworks/tree/master/chicago-regional) and generated random data. The results are as follows
| version | data | node | average out-degree | labels | time (s) |
| - | - | - | - | - | - |
| origin | real | 12982 | 3 | 5 | 502.802 |
| SIMD | real | 12982 | 3 | 5 | 462.585 |
| multithreading | real | 12982 | 3 | 5 | 49.0031 |
| SIMD+multithreading | real | 12982 | 3 | 5 | 43.9188 |
| origin | generated | 2000 | 10 | 30 | 180.881 |
| SIMD | generated | 2000 | 10 | 30 | 163.313 |
| multithreading | generated | 2000 | 10 | 30 | 16.3271 |
| SIMD+multithreading | generated | 2000 | 10 | 30 | 15.2957 |

