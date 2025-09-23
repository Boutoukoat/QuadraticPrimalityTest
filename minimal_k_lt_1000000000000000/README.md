# Quadratic test verification code

WORK IN PROGRESS !!!

Test code to verify a quadratic primality test, based on linear recurrences

So far, no counterexample (false positive, false negative) has been found. Exhaustive tests completed to 10^15

# Simple Test user's guide

```
$ make

$ ./lnrc -st
  run the self tests 

$ ./lnrc -server -e 1000000
  run a server, to dispatch the work of verifying the primality of numbers starting at 1000000.
  (check for false positives and false negatives)

$ ./lnrc -t 12 -s 127.0.0.1
  run 12 worker threads, connecting to the server, get about 2 minutes of work, and report completion, and again, and again, and again .....

```

# Distributed Test user's guide

on host 192.168.1.2

```
$ ./lnrc -server -e 1000000 | tee -a mytest.log
  run a server, to dispatch the work of verifying the primality of numbers starting at 1000000.
  (check for false positives and false negatives)
```

on each host 192.168.1.3 192.168.1.4 ...

```
$ ./lnrc -t 12 -proxy -s 192.168.1.2 -proxy
  run 12 worker threads, connecting to a proxy thread, and the proxy connects to the server, each worker thread get about 2 minutes of work, 
  and report completion, and again, and again, and again .....
```

Linux networking is not at ease with too many connections (there is a limit per process at 1024 in libc fd_sets, select() ....). 
Many worker threads can connect to a proxy, and many proxies can connect to the server. This raises the maximum number of
worker threads to theoretically 1024x1024 without changing system-wide parameters such as ulimit.


The maximum recommended number of worker threads is about the number of logical cores and can be computed with the command
```
$ grep siblings /proc/cpuinfo | sort -u | cut -d: -f 2 
```

# Where is the useful code ?

The useful code is in inner_loop.cpp, where people with a little experiem=nce in the domain will recognize functions 
like gcd(), is_prime(), mod_inverse(), mod_pow(), jacobi_symbol(), kronecker_symbol(), 
is_perfect_square() ... and similar utilities common to many tools running primality tests.




