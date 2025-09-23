# Quadratic primality test (multithreaded)

!!! WORK IN PROGRESS !!!


# Simple Test

```
$ g++ -O3 -o quadratic quadratic.cpp -lgmp
$ echo 17000000000000000000000000000000000007 > testfile.txt
$ ./quadratic testfile.txt
  Has factors!
$ echo 170000000000000000000000000000000000000000000003 > testfile.txt
$ ./quadratic testfile.txt
  Is prime!
```



