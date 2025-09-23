
all:
	clang++ -O3 -fno-stack-protector -fomit-frame-pointer -march=native -o lnrc inner_loop.cpp proxy_loop.cpp client_loop.cpp outer_loop.cpp tlv.cpp -lpthread -lm

clean:
	rm -f ./lnrc nohup.out

run:
	nohup ./lnrc -t `grep siblings /proc/cpuinfo | sort -u | cut -d: -f 2` -s 192.168.1.103 -proxy &

