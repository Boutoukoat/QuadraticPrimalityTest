
#ifndef SETUP_H_H
#define SETUP_H_H

#define PROXY_PORT 15001
#define SERVER_PORT 15002

#define MAX_CID 2048 // a power of 2

void *server_thread(void *arg);
void *proxy_thread(void *arg);
void *client_thread(void *arg);

#endif
