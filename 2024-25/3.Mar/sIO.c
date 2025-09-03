#include <stdio.h> 
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#define BUFSIZE 10000000 

int main(int argc, char *argv[]) {    

	int i, myrank, buf[BUFSIZE];
    	char filename[128];
    	FILE *myfile;

    	strcpy(filename, "testfile");
	myfile = fopen(filename, "w");    
        
     	for (i=0; i<BUFSIZE; i++) {  
	   buf[i] = myrank + i;
	   fprintf(myfile, "%d\n", buf[i]);
	}   

	fclose(myfile);    

	return 0; 

} 
