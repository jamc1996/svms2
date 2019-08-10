#include <stdio.h>

int main(){

	int A[2][2];
	printf("A = %p\n",A);
	printf("A[0] = %p\n",A[0]);
	printf("&A[0][0] = %p\n", &(A[0][0]));
	return 0;
}
