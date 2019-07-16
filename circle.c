#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
  FILE *file = fopen("circleData.txt", "w");



  double num, num2, num3;
  for (int i = 0; i < 10; i++) {
      	num = (drand48()*2)-1;
      	fprintf(file, "%lf\t",num );
      	num2 = sqrt(drand48()*(1-(num*num)));
     	num3 = drand48();
	if(num3 < 0.5){
	      fprintf(file, "%lf\t",-num2 );
	}
	else
	{
	      fprintf(file, "%lf\t",num2 );

	}
      	fprintf(file, "%lf\n",-1.0 );
  }

  for (int i = 0; i < 10; i++) {
    	num = (drand48()*4)-2;
    	fprintf(file, "%lf\t",num );
	if(num*num < 1){
		num2 = (drand48())+1.0;
		num3 = drand48();
		if(num3 < 0.5){
			fprintf(file, "%lf\t",-num2 );
		}else{
			fprintf(file, "%lf\t",num2 );
		}
	}
	else{
		printf("this uno\n");
		num2 = (drand48()*2)-1;
		fprintf(file, "%lf\t",num2 );
	}
    fprintf(file, "%lf\n",1.0 );
  }

  fclose(file);
}

