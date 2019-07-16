#include <stdio.h>
#include <stdlib.h>


int main(){
  FILE *file = fopen("bigData.txt", "w");



  double num;
  for (int i = 0; i < 10000; i++) {
    for (int i = 0; i < 40; i++) {
      fprintf(file, "%lf\t",num );
      num = (drand48()*20)-19;
    }
    fprintf(file, "%lf\n",-1.0 );

  }

  for (int i = 0; i < 10000; i++) {
    for (int i = 0; i < 40; i++) {
      fprintf(file, "%lf\t",num );
      num = 15+(drand48()*20);
    }
    fprintf(file, "%lf\n",1.0 );
  }

  for (int i = 0; i < 40; i++) {
	for (int j= 0; j < 40; j++){
		if(i==j){
			fprintf(file, "4\t");
		}
		else{
			fprintf(file, "1\t");

		} 
 	}
  	fprintf(file, "%lf\n",-1.0 );
  }

  for (int i = 0; i < 40; i++) {
	for (int j= 0; j < 40; j++){
		if(i==j){
			fprintf(file, "11\t");
		}
		else{
			fprintf(file, "21\t");

		} 
 	}
  	fprintf(file, "%lf\n",1.0 );
  }



  fclose(file);
  return 0;
}
