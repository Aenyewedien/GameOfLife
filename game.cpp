# include <iostream>
#include <sstream>
#include "bitmap_image.hpp"
# include <mpi.h>

using namespace std;
int main(int argc, char *argv[])

{
 	int iterations= atoi(argv[1]);
	int it=0;
	int wym = atoi(argv[2]);
	int size, rank; // size - ilosc procesow
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat;
	int row_per_proc = wym / size;
	int per_proc = wym * row_per_proc;

		int* table = new int[wym*wym];



	// inicjalizacja "T"

	if (rank == 0) {	
		for (int i = 0; i < wym; ++i) {
			for (int j = 0; j < wym; ++j) {
				if (i < wym / 3 || (j >= wym / 3 && j < 2 * wym / 3)) {
					table[i+ j*wym] = 1;

				} else {

					table[i+ j*wym] = 0;

				}
			//	cout<<table[i+j*wym];
			}
//cout<< endl;
		}

	}

	

	//wycinek dla jednego procesu <-- wektor
		int* smallslice = new int[row_per_proc*wym]; 
		int** pizzaslice = new int*[wym];
			for (int i = 0; i < wym; i++) {
				pizzaslice[i] = new int[row_per_proc];
			}
while(it<iterations){
		MPI_Scatter(table, per_proc, MPI_INT, smallslice, per_proc, MPI_INT, 0, MPI_COMM_WORLD);

for(int i =0;i<wym;i++){
	for (int j=0;j<row_per_proc;j++){
			pizzaslice[i][j]=smallslice[i+j*wym];
		//	cout<<pizzaslice[i][j];
		}
	//	cout<<endl;
	}

	int* getright = new int[wym];
	int* getleft = new int[wym];
	int* sendright = new int[wym];
	int* sendleft = new int[wym];
	if (rank != size - 1) // all except for last send down
		{
			for (int j = 0; j < wym; j++)
			{
				sendright[j] = pizzaslice[j][row_per_proc-1];
			}


MPI_Send(sendright, wym, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);			
		}
	else {
			for (int k = 0; k < wym; k++)
			{
				getright[k] = 0;
			}
		}
			if (rank != 0)
		{
			MPI_Recv(getleft, wym, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, &stat);
		}
		else {
			for (int k = 0; k < wym; k++)
			{
			getleft[k] = 0;
			}
		}
		if (rank != 0)
		{
			for (int j = 0; j < wym; j++)
			{
				sendleft[j] = pizzaslice[j][0];			
			}
			MPI_Send(sendleft, wym, MPI_INT, rank - 1, 1, MPI_COMM_WORLD);
		}
		
			if (rank != size - 1)
			{
			MPI_Recv(getright, wym, MPI_INT, rank + 1, 1, MPI_COMM_WORLD, &stat);
			}

	int sum = 0; // sum of neighbours
//TUTAJ OBLICZANA JEST SUMA TABLICY;
/*if(rank==0){
for(int y=0; y<wym;y++){
	for (int x = 0;x<row_per_proc;x++){
cout<<pizzaslice[y][x];
}
cout<< endl;
}}
*/
int** lastslice = new int*[wym];
for(int i = 0;i<wym;i++){
lastslice[i] = new int[row_per_proc];
}

// tutaj
	#pragma omp parallel for
for(int y= 0;y<wym;y++){
for(int x = 0; x<row_per_proc;x++){
	
		if(x==0 && y==0){
	sum = pizzaslice[y][x+1]+ pizzaslice[y+1][x+1]+ pizzaslice[y+1][x]+getleft[y]+getleft[y+1];
} else if(x==0 && y == wym-1){
	sum = pizzaslice[y-1][x]+pizzaslice[y-1][x+1] + pizzaslice[y][x+1]+ getleft[y]+getleft[y-1]; 
}else if(x==row_per_proc-1 && y==wym-1){
	sum = pizzaslice[y-1][x-1]+pizzaslice[y-1][x]+ pizzaslice[y][x-1] + getright[y]+getright[y-1];
	}else if(x==row_per_proc-1 && y==0){
	sum = pizzaslice[y+1][x-1]+pizzaslice[y][x-1] + pizzaslice[y+1][x]+ getright[y]+getright[y+1];
}
	else{
	if (y==0){
		sum= pizzaslice[y][x-1]+pizzaslice[y][x+1]+pizzaslice[y+1][x+1]+pizzaslice[y+1][x]+pizzaslice[y+1][x-1]; 
		}else if(y==wym-1){
	sum= pizzaslice[y-1][x-1]+pizzaslice[y-1][x]+ pizzaslice[y-1][x+1]+ pizzaslice[y][x-1]+ pizzaslice[y][x+1];
	}else if(x==row_per_proc-1){
	sum = pizzaslice[y-1][x-1]+ pizzaslice[y-1][x]+pizzaslice[y][x-1]+ pizzaslice[y+1][x-1]+ pizzaslice[y+1][x-1]+ getright[y]+getright[y-1]+getright[y+1];
	}else if(x==0){
	sum = pizzaslice[y-1][x]+pizzaslice[y-1][x+1]+pizzaslice[y][x+1]+pizzaslice[y+1][x]+ pizzaslice[y+1][x+1]+ getleft[y]+ getleft[y-1]+getleft[y+1];
	 } else{sum=pizzaslice[y-1][x-1]+pizzaslice[y-1][x]+pizzaslice[y-1][x+1]+pizzaslice[y][x-1]+pizzaslice[y][x+1]+pizzaslice[y+1][x-1]+pizzaslice[y+1][x]+pizzaslice[y+1][x+1];}
}
				if (pizzaslice[y][x] == 1 && (sum == 2 || sum == 3))
				{
					lastslice[y][x] = 1;
				}
				else if (pizzaslice[y][x] == 1 && sum > 3)
				{
					lastslice[y][x] = 0;;
				}
				else if (pizzaslice[y][x] == 1 && sum < 1)
				{
					lastslice[y][x] = 0;
				}
				else if (pizzaslice[y][x] == 0 && sum == 3)
				{
					lastslice[y][x] = 1;
				}
				else
				{
					lastslice[y][x] = 0;
				}
}
}

for(int y= 0;y<wym;y++){
for(int x = 0; x<row_per_proc;x++){
smallslice[y+x*wym]= lastslice[y][x];
}
}
MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(smallslice, per_proc, MPI_INT, table, per_proc, MPI_INT, 0, MPI_COMM_WORLD); 

	if (rank == 0){
bitmap_image something(wym, wym);
something.clear();
		for (int i = 0; i < wym; ++i) {
			for (int j = 0; j < wym; ++j) {
			if (table[i+j*wym] == 1) {
				rgb_t white = make_colour(255, 255, 255);
				something.set_pixel(j, i, white);
			}else{
				rgb_t black = make_colour(0, 0, 0);
				something.set_pixel(j, i, black);
			}
		}
		}
	
				stringstream ss;
				ss << it;//<===== tu zmienna od iteracji;
				string i = ss.str();
				string filename = "file_" + i + ".bmp";
				const char* c = filename.c_str();
				something.save_image(c);
}
it++;
free(getleft);
free(getright);
free(sendleft);
free(sendright);
}
	free(smallslice);
	free(table);
	free(pizzaslice);

	MPI_Finalize();

	return 0;

}
