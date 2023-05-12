#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"
#include <stdio.h>

__global__ void ComputeKernel(vector3** accels,vector3 *hPos, double *mass){
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	//for (i=0;i<NUMENTITIES;i++){
	//	for (j=0;j<NUMENTITIES;j++){
	//cuPrintf("Curr Row: %d, Curr Col: %d",i,j);
	if(i<NUMENTITIES && j<NUMENTITIES){
		if (i==j) {
				FILL_VECTOR(accels[i][j],0,0,0);
			}
			else{
				vector3 distance;
				for (int k=0;k<3;k++) distance[k]=hPos[i][k]-hPos[j][k];
				double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
				double magnitude=sqrt(magnitude_sq);
				double accelmag=-1*GRAV_CONSTANT*mass[j]/magnitude_sq;
				FILL_VECTOR(accels[i][j],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
			}
	}
	//	}
	//}
}

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
void compute(){
	//make an acceleration matrix which is NUMENTITIES squared in size;
	int i,j,k;
	vector3* values=(vector3*)malloc(sizeof(vector3)*NUMENTITIES*NUMENTITIES);
	vector3** accels=(vector3**)malloc(sizeof(vector3*)*NUMENTITIES);
	//vector3* d_values;
	vector3** d_accels;
	//cudaMalloc(&d_values, sizeof(vector3)*NUMENTITIES*NUMENTITIES);
	cudaMalloc(&d_accels, sizeof(vector3*)*NUMENTITIES);

	for (i=0;i<NUMENTITIES;i++){
		accels[i]=&values[i*NUMENTITIES];
		//cudaMalloc(&accels[i], sizeof(vector3)*NUMENTITIES);
		//cudaMemcpy(accels[i], &values[i*NUMENTITIES], sizeof(vector3)*NUMENTITIES, cudaMemcpyHostToDevice);
	}
	cudaMemcpy(d_accels, accels, sizeof(vector3*)*NUMENTITIES, cudaMemcpyHostToDevice);
	//cudaMemcpy(d_values, values, NUMENTITIES*NUMENTITIES*sizeof(vector3), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_accels, accels, sizeof(vector3*)*NUMENTITIES, cudaMemcpyDeviceToHost);
	printf("Last Error: %d",cudaGetLastError());

	dim3 dimBlock(16, 16);
	dim3 dimGrid((NUMENTITIES / dimBlock.x)+1, (NUMENTITIES / dimBlock.y)+1);
	
	ComputeKernel<<<dimGrid, dimBlock>>>(d_accels,hPos,mass);
	//first compute the pairwise accelerations.  Effect is on the first argument.
	
	cudaMemcpy(accels, d_accels, sizeof(vector3*)*NUMENTITIES, cudaMemcpyDeviceToHost);
	//printf("Test: %f", accels[10][20][0]);

	/*
	cudaMemcpy(accels, d_accels, sizeof(vector3*)*NUMENTITIES, cudaMemcpyDeviceToHost);
	for(int l = 0; l< NUMENTITIES; l++){
		cudaMemcpy(values + l, accels[l], sizeof(vector3), cudaMemcpyDeviceToHost);
	}
	*/
	//printf("accles Size: %lu, d_accels Size: %lu, Correct Size: %lu",sizeof(accels),sizeof(d_accels),sizeof(vector3*)*NUMENTITIES);
	
	//sum up the rows of our matrix to get effect on each entity, then update velocity and position.
	for (i=0;i<NUMENTITIES;i++){
		vector3 accel_sum={0,0,0};
		for (j=0;j<NUMENTITIES;j++){
			for (k=0;k<3;k++)
				accel_sum[k]+=accels[i][j][k];
		}
		//compute the new velocity based on the acceleration and time interval
		//compute the new position based on the velocity and time interval
		for (k=0;k<3;k++){
			hVel[i][k]+=accel_sum[k]*INTERVAL;
			hPos[i][k]=hVel[i][k]*INTERVAL;
		}
	}

	free(accels);
	free(values);
	//cudaFree(d_values);
	cudaFree(d_accels);
}
