/*
 * Voronoi implementation
 *
 * Code author: Botond Molnár
 *
 * Based on:
 * Molnár, B., Márton, I.B., Horvát, S. et al.
 * "Community detection in directed weighted networks using Voronoi partitioning"
 * Scientific Reports 14, 8124 (2024)
 * https://doi.org/10.1038/s41598-024-58624-4
 *
 * License: Creative Commons Attribution 4.0 International License (CC BY 4.0)
 */

//v2023.1
//February 8, 2023
//adapted to work with edge list instead of adjacency matrix
//speed up by removing many exports to disk : only the results are written to file
//includes some optimization tweaks

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define SIZE_S 256
#define MAX_PATH 1000000000.0

void checkFileOpening(FILE *f, char *filename)
{
	if(f == NULL)
	{
		printf("ERROR opening file %s...\n", filename);
		exit(1);
	}
}

int readMatrixFromFile(FILE *f, double **matrix, int size)
{
	int count = 0;
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			count += fscanf(f, "%lf", &matrix[i][j]);
		}
	}
	return count;
}

void buildAdjacencyMatrixFromEdgeListFile(FILE *f, double **matrix)
{
	int count = 0;
	int _src, _tgt;
	double _w;

	while(1)
	{
		if(feof(f))
		{
			break;
		}

		count = fscanf(f, "%d%d%lf", &_src, &_tgt, &_w);
		
		if(count%3 == 0)
		{
			//in the edgelist they are numbered from 1
			_src--;
			_tgt--;

			matrix[_src][_tgt] = _w;
		}
	}
}

void printVector(double *vector, int size)
{
	for(int i = 0; i < size; i++)
	{
		printf("%d: %lf\n", i, vector[i]);
	}
}

void printVector(int *vector, int size)
{
	for(int i = 0; i < size; i++)
	{
		printf("%d: %d\n", i, vector[i]);
	}
}

void printMatrix(double **matrix, int size)
{
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			printf("%.12lf ", matrix[i][j]);	
		}
		printf("\n");
	}
}

void printMatrix(double **matrix, int size, FILE *f)
{
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			fprintf(f, "%.6lf ", matrix[i][j]);	
		}
		fprintf(f, "\n");
	}
}

void getFirstOrderneighborhood(double **matrix, int size, int node, int *neighborhood)
{
	for(int i = 0; i < size; i++)
	{
		if(node != i)
		{
			if(matrix[node][i] != 0 || matrix[i][node] != 0)
			{
				neighborhood[i] = 1;
			}
		}
	}
}

double calculateM(double **matrix, int size, int node, int *neighborhood)
{
	double m = 0;
	
	for(int i = 0; i < size; i++)
	{
		if(neighborhood[i] == 1)
		{
			if(matrix[node][i] != 0)
			{
				m++;
			}
			if(matrix[i][node] != 0)
			{
				m++;
			}
		}

		//testing optizimation - remove this and uncomment next cycle if fails
		for(int j = i + 1; j < size; j++)
		{
			if(i != j && neighborhood[i] == 1 && neighborhood[j] == 1 && node != i && node != j)
			{
				if(matrix[i][j] != 0)
				{
					m++;
				}
				if(matrix[j][i] != 0)
				{
					m++;
				}
			}
		}
	}

	// for(int i = 0; i < size; i++)
	// {
	// 	for(int j = i + 1; j < size; j++)
	// 	{
	// 		if(i != j && neighborhood[i] == 1 && neighborhood[j] == 1 && node != i && node != j)
	// 		{
	// 			if(matrix[i][j] != 0 || matrix[j][i] != 0)
	// 			{
	// 				m += 1;
	// 			}
	// 		}
	// 	}
	// }

	return m;
}

double calculateK(double **matrix, int size, int node, int *neighborhood)
{
	double k = 0;
	int *neighbors = (int *)calloc(size, sizeof(int));

	for(int i = 0; i < size; i++)
	{
		if(neighborhood[i] == 1)
		{
			getFirstOrderneighborhood(matrix, size, i, neighbors);

			for(int j = 0; j < size; j++)
			{
				if(j != node && i != j && neighbors[j] == 1 && neighborhood[j] == 0)
				{
					if(matrix[i][j] != 0)
					{
						k++;
					}
					if(matrix[j][i] != 0)
					{
						k++;
					}
				}
			}
		}
	}

	free(neighbors);

	return k;
}

double nodeStrength(double **matrix, int size, int node)
{
	double degree = 0;

	for(int i = 0; i < size; i++)
	{
		if(matrix[node][i] > 0)
		{
			degree += matrix[node][i];
		}
		if(matrix[i][node] > 0)
		{
			degree += matrix[i][node];
		}
	}

	return degree;
}

double calculateNodeDensity(double **matrix, int size, int node)
{
	//using formula rho = m/(m + k)
	//m - inner links of first order neighborhood of node i
	//k - outer links of first order neighborhood of node i

	int *neighborhood = (int *)calloc(size, sizeof(int));

	getFirstOrderneighborhood(matrix, size, node, neighborhood);

	double m = calculateM(matrix, size, node, neighborhood);
	double k = calculateK(matrix, size, node, neighborhood);
	double strength = nodeStrength(matrix, size, node);

	free(neighborhood);

	return strength * m/(m + k);
}

int numberOfTriangles(double **weight_matrix, int size, int node_i, int node_j)
{
	int triangles = 0;

	for(int i = 0; i < size; i++)
	{
		if((weight_matrix[node_i][i] > 0 || weight_matrix[i][node_i] > 0) && (weight_matrix[node_j][i] > 0 || weight_matrix[i][node_j] > 0) && (weight_matrix[node_i][node_j] > 0 || weight_matrix[node_j][node_i] > 0))
		{
			triangles += 1;
		}
	}

	return triangles;
}

int nodeFullDegree(double **weight_matrix, int size, int node)
{
	int degree = 0;

	for(int i = 0; i < size; i++)
	{
		if(weight_matrix[node][i] > 0 || weight_matrix[i][node] > 0) //TODO: comment properly
		{
			degree += 1;
		}
	}

	return degree;
}

double calculateEdgeClusteringCoefficient(double **weight_matrix, int size, int node_i, int node_j)
{
	int triangles = numberOfTriangles(weight_matrix, size, node_i, node_j);
	int degree_i = nodeFullDegree(weight_matrix, size, node_i);
	int degree_j = nodeFullDegree(weight_matrix, size, node_j);
	int min_degree = ((degree_i - 1) < (degree_j - 1)) ? (degree_i - 1) : (degree_j - 1);

	if(min_degree == 0)
	{
		return (triangles + 1) * 1.0/(min_degree + 0.000001);
	}

	return (triangles + 1) * 1.0/min_degree;
}

void transformWeights2Distances(double **weight_matrix, int size, double **distance_matrix)
{
	double ecc;

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			if(i != j)
			{
				if(weight_matrix[i][j] > 0)
				{
					ecc = calculateEdgeClusteringCoefficient(weight_matrix, size, i, j);
					distance_matrix[i][j] = 1.0/(weight_matrix[i][j] * ecc);
				}
				else
				{
					distance_matrix[i][j] = MAX_PATH;
				}
			}
			else
			{
				distance_matrix[i][j] = MAX_PATH;
			}
		}
	}
}

double calculateShortestPath(double **distance_matrix, int size, double **shortest_path)
{
	double **path = (double **)calloc(size, sizeof(double *));
	for(int i = 0; i < size; i++)
	{
		path[i] = (double *)calloc(size, sizeof(double));
	}

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			if(distance_matrix[i][j] != 0)
			{
				shortest_path[i][j] = distance_matrix[i][j];
			}
			else
			{
				shortest_path[i][j] = MAX_PATH;
			}
		}
	}

	for(int k = 0; k < size; k++)
	{
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < size; j++)
			{
				if(shortest_path[i][j] < shortest_path[i][k] + shortest_path[k][j])
				{
					path[i][j] = shortest_path[i][j];
				}
				else
				{
					path[i][j] = shortest_path[i][k] + shortest_path[k][j];
				}
			}
		}

		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < size; j++)
			{
				if(i != j)
				{
					shortest_path[i][j] = path[i][j];
				}
				else
				{
					shortest_path[i][j] = 0;
				}
			}
		}
	}

	for(int i = 0; i < size; i++)
	{
		free(path[i]);
	}
	free(path);

	double average = 0;
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			average += shortest_path[i][j];
		}
	}

	return average/(size * (size - 1));
}

double findMinimum(double **matrix, int size)
{
	double min = MAX_PATH;
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			if(i != j && matrix[i][j] < min)
			{
				min = matrix[i][j];
			}
		}
	}

	return min;
}

double findMaximum(double **matrix, int size)
{
	double max = 0;
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			if(i != j && matrix[i][j] > max && matrix[i][j] != MAX_PATH)
			{
				max = matrix[i][j];
			}
		}
	}

	return max;
}

void reorder(double *rho, int size)
{
	double tmp;

	for(int i = 0; i < size - 1; i++)
	{
		for(int j = i + 1; j < size; j++)
		{
			if(rho[i] < rho[j])
			{
				tmp = rho[i];
				rho[i] = rho[j];
				rho[j] = tmp;
			}
		}
	}
}

int findIndex(double *rho, int size, double value)
{
	for(int i = 0; i < size; i++)
	{
		if(rho[i] == value)
		{
			return i;
		}
	}

	return -1;
}

void findVoronoiGeneratorPoints(double *rho, double radius, double **shortest_path, int size, int *rho_indices, int *generator_points)
{
	bool largest;
	int cluster = 1;

	for(int node_i = 0; node_i < size; node_i++)
	{
		largest = true;
		for(int node_j = 0; node_j < size; node_j++)
		{
			if(rho_indices[node_i] != rho_indices[node_j])
			{
				if(shortest_path[rho_indices[node_i]][rho_indices[node_j]] <= radius) //taking into account only the links originating in node_i
				{
					if(rho[rho_indices[node_i]] < rho[rho_indices[node_j]] || generator_points[rho_indices[node_j]] != 0)
					{
						largest = false;
						break;
					}
				}
			}
		}

		if(largest)
		{
			generator_points[rho_indices[node_i]] = cluster;
			cluster++;
		}
	}
}

int obtainVoronoiCells(double **shortest_path, int size, double max_path, int *generator_points, int *clustering)
{
	double min;
	int min_index;
	int counter = 0;

	for(int node = 0; node < size; node++)
	{
		if(generator_points[node] == 0)
		{
			min = max_path;
			min_index = 0;

			for(int gp = 0; gp < size; gp++)
			{
				if(generator_points[gp] != 0 && shortest_path[node][gp] <= min) //previous version was < only, so if the shortest path was equal to the radius the node did not get assign into a cluster
				{
					min = shortest_path[node][gp];
					min_index = gp;
				}
			}

			clustering[node] = generator_points[min_index];
		}
		else
		{
			clustering[node] = -1.0 * generator_points[node];
			counter++;
		}
	}

	return counter;
}

double sumOfWeights(double **matrix, int size)
{
	double m = 0;

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			if(matrix[i][j] > 0)
			{
				m += matrix[i][j];
			}
		}
	}

	return m;
}

double nodeDegree(double **matrix, int size, int node)
{
	double degree = 0;
	double tmp = 0;

	for(int i = 0; i < size; i++)
	{
		tmp = 0;
		if(matrix[node][i] > 0)
		{
			tmp += matrix[node][i];
		}
		if(matrix[i][node] > 0)
		{
			tmp += matrix[i][node];
		}

		degree += tmp/2;
	}

	return degree;
}

void outDegree(double **matrix, int size, double *degree)
{
	double tmp;

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			tmp = matrix[i][j];
			if(tmp > 0)
			{
				degree[i] += tmp;
			}
		}
	}
}

void inDegree(double **matrix, int size, double *degree)
{
	double tmp;

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			tmp = matrix[j][i];
			if(tmp > 0)
			{
				degree[i] += tmp;
			}
		}
	}
}

double modularity(double **matrix, int size, int *clustering)
{
	//calculating weighted Newman modularity
	double q = 0;
	double k_i, k_j;

	double m = sumOfWeights(matrix, size);

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			if(i != j)
			{
				k_i = nodeDegree(matrix, size, i);
				k_j = nodeDegree(matrix, size, j);
				q += (matrix[i][j] - (k_i * k_j)/(2 * m)) * (abs(clustering[i]) == abs(clustering[j]) ? 1 : 0); //generator points appear with minus sign, thus abs() is needed
			}
		}
	}
	return q/(2 * m);
}

double directed_modularity(double **matrix, int size, int *clustering)
{
	double q = 0;
	double *k_out = (double *)calloc(size, sizeof(double));
	double *k_in = (double *)calloc(size, sizeof(double));

	double m = sumOfWeights(matrix, size);

	outDegree(matrix, size, k_out);
	inDegree(matrix, size, k_in);

	double weight;

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			if(i != j)
			{
				weight = matrix[i][j];
			}
			else
			{
				weight = 0;
			}
			q += (weight/m - (k_out[i] * k_in[j])/(m * m)) * (abs(clustering[i]) == abs(clustering[j]) ? 1 : 0); //generator points appear with minus sign, thus abs() is needed
		}
	}

	free(k_out);
	free(k_in);

	return q;
}

int getMaxClusterNumber(int *clustering, int size)
{
	int max = abs(clustering[0]);

	for(int i = 1; i < size; i++)
	{
		if(abs(clustering[i]) > max)
		{
			max = abs(clustering[i]);
		}
	}

	return max;
}

void buildClusterMap(int *clustering, int size, int **clustering_map)
{
	int cluster;
	for(int i = 0; i < size; i++)
	{
		cluster = abs(clustering[i]) - 1;
		clustering_map[cluster][i] = 1;
	}
}

double mutual_information(int *clustering, int *clustering_base, int size)
{
	//notation based on Standardized Mutual Information for Clustering Comparisons: One Step Further in Adjustment for Chance by Simone Romano (romano14.pdf)
	double mi = 0;

	int r = getMaxClusterNumber(clustering, size);
	int c = getMaxClusterNumber(clustering_base, size);

	int **M = (int **)calloc(r, sizeof(int *));
	int **clustering_map = (int **)calloc(r, sizeof(int *));
	for(int i = 0; i < r; i++)
	{
		M[i] = (int *)calloc(c, sizeof(int));
		clustering_map[i] = (int *)calloc(size, sizeof(int));
	}
	int **clustering_base_map = (int **)calloc(c, sizeof(int *));
	for(int i = 0; i < c; i++)
	{
		clustering_base_map[i] = (int *)calloc(size, sizeof(int));
	}

	buildClusterMap(clustering, size, clustering_map);
	buildClusterMap(clustering_base, size, clustering_base_map);

	for(int i = 0; i < r; i++)
	{
		for(int j = 0; j < c; j++)
		{
			for(int k = 0; k < size; k++)
			{
				M[i][j] += clustering_map[i][k] * clustering_base_map[j][k];
			}
		}
	}

	int *A = (int *)calloc(r, sizeof(int));
	int *B = (int *)calloc(c, sizeof(int));

	double Ha = 0, Hb = 0;

	for(int i = 0; i < r; i++)
	{
		for(int j = 0; j < c; j++)
		{
			A[i] += M[i][j];
		}
		if(A[i] != 0)
		{			
			Ha += A[i] * 1.0/size * log2(A[i] * 1.0/size);
		}
	}

	for(int j = 0; j < c; j++)
	{
		for(int i = 0; i < r; i++)
		{
			B[j] += M[i][j];
		}
		if(B[j] != 0)
		{			
			Hb += B[j] * 1.0/size * log2(B[j] * 1.0/size);
		}
	}

	for(int i = 0; i < r; i++)
	{
		for(int j = 0; j < c; j++)
		{
			if(M[i][j] != 0)
			{
				mi += M[i][j] * 1.0/size * log2(M[i][j] * size * 1.0/(A[i] * B[j]));
			}
		}
	}

	Ha *= -1.0;
	Hb *= -1.0;

	for(int i = 0; i < c; i++)
	{
		free(clustering_base_map[i]);
	}
	for(int i = 0; i < r; i++)
	{
		free(clustering_map[i]);
		free(M[i]);
	}
	free(clustering_base_map);
	free(clustering_map);
	free(A);
	free(B);
	free(M);

	double max = (Ha <= Hb) ? Hb : Ha;
	if(max == 0)
	{
		return mi;
	}
	
	return mi/max;
}

#ifdef __cplusplus
extern "C" {
#endif

__attribute__((visibility("default")))
int* community_voronoi(double** matrix, int size) {
    if (!matrix || size <= 0) {
        return NULL;
    }

    // Initialize all pointers to NULL
    double** distanceMatrix = NULL;
    double** shortestPath = NULL;
    double* rho = NULL;
    double* reorderedRho = NULL;
    int* rho_indices = NULL;
    int* best_clustering = NULL;
    int* result = NULL;

    try {
        // Allocate result array first
        result = (int*)calloc(size, sizeof(int));
        if (!result) {
        //     goto cleanup;
		
			// Clean up all allocated memory
			if (distanceMatrix) {
				for (int i = 0; i < size; i++) {
					if (distanceMatrix[i]) free(distanceMatrix[i]);
				}
				free(distanceMatrix);
			}
			if (shortestPath) {
				for (int i = 0; i < size; i++) {
					if (shortestPath[i]) free(shortestPath[i]);
				}
				free(shortestPath);
			}
			free(rho);
			free(reorderedRho);
			free(rho_indices);
			free(best_clustering);
        }

        // Allocate matrices and vectors
        distanceMatrix = (double**)calloc(size, sizeof(double*));
        shortestPath = (double**)calloc(size, sizeof(double*));
        rho = (double*)calloc(size, sizeof(double));
        reorderedRho = (double*)calloc(size, sizeof(double));
        rho_indices = (int*)calloc(size, sizeof(int));
        
        // if (!distanceMatrix || !shortestPath || !rho || !reorderedRho || !rho_indices) {
        //     goto cleanup;
        // }

        // Allocate matrix rows
        for (int i = 0; i < size; i++) {
            if (distanceMatrix) distanceMatrix[i] = (double*)calloc(size, sizeof(double));
            if (shortestPath) shortestPath[i] = (double*)calloc(size, sizeof(double));
            if (!distanceMatrix[i] || !shortestPath[i]) {
            //     goto cleanup;

				// Clean up all allocated memory
				if (distanceMatrix) {
					for (int i = 0; i < size; i++) {
						if (distanceMatrix[i]) free(distanceMatrix[i]);
					}
					free(distanceMatrix);
				}
				if (shortestPath) {
					for (int i = 0; i < size; i++) {
						if (shortestPath[i]) free(shortestPath[i]);
					}
					free(shortestPath);
				}
				free(rho);
				free(reorderedRho);
				free(rho_indices);
				free(best_clustering);
            }
        }

        // Calculate node densities
        for (int i = 0; i < size; i++) {
            rho[i] = calculateNodeDensity(matrix, size, i);
            reorderedRho[i] = rho[i];
        }

        // Order densities and get indices
        reorder(reorderedRho, size);
        for (int i = 0; i < size; i++) {
            rho_indices[i] = findIndex(rho, size, reorderedRho[i]);
        }

        // Transform weights to distances
        transformWeights2Distances(matrix, size, distanceMatrix);
        double avg_path = calculateShortestPath(distanceMatrix, size, shortestPath);

        // Calculate paths
        double min_path = findMinimum(shortestPath, size);
        double max_path = (avg_path > 200) ? 200 : avg_path;
        double min_R = min_path;
        double max_R = (avg_path > 200) ? 200 : avg_path;

        // Initialize best values
        double best_modularity = -1;
        best_clustering = (int*)calloc(size, sizeof(int));
        if (!best_clustering) {
        //     goto cleanup;

			// Clean up all allocated memory
			if (distanceMatrix) {
				for (int i = 0; i < size; i++) {
					if (distanceMatrix[i]) free(distanceMatrix[i]);
				}
				free(distanceMatrix);
			}
			if (shortestPath) {
				for (int i = 0; i < size; i++) {
					if (shortestPath[i]) free(shortestPath[i]);
				}
				free(shortestPath);
			}
			free(rho);
			free(reorderedRho);
			free(rho_indices);
			free(best_clustering);
        }

        // Iterate through different radii
        for (double r = min_R; r <= max_R; r += 0.25) {
            int* generator_points = (int*)calloc(size, sizeof(int));
            int* clustering = (int*)calloc(size, sizeof(int));
            if (!generator_points || !clustering) {
                free(generator_points);
                free(clustering);
                continue;
            }

            findVoronoiGeneratorPoints(rho, r, shortestPath, size, rho_indices, generator_points);
            int clusterNumber = obtainVoronoiCells(shortestPath, size, max_path, generator_points, clustering);

            double current_modularity = directed_modularity(matrix, size, clustering);

            if (current_modularity > best_modularity) {
                best_modularity = current_modularity;
                memcpy(best_clustering, clustering, size * sizeof(int));
            }

            free(generator_points);
            free(clustering);

            if (clusterNumber == 1) break;
        }

        // Copy best clustering to result
        memcpy(result, best_clustering, size * sizeof(int));

// cleanup:
//         // Clean up all allocated memory
//         if (distanceMatrix) {
//             for (int i = 0; i < size; i++) {
//                 if (distanceMatrix[i]) free(distanceMatrix[i]);
//             }
//             free(distanceMatrix);
//         }
//         if (shortestPath) {
//             for (int i = 0; i < size; i++) {
//                 if (shortestPath[i]) free(shortestPath[i]);
//             }
//             free(shortestPath);
//         }
//         free(rho);
//         free(reorderedRho);
//         free(rho_indices);
//         free(best_clustering);

        return result;

    } catch (...) {
        if (result) free(result);
        return NULL;
    }
}

#ifdef __cplusplus
}
#endif




int main(int argc, const char *argv[])
{
	if(argc != 10)
	{
		printf("> usage %s [number of nodes] [average degree] [mixing parameter] [type] [ins] [outs] [basename] [start index] [end index]\n", argv[0]);
		printf(">> [number of nodes]: number of nodes of the network\n");
		printf(">> [average degree]: average degree of a node in the network\n");
		printf(">> [mixing parameter]: mixing parameter of the network\n");
		printf(">> [type]: distribution type of weights: normal (for normal distribution) or power (for power-law distribution)\n");
		printf(">> [ins]: mu of the inner-cluster link weight-distribution (w/o decimal dot)\n");
		printf(">> [outs]: mu of the inter-cluster link weight-distribution (w/o decimal dot)\n");
		printf(">> [basename]: basename of adjacency matrix files (extension is condidered to be txt\n");
		printf(">> [start index]: start of iteration\n");
		printf(">> [end index]: end of iteration\n");
		exit(0);
	}

	time_t t = time(NULL);
  	struct tm tm;
  	struct stat sb = {0};

  	char *foldername = (char *)malloc(SIZE_S * sizeof(char));
	char *filename = (char *)malloc(SIZE_S * sizeof(char));
	char *base = (char *)malloc(SIZE_S * sizeof(char));
	char *rtype = (char *)malloc(SIZE_S * sizeof(char));
	char *mu = (char *)malloc(SIZE_S * sizeof(char));
	char *ins = (char *)malloc(SIZE_S * sizeof(char));
	char *outs = (char *)malloc(SIZE_S * sizeof(char));

	int N = atoi(argv[1]);
	int k = atoi(argv[2]);
	strcpy(mu, argv[3]);
	strcpy(rtype, argv[4]);
	strcpy(ins, argv[5]);
	strcpy(outs, argv[6]);
	strcpy(base, argv[7]);
	int start = atoi(argv[8]);
	int end = atoi(argv[9]);

	sprintf(foldername, "N%d_k%d_mu%s_%s%s%s", N, k, mu, rtype, ins, outs);

	char *resultsFolder = (char *)malloc(SIZE_S * sizeof(char));
	sprintf(resultsFolder, "results/%s", foldername);
	if(stat(resultsFolder, &sb) == -1) // && S_ISDIR(sb.st_mode))
	{
		int dir_err = mkdir(resultsFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if(-1 == dir_err)
		{
			printf("> [E] : error creating folder : %s\n", resultsFolder);
		}
	}
	free(resultsFolder);

	for(int index = start; index <= end; index++)
	{
		double **adjacencyMatrix = (double **)calloc(N, sizeof(double *));
		double **distanceMatrix = (double **)calloc(N, sizeof(double *));
		double **shortestPath = (double **)calloc(N, sizeof(double *));
		for(int i = 0; i < N; i++)
		{
			adjacencyMatrix[i] = (double *)calloc(N, sizeof(double));
			distanceMatrix[i] = (double *)calloc(N, sizeof(double));
			shortestPath[i] = (double *)calloc(N, sizeof(double));
		}
		char *filenameWPath = (char *)malloc(SIZE_S * sizeof(char));

		sprintf(filename, "%s_%d", base, index);
		sprintf(filenameWPath, "benchmarks/%s/%s.txt", foldername, filename);
		
		FILE *f = fopen(filenameWPath, "r");
		checkFileOpening(f, filenameWPath);
		
		buildAdjacencyMatrixFromEdgeListFile(f, adjacencyMatrix);

		fclose(f);

		int *clustering_base = (int *)calloc(N, sizeof(int));
		char *baseClusteringFile = (char *)malloc(SIZE_S * sizeof(char));

		sprintf(baseClusteringFile, "generator/N%d_k%d_base_mu%s/community_%d.txt", N, k, mu, index);
		f = fopen(baseClusteringFile, "r");
		checkFileOpening(f, baseClusteringFile);
		
		int idummy, icount;
		for(int i = 0; i < N; i++)
		{
			icount = fscanf(f, "%d%d", &idummy, &clustering_base[i]);
		}
		fclose(f);
		free(baseClusteringFile);

		double *rho = (double *)calloc(N, sizeof(double));
		double *reorderedRho = (double *)calloc(N, sizeof(double));
		int *rho_indices = (int *)calloc(N, sizeof(int));
		for(int i = 0; i < N; i++)
		{
			rho[i] = calculateNodeDensity(adjacencyMatrix, N, i);
			reorderedRho[i] = rho[i];
		}
		reorder(reorderedRho, N);
		for(int i = 0; i < N; i++)
		{
			rho_indices[i] = findIndex(rho, N, reorderedRho[i]);
		}

		free(reorderedRho);

		transformWeights2Distances(adjacencyMatrix, N, distanceMatrix);

		double avg_path = calculateShortestPath(distanceMatrix, N, shortestPath);

		for(int i = 0; i < N; i++)
		{
			free(distanceMatrix[i]);
		}
		free(distanceMatrix);

		double min_path = findMinimum(shortestPath, N);
		double max_path = findMaximum(shortestPath, N);

		double min_R = min_path;
		double max_R = (avg_path > 200) ? 200 : avg_path;

		double best_radius;
		int best_cluster_nr;
		double best_modularity = -1;
		double best_mutual_information = -1;
		double *best_rho = (double *)calloc(N, sizeof(double));
		int *best_clustering = (int *)calloc(N, sizeof(int));

		int s = 0;
		double step = 0.25;

		int totalSteps = 1; //(int)((max_R - min_R)/step) + 1;
		double *Q = (double *)calloc(totalSteps, sizeof(double));
		double *MI = (double *)calloc(totalSteps, sizeof(double));
		
		for(double r = min_R; r <= max_R; r += step)
		{
			int *generator_points = (int *)calloc(N, sizeof(int));
			int *clustering = (int *)calloc(N, sizeof(int));

			findVoronoiGeneratorPoints(rho, r, shortestPath, N, rho_indices, generator_points);

			int clusterNumber = obtainVoronoiCells(shortestPath, N, max_path, generator_points, clustering);
			
			Q[s] = directed_modularity(adjacencyMatrix, N, clustering);
			MI[s] = mutual_information(clustering, clustering_base, N);

			totalSteps++;
			Q = (double *)realloc(Q, totalSteps * sizeof(double));
			MI = (double *)realloc(MI, totalSteps * sizeof(double));

			if(Q[s] > best_modularity)
			{
				best_modularity = Q[s];
				best_radius = r;
				best_cluster_nr = clusterNumber;
				best_mutual_information = MI[s];

				for(int i = 0; i < N; i++)
				{
					best_rho[i] = rho[i];
					best_clustering[i] = clustering[i];
				}
			}
			else
			{
				if(clusterNumber == 1)
				{
					break;
				}
			}

			free(generator_points);
			free(clustering);

			s++;
		}

		sprintf(filenameWPath, "results/N%d_k%d_mu%s_%s%s%s/%s.txt", N, k, mu, rtype, ins, outs, filename);
		f = fopen(filenameWPath, "w");
		//steps, radius, cluster count, modularity, mutual information
		fprintf(f, "%d\t%.16lf\t%d\t%lf\t%lf\n", s - 1, best_radius, best_cluster_nr, best_modularity, best_mutual_information);
		fclose(f);

		sprintf(filenameWPath, "results/N%d_k%d_mu%s_%s%s%s/%s_clustering.txt", N, k, mu, rtype, ins, outs, filename);
		f = fopen(filenameWPath, "w");
		for(int i = 0; i < N; i++)
		{
			fprintf(f, "%d\n", best_clustering[i]);
		}
		fclose(f);

		sprintf(filenameWPath, "results/N%d_k%d_mu%s_%s%s%s/%s_mmi.txt", N, k, mu, rtype, ins, outs, filename);
		f = fopen(filenameWPath, "w");
		for(int i = 0; i < totalSteps; i++)
		{
			fprintf(f, "%lf\t%lf\t%lf\n", min_R + (step * i), Q[i], MI[i]);
		}
		fclose(f);

		for(int i = 0; i < N; i++)
		{			
			free(adjacencyMatrix[i]);
			free(shortestPath[i]);
		}
		free(adjacencyMatrix);
		free(shortestPath);
		free(rho);	
		free(rho_indices);
		free(best_rho);
		free(best_clustering);
		free(filenameWPath);
		free(clustering_base);
		free(Q);
		free(MI);
	}

	free(filename);
	free(base);
	free(rtype);
	free(ins);
	free(outs);
	free(mu);

	return 0;
}