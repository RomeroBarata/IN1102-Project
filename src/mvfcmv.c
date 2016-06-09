#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "util.h"
#include "stex.h"

#define BUFF_SIZE 1024

// Set this to make the program give more output.
bool verbose;
// Maximum number of iteration for a given instance.
size_t max_iter;
// Number of objects.
int objc;
// Number of clusters.
size_t clustc;
// Parameter 'm' of the algorithm.
double mfuz;
// mfuzval = 1.0 / (1.0 - 'mfuz')
double mfuzval;
// Number of dissimilarity matrices.
int dmatrixc;
// 'dmatrixc' dissimilarity matrices.
double ***dmatrix;
// Cardinality of medoids.
size_t medoids_card;
// Medoids for each cluster for each dissimilarity matrix.
size_t ***medoids;
// Relevance weights for each cluster for each dissimilarity matrix.
double **weights;
// Membership matrix for each object.
double **memb;
// Parameter 'epsilon' of the algortihm. Stopping criteria.
double epsilon;
// Parameter 'theta' of the algorithm.
double theta;
// Adequacy broken down by each object.
double *parc_obj_adeq;
// Adequacy broken down by each cluster.
double *parc_cluster_adeq;
// Previous computed adequacy.
double prev_adeq;

// Prints and checks 'weights'.
void print_weights(double **weights) {
	printf("Weights:\n");
	size_t j;
	size_t k;
	double prod;
	for(k = 0; k < clustc; ++k) {
		prod = 1.0;
		for(j = 0; j < dmatrixc; ++j) {
			if(weights[k][j] < 0.0) {
				printf("!");
			}
			printf("%lf ", weights[k][j]);
			prod *= weights[k][j];
		}
		printf("[%lf]", prod);
		if(!deq(prod, 1.0)) {
			printf(" =/= 1.0?\n");
		} else {
			printf("\n");
		}
	}
}

// Randomly initialize 'medoids'.
void init_medoids() {
	size_t i;
    size_t j;
	size_t e;
	size_t k;
	int obj;
	bool chosen[objc];
	for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            for(i = 0; i < objc; ++i) {
                chosen[i] = false;
            }
            for(e = 0; e < medoids_card; ++e) {
                do {
                    obj = rand() % objc;
                } while(chosen[obj]);
                medoids[k][j][e] = obj;
                chosen[obj] = true;
            }
        }
	}
}

// Prints and checks 'medoids'.
void print_medoids(size_t ***medoids) {
    printf("Medoids:\n");
    size_t e;
    size_t j;
    size_t k;
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            printf("{ ");
            for(e = 0; e < medoids_card; ++e) {
                printf("%d ", medoids[k][j][e]);
            }
            printf("} ");
        }
        printf("\n");
    }
}

// Prints and checks 'memb'.
void print_memb(double **memb) {
	printf("Membership:\n");
	size_t i;
	size_t k;
	double sum;
	for(i = 0; i < objc; ++i) {
		printf("%u: ", i);
		sum = 0.0;
		for(k = 0; k < clustc; ++k) {
			printf("%lf ", memb[i][k]);
			sum += memb[i][k];
		}
		printf("[%lf]", sum);
		if(!deq(sum, 1.0)) {
			printf("*\n");
		} else {
			printf("\n");
		}
	}
}

// Updates 'memb'.
void update_memb() {
	size_t e;
	size_t h;
	size_t i;
	size_t j;
	size_t k;
	size_t a_card;
	double sums[clustc];
	double sumd;
	double val;
	for(i = 0; i < objc; ++i) {
		a_card = 0;
		for(k = 0; k < clustc; ++k) {
			sums[k] = 0.0;
			for(j = 0; j < dmatrixc; ++j) {
				sumd = 0.0;
				for(e = 0; e < medoids_card; ++e) {
					sumd += dmatrix[j][i][medoids[k][j][e]];
				}
				sums[k] += weights[k][j] * sumd;
			}
			if(deq(sums[k], 0.0)) {
				++a_card;
			}
		}
		if(a_card) {
			printf("Object %u has zero val.\n", i);
			val = 1.0 / ((double) a_card);
			for(k = 0; k < clustc; ++k) {
				if(deq(sums[k], 0.0)) {
					memb[i][k] = val;
				} else {
					memb[i][k] = 0.0;
				}
			}
		} else {
			for(k = 0; k < clustc; ++k) {
				memb[i][k] = 0.0;
				for(h = 0; h < clustc; ++h) {
					memb[i][k] += pow(sums[k] / sums[h], mfuzval);
				}
				memb[i][k] = 1.0 / memb[i][k];
			}
		}
	}
}

// Computes the adequacy and stores the contribution of each cluster
// into 'parc_cluster_adeq'.
double adequacy_cluster(bool check) {
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    double sum;
    double sumd;
	double sumw;
    double adeq = 0.0;
	for(k = 0; k < clustc; ++k) {
		sum = 0.0;
		for(i = 0; i < objc; ++i) {
			sumw = 0.0;
            for(j = 0; j < dmatrixc; ++j) {
                sumd = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    sumd += dmatrix[j][i][medoids[k][j][e]];
                }
                sumw += weights[k][j] * sumd;
            }
            sum += pow(memb[i][k], mfuz) * sumw;
        }
		if(check) {
			if(dlt(parc_cluster_adeq[k], sum)) {
				printf("Msg: adequacy for cluster %d is greater "
                        "than previous (%.15lf).\n", k,
						sum - parc_cluster_adeq[k]);
			}
		}
		parc_cluster_adeq[k] = sum;
		adeq += sum;
    }
    if(check) {
        if(dlt(prev_adeq, adeq)) {
            printf("Msg: overall adequacy is greater than "
                    "previous (%.15lf).\n", adeq - prev_adeq);
        }
    }
    prev_adeq = adeq;
    return adeq;
}

// Computes the adequacy and stores the contribution of each object 
// into 'parc_obj_adeq'.
double adequacy_obj(bool check) {
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    double sum;
    double sumd;
	double sumw;
    double adeq = 0.0;
	for(i = 0; i < objc; ++i) {
		sum = 0.0;
		for(k = 0; k < clustc; ++k) {
			sumw = 0.0;
            for(j = 0; j < dmatrixc; ++j) {
                sumd = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    sumd += dmatrix[j][i][medoids[k][j][e]];
                }
                sumw += weights[k][j] * sumd;
            }
            sum += pow(memb[i][k], mfuz) * sumw;
        }
		if(check) {
			if(dlt(parc_obj_adeq[i], sum)) {
				printf("Msg: adequacy for object %d is greater "
                        "than previous (%.15lf).\n", i,
						sum - parc_obj_adeq[i]);
			}
		}
		parc_obj_adeq[i] = sum;
		adeq += sum;
    }
    if(check) {
        if(dlt(prev_adeq, adeq)) {
            printf("Msg: overall adequacy is greater than "
                    "previous (%.15lf).\n", adeq - prev_adeq);
        }
    }
    prev_adeq = adeq;
    return adeq;
}

// Structure used in 'update_medoids' to store the object and its
// computed value.
typedef struct objnval {
	size_t obj;
	double val;
} objnval;

// Comparator for 'objnval', used in the sort function by
// 'update_medoids'.
static int objnval_cmp(const void *p1, const void *p2) {
	const objnval *a = (const objnval *) p1;
	const objnval *b = (const objnval *) p2;
	return (a->val > b->val) - (a->val < b->val);
}

// Updates 'medoids' by selecting the best candidates.
void update_medoids() {
	size_t h;
	size_t i;
	size_t j;
	size_t k;
	objnval candidates[objc];
	for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            for(h = 0; h < objc; ++h) {
                candidates[h].obj = h;
                candidates[h].val = 0.0;
                for(i = 0; i < objc; ++i) {
                    candidates[h].val += pow(memb[i][k], mfuz) * dmatrix[j][i][h];
                }
            }
            qsort(candidates, objc, sizeof(objnval), objnval_cmp);
            for(h = 0; h < medoids_card; ++h) {
                medoids[k][j][h] = candidates[h].obj;
            }
        }
	}
}

// Updates 'weights'.
void update_weights() {
	size_t e;
	size_t i;
	size_t j;
	size_t k;
	size_t v_card;
	double sums[dmatrixc];
	double sumd;
	double chi;
	double prod_num;
	double exp;
	for(k = 0; k < clustc; ++k) {
		v_card = 0;
		for(j = 0; j < dmatrixc; ++j) {
			sums[j] = 0.0;
			for(i = 0; i < objc; ++i) {
				sumd = 0.0;
				for(e = 0; e < medoids_card; ++e) {
					sumd += dmatrix[j][i][medoids[k][j][e]];
				}
				sums[j] += pow(memb[i][k], mfuz) * sumd;
			}
			if(sums[j] <= theta) {
				++v_card;
			}
		}
		exp = 1.0 / ((double)(dmatrixc - v_card));
		chi = 1.0;
		prod_num = 1.0;
		for(j = 0; j < dmatrixc; ++j) {
			if(sums[j] <= theta) {
				chi *= pow(weights[k][j], exp);
			} else {
				prod_num *= pow(sums[j], exp);
			}
		}
		prod_num = (1.0 / chi) * prod_num;
		for(j = 0; j < dmatrixc; ++j) {
			if(sums[j] > theta) {
				weights[k][j] = prod_num / sums[j];
			}
		}
	}
}

// Main loop for an instance of the algorithm.
double run() {
	size_t i;
	size_t j;
	size_t k;
	printf("Initialization.\n");
	init_medoids();
    if(verbose) print_medoids(medoids);
	for(k = 0; k < clustc; ++k) {
		for(j = 0; j < dmatrixc; ++j) {
			weights[k][j] = 1.0;
		}
	}
	if(verbose) print_weights(weights);
	update_memb();
	if(verbose) print_memb(memb);
	double prev_adeq = 0.0;
	double adeq = adequacy_obj(false);
	printf("Adequacy: %.20lf\n", adeq);
    double diff = fabs(adeq - prev_adeq);
	for(i = 1; i <= max_iter && diff > epsilon; ++i) {
        printf("Iteration %d.\n", i);
        prev_adeq = adeq;
		adequacy_cluster(false);
        update_medoids();
		adeq = adequacy_cluster(true);
        if(verbose) {
            print_medoids(medoids);
            printf("Adequacy1: %.20lf\n", adeq);
        }
		adequacy_cluster(false);
        update_weights();
		adeq = adequacy_cluster(true);
        if(verbose) {
            print_weights(weights);
            printf("Adequacy2: %.20lf\n", adeq);
        }
		adequacy_obj(false);
        update_memb();
		adeq = adequacy_obj(true);
        if(verbose) print_memb(memb);
        printf("Adequacy: %.20lf\n", adeq);
        if(dgt(adeq, prev_adeq)) {
            printf("Warn: current adequacy is greater than "
                    "previous iteration (%.20lf)\n",
                    adeq - prev_adeq);
        }
        diff = fabs(adeq - prev_adeq);
	}
    printf("Adequacy difference threshold reached (%.20lf).\n",
            diff);
    return adeq;
}

// Computes and prints the global energy.
// lenergy - J value, calculated using the medoids used by the
// algorithm
// genergy - T value, calculated using the global medoids computed
// in this function
void global_energy() {
    size_t e;
    size_t h;
    size_t i;
    size_t j;
    size_t k;
    // Defining global medoids
    objnval candidates[objc];    
    size_t global_medoids[dmatrixc][medoids_card];
    for(j = 0; j < dmatrixc; ++j) {
        for(h = 0; h < objc; ++h) {
            candidates[h].obj = h;
            candidates[h].val = 0.0;
            for(k = 0; k < clustc; ++k) {
                for(i = 0; i < objc; ++i) {
                    candidates[h].val += pow(memb[i][k], mfuz) *
                        weights[k][j] * dmatrix[j][i][h];
                }
            }
        }
        qsort(candidates, objc, sizeof(objnval), objnval_cmp);
        for(h = 0; h < medoids_card; ++h) {
            global_medoids[j][h] = candidates[h].obj;
        }
    }
    printf("Global medoids:\n");
    for(j = 0; j < dmatrixc; ++j) {
        printf("%d:", j + 1);
        for(e = 0; e < medoids_card; ++e) {
            printf(" %d", global_medoids[j][e]);
        }
        printf("\n");
    }
    printf("\n");
    // Global energy per cluster and matrix
    double genergy[clustc][dmatrixc];
    double lenergy[clustc][dmatrixc];
    double sumd_global;
    double sumd_local;
    double val;
    double mtx_genergy[dmatrixc];
    double mtx_lenergy[dmatrixc];
    for(j = 0; j < dmatrixc; ++j) {
        mtx_genergy[j] = 0.0;
        mtx_lenergy[j] = 0.0;
    }
    double clust_genergy[clustc];
    double clust_lenergy[clustc];
    for(k = 0; k < clustc; ++k) {
        clust_genergy[k] = 0.0;
        clust_lenergy[k] = 0.0;
    }
    double total_genergy = 0.0;
    double total_lenergy = 0.0;
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            genergy[k][j] = 0.0;
            lenergy[k][j] = 0.0;
            for(i = 0; i < objc; ++i) {
                sumd_global = 0.0;
                sumd_local = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    sumd_global +=
                        dmatrix[j][i][global_medoids[j][e]];
                    sumd_local += dmatrix[j][i][medoids[k][j][e]];
                }
                val = pow(memb[i][k], mfuz) * weights[k][j];
                genergy[k][j] += val * sumd_global;
                lenergy[k][j] += val * sumd_local;
            }
            clust_genergy[k] += genergy[k][j];
            clust_lenergy[k] += lenergy[k][j];
            mtx_genergy[j] += genergy[k][j];
            mtx_lenergy[j] += lenergy[k][j];
            total_genergy += genergy[k][j];
            total_lenergy += lenergy[k][j];
        }
    }
    printf("Global energy matrix:\n");
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            printf("%lf ", genergy[k][j]);
        }
        printf("\n");
    }
    printf("Local energy matrix:\n");
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            printf("%lf ", lenergy[k][j]);
        }
        printf("\n");
    }
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            if(dgt(lenergy[k][j], genergy[k][j])) {
                printf("Msg: lenergy > genergy for cluster %d and matrix %d.\n", k, j);
            }
        }
    }
    printf("\n");
    double sum = 0.0;
    printf("Cluster global energy:\n");
    for(k = 0; k < clustc; ++k) {
        sum += clust_genergy[k];
        printf("%lf ", clust_genergy[k]);
    }
    printf("[%lf]\n", sum);
    sum = 0.0;
    printf("Cluster local energy:\n");
    for(k = 0; k < clustc; ++k) {
        sum += clust_lenergy[k];
        printf("%lf ", clust_lenergy[k]);
    }
    printf("[%lf]\n", sum);
    for(k = 0; k < clustc; ++k) {
        if(dgt(clust_lenergy[k], clust_genergy[k])) {
            printf("Warn: clust_lenergy > clust_genergy for cluster %d.\n", k);
        }
    }
    printf("\n");
    sum = 0.0;
    printf("Matrix global energy:\n");
    for(j = 0; j < dmatrixc; ++j) {
        sum += mtx_genergy[j];
        printf("%lf ", mtx_genergy[j]);
    }
    printf("[%lf]\n", sum);
    sum = 0.0;
    printf("Matrix local energy:\n");
    for(j = 0; j < dmatrixc; ++j) {
        sum += mtx_lenergy[j];
        printf("%lf ", mtx_lenergy[j]);
    }
    printf("[%lf]\n", sum);
    for(j = 0; j < dmatrixc; ++j) {
        if(dgt(mtx_lenergy[j], mtx_genergy[j])) {
            printf("Warn: mtx_lenergy > mtx_genergy for matrix %d.\n", j);
        }
    }
    printf("\n");
    printf("Global energy: %lf\n", total_genergy);
    printf("Local energy: %lf\n", total_lenergy);
    if(dgt(total_lenergy, total_genergy)) {
        printf("Warn: total_lenergy > total_genergy.\n");
    }
    printf("\n");
    printf("Matrix global heterogeneity index:\n");
    for(j = 0; j < dmatrixc; ++j) {
        printf("%lf ", (1.0 - (mtx_lenergy[j] / mtx_genergy[j])));
    }
    printf("\n");
    printf("\n");
    printf("Global heterogeneity index: %lf\n",
            (1.0 - (total_lenergy / total_genergy)));
    printf("\n");
    printf("Cluster heterogeinety indexes:\n");
    printf("T:\n");
    for(k = 0; k < clustc; ++k) {
        printf("%lf ", clust_genergy[k] / total_genergy);
    }
    printf("\n");
    printf("J:\n");
    for(k = 0; k < clustc; ++k) {
        printf("%lf ", clust_lenergy[k] / total_lenergy);
    }
    printf("\n");
    printf("B:\n");
    for(k = 0; k < clustc; ++k) {
        printf("%lf ", (clust_genergy[k] - clust_lenergy[k]) /
                    (total_genergy - total_lenergy));
    }
    printf("\n");
    printf("Q:\n");
    for(k = 0; k < clustc; ++k) {
        printf("%lf ", (1.0 - (clust_lenergy[k] / clust_genergy[k])));
    }
    printf("\n");
    printf("\n");
    printf("Cluster heterogeinety indexes for matrices:\n");
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            printf("%lf ", (1.0 - (lenergy[k][j] / genergy[k][j])));
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
	verbose = false;
	int insts;
    // Opening and reading config file.
    FILE *cfgfile = fopen(argv[1], "r");
    if(!cfgfile) {
        printf("Error: could not open config file.\n");
        return 1;
    }
    fscanf(cfgfile, "%d", &objc);
    if(objc <= 0) {
        printf("Error: objc <= 0.\n");
        return 2;
    }
    int labels[objc];
    // reading labels
    int classc;
    fscanf(cfgfile, "%d", &classc);
	size_t i;
    for(i = 0; i < objc; ++i) {
        fscanf(cfgfile, "%d", &labels[i]);
    }
    // reading labels end
    fscanf(cfgfile, "%d", &dmatrixc);
    if(dmatrixc <= 0) {
        printf("Error: dmatrixc <= 0.\n");
        return 2;
    }
    char dmtx_file_name[dmatrixc][BUFF_SIZE];
	size_t j;
    for(j = 0; j < dmatrixc; ++j) {
        fscanf(cfgfile, "%s", dmtx_file_name[j]);
    }
    char out_file_name[BUFF_SIZE];
    fscanf(cfgfile, "%s", out_file_name);
    fscanf(cfgfile, "%d", &clustc);
    if(clustc <= 0) {
        printf("Error: clustc <= 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%d", &medoids_card);
    if(medoids_card <= 0) {
        printf("Error: medoids_card <= 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%d", &insts);
    if(insts <= 0) {
        printf("Error: insts <= 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%lf", &theta);
    if(dlt(theta, 0.0)) {
        printf("Error: theta < 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%d", &max_iter);
    fscanf(cfgfile, "%lf", &epsilon);
    if(dlt(epsilon, 0.0)) {
        printf("Error: epsilon < 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%lf", &mfuz);
    if(!dgt(mfuz, 0.0)) {
        printf("Error: mfuz <= 0.\n");
        return 2;
    }
    fclose(cfgfile);
    // Done reading config file.
    freopen(out_file_name, "w", stdout);
	mfuzval = 1.0 / (mfuz - 1.0);
    printf("######Config summary:######\n");
    printf("Number of clusters: %d.\n", clustc);
    printf("Medoids cardinality: %d.\n", medoids_card);
    printf("Number of iterations: %d.\n", max_iter);
    printf("Epsilon: %.15lf.\n", epsilon);
    printf("Theta: %.15lf.\n", theta);
    printf("Parameter m: %.15lf.\n", mfuz);
    printf("Number of instances: %d.\n", insts);
    printf("###########################\n");
	size_t k;
	// Allocating memory start
	parc_cluster_adeq = malloc(sizeof(double) * clustc);
    parc_obj_adeq = malloc(sizeof(double) * objc);
	dmatrix = malloc(sizeof(double **) * dmatrixc);
	for(j = 0; j < dmatrixc; ++j) {
		dmatrix[j] = malloc(sizeof(double *) * objc);
		for(i = 0; i < objc; ++i) {
			dmatrix[j][i] = malloc(sizeof(double) * objc);
		}
	}
	medoids = malloc(sizeof(size_t **) * clustc);
	size_t ***best_medoids = malloc(sizeof(size_t **) * clustc);
	for(k = 0; k < clustc; ++k) {
		medoids[k] = malloc(sizeof(size_t *) * dmatrixc);
		best_medoids[k] = malloc(sizeof(size_t *) * dmatrixc);
        for(j = 0; j < dmatrixc; ++j) {
            medoids[k][j] = malloc(sizeof(size_t) * medoids_card);
            best_medoids[k][j] = malloc(sizeof(size_t) * medoids_card);
        }
	}
	weights = malloc(sizeof(double *) * clustc);
	double **best_weights = malloc(sizeof(double *) * clustc);
	for(k = 0; k < clustc; ++k) {
		weights[k] = malloc(sizeof(double) * dmatrixc);
		best_weights[k] = malloc(sizeof(double) * dmatrixc);
	}
	memb = malloc(sizeof(double *) * objc);
	double **best_memb = malloc(sizeof(double *) * objc);
	for(i = 0; i < objc; ++i) {
		memb[i] = malloc(sizeof(double) * clustc);
		best_memb[i] = malloc(sizeof(double) * clustc);
	}
	// Allocating memory end
	for(j = 0; j < dmatrixc; ++j) {
		if(!load_data(dmtx_file_name[j], dmatrix[j], objc, objc)) {
			printf("Error: could not load %s.\n", dmtx_file_name[j]);
			goto END;
		}
	}
    size_t best_inst;
    double best_inst_adeq;
    double cur_inst_adeq;
	srand(time(NULL)); // Random seed.
    // Start main program loop.
	for(i = 1; i <= insts; ++i) {
		printf("Instance %u:\n", i);
		cur_inst_adeq = run();
        if(i == 1 || cur_inst_adeq < best_inst_adeq) {
            // Saves the best configuration based on the adequacy.
            mtxcpy_d(best_memb, memb, objc, clustc);
            mtxcpy_d(best_weights, weights, clustc, dmatrixc);
            for(k = 0; k < clustc; ++k) {
                mtxcpy_size_t(best_medoids[k], medoids[k], dmatrixc,
                        medoids_card);
            }
            best_inst_adeq = cur_inst_adeq;
            best_inst = i;
        }
	}
    // Print the best configuration, confusion matrix and the CR
    // index.
	printf("\n");
    printf("Best adequacy %.15lf on instance %d.\n",
            best_inst_adeq, best_inst);
    printf("\n");
    print_medoids(best_medoids);
    printf("\n");
	print_memb(best_memb);
	printf("\n");
	print_weights(best_weights);
	printf("\n");
    global_energy();
    printf("\n");
    int *pred = defuz(memb, objc, clustc);
    print_groups(pred, objc, clustc);
    double **confmtx = confusion(pred, labels, objc);
    printf("\nConfusion matrix (class x predicted):\n");
    print_mtx_d(confmtx, classc, classc, 0);
    ++classc;
    for(i = 0; i < classc; ++i) {
        free(confmtx[i]);
    }
    free(confmtx);
    printf("Corrected Rand: %.7lf\n", corand(pred, labels, objc));
    free(pred);
    // Freeing memory.
END:
    fclose(stdout);
	for(i = 0; i < dmatrixc; ++i) {
		for(j = 0; j < objc; ++j) {
			free(dmatrix[i][j]);
		}
		free(dmatrix[i]);
	}
	free(dmatrix);
	for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            free(medoids[k][j]);
            free(best_medoids[k][j]);
        }
		free(medoids[k]);
		free(best_medoids[k]);
		free(weights[k]);
		free(best_weights[k]);
	}
	free(medoids);
	free(best_medoids);
	free(weights);
	free(best_weights);
	for(i = 0; i < objc; ++i) {
		free(memb[i]);
		free(best_memb[i]);
	}
	free(memb);
	free(best_memb);
	free(parc_cluster_adeq);
    free(parc_obj_adeq);
	return 0;
}
