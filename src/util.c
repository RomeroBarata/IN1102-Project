#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "util.h"

bool load_data(char *fname, double **matrix, size_t nrow,
        size_t ncol) {
	FILE *ifile = fopen(fname, "r");
	if(!ifile) {
		return false;
	}
	size_t i;
	size_t j;
	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < ncol; ++j) {
			if(fscanf(ifile, "%lf", &matrix[i][j]) == EOF) {
				fclose(ifile);
				return false;
			}
		}
	}
	fclose(ifile);
	return true;
}

void print_mtx_d(double **matrix, size_t nrow, size_t ncol,
        size_t prec) {
	size_t i;
	size_t j;
	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < ncol - 1; ++j) {
			printf("%5.*lf ", prec, matrix[i][j]);
		}
		printf("%5.*lf\n", prec, matrix[i][j]);
	}
}

bool deq(double a, double b) {
    return (a < (b + FPOINT_OFFSET) && a > (b - FPOINT_OFFSET));
}

bool dgt(double a, double b) {
    return a > (b + FPOINT_OFFSET);
}

bool dlt(double a, double b) {
    return a < (b - FPOINT_OFFSET);
}

int max(int *vec, size_t size) {
    int ret;
    size_t i;
    for(i = 0; i < size; ++i) {
        if(!i || vec[i] > ret) {
            ret = vec[i];
        }
    }
    return ret;
}

void mtxcpy_d(double **destination, double **source, size_t nrow,
        size_t ncol) {
    size_t i;
    for(i = 0; i < nrow; ++i) {
        memcpy(destination[i], source[i], sizeof(double) * ncol);
    }
}

void mtxcpy_size_t(size_t **destination, size_t **source, size_t nrow,
        size_t ncol) {
    size_t i;
    for(i = 0; i < nrow; ++i) {
        memcpy(destination[i], source[i], sizeof(size_t) * ncol);
    }
}
