#ifndef _UTIL_H_
#define _UTIL_H_

#pragma once

#include <stdbool.h>

#define FPOINT_OFFSET 1e-10

bool load_data(char *fname, double **matrix, size_t nrow, size_t ncol);

void print_mtx_d(double **matrix, size_t nrow, size_t ncol,
        size_t prec);

bool deq(double a, double b);

bool dgt(double a, double b);

bool dlt(double a, double b);

int max(int *vec, size_t size);

void mtxcpy_d(double **destination, double **source, size_t nrow,
        size_t ncol);

void mtxcpy_size_t(size_t **destination, size_t **source, size_t nrow,
        size_t ncol);

#endif /* _UTIL_H_ */
