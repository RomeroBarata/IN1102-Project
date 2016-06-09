#ifndef _STEX_H_
#define _STEX_H_

#pragma once

int* defuz(double **fuzmtx, size_t nrow, size_t ncol);

double** confusion(int *labels, int *pred, size_t size);

double corand(int *labels, int *pred, size_t size);

void print_groups(int *labels, size_t size, size_t card);

#endif /* _STEX_H_ */
