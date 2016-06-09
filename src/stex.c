#include <stdlib.h>
#include <stdio.h>

#include "util.h"

// Defuzzifies fuzzy matrix 'fuzmtx' by using the first
// maxima method.
// Params:
//  fuzmtx - an object x class fuzzy matrix.
// Return:
//  An array of size 'fuzmtx->nrow' with the crip values for each
//  object.
int* defuz(double **fuzmtx, size_t nrow, size_t ncol) {
    size_t i;
    size_t j;
    double maxval;
    double val;
    int *labels = malloc(sizeof(int) * nrow);
    for(i = 0; i < nrow; ++i) {
        maxval = fuzmtx[i][0];
        labels[i] = 0;
        for(j = 1; j < ncol; ++j) {
            val = fuzmtx[i][j];
            if(val > maxval) {
                maxval = val;
                labels[i] = j;
            }
        }
    }
    return labels;
}

// Creates a confusion matrix.
// Params:
//  labels - the labels of each object.
//  pred - the predicted labels of each object.
//  size - the size of the 'labels' and 'pred' arrays.
// Return:
//  A confusion matrix.
double** confusion(int *labels, int *pred, size_t size) {
    int classc = max(labels, size) + 2;
	double **confmtx = malloc(sizeof(double *) * classc);
    size_t i;
	for(i = 0; i < classc; ++i) {
		confmtx[i] = calloc(sizeof(double), classc);
	}
    size_t last = classc - 1;
    for(i = 0; i < size; ++i) {
        ++(confmtx[labels[i]][pred[i]]);
        ++(confmtx[pred[i]][labels[i]]);
        ++(confmtx[labels[i]][last]);
        ++(confmtx[last][pred[i]]);
    }
    confmtx[last][last] = size;
    return confmtx;
}

// Computes the corrected Rand index.
// Params:
//  labels - the labels of each object.
//  pred - the predicted labels of each object.
//  size - the size of the 'labels' and 'pred' arrays.
// Return:
//  The corrected Rand index.
double corand(int *labels, int *pred, size_t size) {
    size_t i;
    size_t j;
    size_t last = size - 1;
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double d = 0.0;
    double val_labels;
    double val_pred;
    double val_labels_neg;
    double val_pred_neg;
    for(i = 0; i < last; ++i) {
        for(j = i + 1; j < size; ++j) {
            val_labels = labels[i] == labels[j];
            val_pred = pred[i] == pred[j];
            val_labels_neg = 1.0 - val_labels;
            val_pred_neg = 1.0 - val_pred;
            a += val_labels * val_pred;
            b += val_labels_neg * val_pred;
            c += val_labels * val_pred_neg;
            d += val_labels_neg * val_pred_neg;
        }
    }
    double p = a + b + c + d;
    double term = ((a + b) * (a + c) + (c + d) * (b + d)) * (1.0 / p);
    return ((a + d) - term) / (p - term);
}

// Prints the objects by group.
// Params:
//  labels - the labels of each object.
//  size - the size of the 'labels' and 'pred' arrays.
//  card - the number of different groups in 'labels'.
void print_groups(int *labels, size_t size, size_t card) {
    size_t groups[card][size - 1];
    size_t groupsc[card];
    size_t i;
    for(i = 0; i < card; ++i) {
        groupsc[i] = 0;
    }
    int label;
    for(i = 0; i < size; ++i) {
        label = labels[i];
        groups[label][groupsc[label]] = i;
        groupsc[label] = groupsc[label] + 1;
    }
    size_t j;
    for(i = 0; i < card; ++i) {
        printf("Group %u (%u members):\n", i, groupsc[i]);
        for(j = 0; j < groupsc[i]; ++j) {
            printf("%u ", groups[i][j]);
        }
        printf("\n");
    }
}

