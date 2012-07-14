#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define WORD int
#define WLEN (sizeof(WORD)*8)

WORD bits_to_word(int *bits) {
    WORD result = 0;
    WORD mask = 1;
    int i;
    printf("bits to word\n");
    for (i = 0; i < WLEN; i ++) {
        printf("%d ", bits[i]);
        if (bits[i]) {
            result += mask;
        }
        mask <<= 1;
    }
    printf(" res:%d\n", result);
    return result;
}

void print_word(WORD word) {
    WORD mask = 1;
    int i;
    for (i = 0; i < WLEN; i ++) {
        printf("%d " , (word & mask)?1:0);
        mask <<= 1;
    }
    printf("\n");

}


int count_ones(WORD *seed, int seed_length) {
    int ones = 0;
    int k;
    WORD num;
    for (k = 0; k < seed_length; k++) {
        num = seed[k];
        while (num) {
            num &= num - 1;
            ones ++;
        }
    }
  return ones;
}

void process_chunk(WORD **chunk, int chunk_size, int total_words, int min_inds_per_bicluster){
    WORD *seed = (WORD *)malloc(total_words*sizeof(WORD));
    int i, j, k;
    for (i = 0; i < chunk_size; i ++) {
        for (j = 0; j < chunk_size; j ++) {
            for (k = 0; k < total_words; k++) {
                seed[k] = chunk[i][k] & chunk[j][k];
            }

            if (count_ones(seed, total_words) >= min_inds_per_bicluster) {


            }

        }

    }
    free(seed);
}

int *find_potential_snps (int **matrix, int total_snps, int total_inds) {
    int i, j;

    int total_words = 1 + ((total_inds -1) / WLEN);

    int *bits = (int *)calloc(WLEN, sizeof(int));

    int *potential_snps_table = (int *)calloc(total_snps, sizeof(int));
    int total_potential_snps = 0;

    WORD **bitmat = (WORD **)malloc(total_snps*sizeof(WORD *));

    printf("wlen: %d\n", WLEN);
    printf("total_words: %d\n", total_words);

    for (i = 0; i < total_snps; i++){

        bitmat[i] = (WORD *)calloc(total_words, sizeof(WORD));

        printf("%d\n",i);

        int word_index;
        for (word_index = 0; word_index < total_words; word_index++) {
            for (j = word_index*WLEN; j < (word_index + 1)*WLEN && j < total_inds; j++) {
                bits[j] = matrix[i][j];
            }
            while (!(j % WLEN)){
                bits[j] = 0;
            }
            printf("wi: %d\n", word_index);
            bitmat[i][word_index] = bits_to_word(bits);
        }
    }

    int chunk_size = 10000;
    int total_chunks = 1 + (total_snps - 1)/chunk_size;
    int chunk_index;



    for (chunk_index = 0; chunk_index < total_chunks; chunk_index ++) {
        for (i = chunk_index*chunk_size; i < (chunk_index+1)*chunk_size && i < total_snps; i ++) {
            for (j = chunk_index*chunk_size; j < (chunk_index+1)*chunk_size && j < total_snps; j ++) {



            }

        }

    }

//
//    while (c_chunk < total_chunks){
//        for (i = 0, j = c_chunk*chunk_size; i < chunk_size; i++, j++) {
//            for (int k = 0; k )
//            bitmat[i]
//
//        }
//
//    }

    for (i = 0; i < total_snps; i++){
        free(bitmat[i]);
    }
    free(bitmat);

    int *potential_snps = (int *)malloc(total_potential_snps * sizeof(int));
    for (i = 0, j = 0; i < total_snps; i ++) {

        if (potential_snps_table[i]) {
            potential_snps[j] = i;
            j++;
        }
    }

    free(potential_snps_table);

    return potential_snps;
}

int main(){
    int total_inds = 5;
    int total_snps = 5;
    srand(time(NULL));
    float MAF = 0.5;

    int **test = (int **)malloc(total_snps*sizeof(int *));
    int i, j;
    for (i = 0; i < total_snps; i ++) {
        test[i] = (int *)calloc(total_inds, sizeof(int));
        for (j = 0; j < total_inds; j ++) {
            test[i][j] = (((float)rand()/(float)RAND_MAX) < MAF )? 1 : 0;
            printf("%d ", test[i][j]);
        }
        printf("\n");
    }

    int *r = find_potential_snps(test, total_snps, total_inds);

    return 0;
}