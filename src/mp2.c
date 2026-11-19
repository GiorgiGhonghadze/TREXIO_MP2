#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <trexio.h>

//ERI OFFSET
static int32_t eri_offset(int32_t mo_num, int32_t p, int32_t q, int32_t r, int32_t s){
    int32_t n = (int32_t)mo_num;
    return ((p * n + q) * n + r) * n + s;
}

//STORE 8-FOLD SYMMETRIC ERIs
static void store_eri(double* eri, int32_t mo_num, int32_t i, int32_t j, int32_t k, int32_t l, double val){
    eri[eri_offset(mo_num, i, j, k, l)] = val;
    eri[eri_offset(mo_num, i, l, k, j)] = val;
    eri[eri_offset(mo_num, k, j, i, l)] = val;
    eri[eri_offset(mo_num, k, l, i, j)] = val;
    eri[eri_offset(mo_num, j, i, l, k)] = val;
    eri[eri_offset(mo_num, l, i, j, k)] = val;
    eri[eri_offset(mo_num, l, k, j, i)] = val;
    eri[eri_offset(mo_num, j, k, l, i)] = val;
}

int main(int argc, char* argv[])
{
    // CHECK THE PROGRAM ARGUMENTS
    if (argc < 2) {
        printf( "Please Provide Trexio .h5 File\n");
        printf( "Usage: ./program trexio_file_path\n");
        exit(1);
    }

    //OPEN THE TREXIO FILE
    const char* filename = argv[1];
    trexio_exit_code rc;
    trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }

    //READ E_nn, MO_num and N_occ 
    double e_nn = 0.0;
    int32_t n_occ = 0;
    int32_t mo_num = 0;

    rc = trexio_read_nucleus_repulsion(trexio_file, &e_nn);
    rc = trexio_read_electron_up_num(trexio_file, &n_occ);
    rc = trexio_read_mo_num(trexio_file, &mo_num);


    //ALLLOCATE PROPERLY SIZED VARIABLES ON THE HEAP, READ ORBITAL ENERGIES AND FLATTENED ONE-ELECTRON INTEGRALS
    double* mo_energy = malloc((double)mo_num * sizeof(double));
    double* h_core = malloc((double)mo_num * (double)mo_num * sizeof(double));

    rc = trexio_read_mo_energy(trexio_file, mo_energy);
    rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, h_core);

    //READ UNIQUE ERIs 
    int64_t n_integrals = 0;
    rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);

    int32_t* const index = malloc(4 * n_integrals * sizeof(int32_t));
    if (index == NULL) {
        fprintf(stderr, "Malloc failed for index");
        exit(1);
    }
    double* const value = malloc(n_integrals * sizeof(double));
    if (value == NULL) {
        fprintf(stderr, "Malloc failed for value");
        exit(1);
    }

    int64_t buffer_size = n_integrals;
    rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &buffer_size, index, value);
    
    //CLOSE TREXIO FILE
    rc = trexio_close(trexio_file);
    trexio_file = NULL;

    //STORE ALL 8-FOLD SYMMETRIC INTEGRALS
    int32_t total_eri_elements = pow(mo_num, 4);
    double* eri = calloc(total_eri_elements, sizeof(double));

    for (int64_t n = 0; n < buffer_size; n++) {

        int32_t i = index[4 * n + 0];
        int32_t j = index[4 * n + 1];
        int32_t k = index[4 * n + 2];
        int32_t l = index[4 * n + 3];
        store_eri(eri, mo_num, i, j, k, l, value[n]);

    }

    free(index);
    free(value);

    //CALCULATE HF ENERGY
    double e_1e = 0.0;
    for (int32_t i = 0; i < n_occ; i++) {

        e_1e += 2.0 * h_core[i * mo_num + i];

    }

    double e_2e = 0.0;
    for (int32_t i = 0; i < n_occ; i++) {
        for (int32_t j = 0; j < n_occ; j++) {

            double j_ij = eri[eri_offset(mo_num, i, j, i, j)];
            double k_ij = eri[eri_offset(mo_num, i, j, j, i)];
            e_2e += 2.0 * j_ij - k_ij;

        }
    }

    double e_hf = e_nn + e_1e + e_2e;

    //CALCULATE MP2 CORRELATION ENERGY
    double e_mp2 = 0.0;
    for (int32_t i = 0; i < n_occ; i++) {
        for (int32_t j = 0; j < n_occ; j++) {
            for (int32_t a = n_occ; a < mo_num; a++) {
                for (int32_t b = n_occ; b < mo_num; b++) {

                    double ijab = eri[eri_offset(mo_num, i, j, a, b)];
                    double ijba = eri[eri_offset(mo_num, i, j, b, a)];
                    double denom = mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b];
                    e_mp2 += (ijab * (2.0 * ijab - ijba)) / denom;

                }
            }
        }
    }

    //PRINT THE RESULTS
    printf("Number of MOs: %d\n", mo_num);
    printf("Number of occupied MOs: %d\n", n_occ);
    printf("\nHartree-Fock energy:       %.10f a.u.\n", e_hf);
    printf("MP2 correlation energy:     %.10f a.u.\n", e_mp2);
    printf("Total MP2 energy:          %.10f a.u.\n", e_hf + e_mp2);


    free(eri);
    free(mo_energy);
    free(h_core);

    return 0;
}