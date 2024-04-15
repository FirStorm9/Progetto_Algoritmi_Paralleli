// mpic++ -o main-no-allegro main-no-allegro.cpp && mpirun -np 4 ./main-no-allegro

// or

// mpic++ -o main-no-allegro main-no-allegro.cpp  //in wsl: -pthread
// mpirun -np 4 ./main-no-allegro

#include <filesystem>
#include <mpi.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include "../Config.hpp"
#include <stdlib.h>
#include <time.h>
using namespace std;

const Config* cfg = new Config();

bool flag = false;

int* writeM = new int[(cfg->n_row/cfg->mpi_threads+2) * cfg->n_col];
int* readM = new int[(cfg->n_row/cfg->mpi_threads+2) * cfg->n_col];
int* combined_matrix = new int[(cfg->n_row) * cfg->n_col];

pthread_barrier_t barrier;

#define v(r,c) ((r)* cfg->n_col +(c))

#define EMPTY 0
#define PREY 1
#define PREDATOR 2

// time
double total_time, generation_time, draw_time, start_time, end_time, sum_total_time, communication_time = 0;

// MPI
MPI_Datatype row_type;
MPI_Comm cart_comm;
int my_rank, size_;
int upper_rank, lower_rank;

void swap_matrix(){
    int * m_temp = readM;
    readM = writeM;
    writeM = m_temp;
}

void print_matrix() {
    printf("Matrix rank%d!\n",my_rank);
    for (int i = 0; i < cfg->n_row / cfg->mpi_threads; i++) {
        for (int j = 0; j < cfg->n_col; j++) {
            printf("%d ", readM[v(i,j)]);
        }
        printf("\n");
    }
    printf("\n");
}


void create_matrix_in_file(){
    ofstream file_write("input.txt", ios::trunc);

    if (!file_write.is_open()) {
        cerr << "Unable to open file in writing mode" << endl;
        exit(1);
    }

    for (int i = 0; i < cfg->n_row / cfg->mpi_threads +2; i++) {
        for (int j = 0; j < cfg->n_col; j++) {
            int randNum = rand() % 100;
            if (randNum < 30) {
                file_write << PREY << " ";
            } else if (randNum < 40) {
                file_write << PREDATOR << " ";
            } else {
                file_write << EMPTY << " ";
            }
        }
        file_write << endl;
    }
    file_write.close();
}

void read_matrix_from_file(){
    ifstream file_read("input.txt");

    if (!file_read.is_open()) {
        cerr << "Unable to open file" << endl;
        exit(1);
    }

    // Lettura dei valori dalla matrice
    for (int i = 0; i < cfg->n_row / cfg->mpi_threads +2; i++) {
        for (int j = 0; j < cfg->n_col; j++) {
            if (!(file_read >> readM[v(i,j)])) {
                cerr << "Error in matrix reading" << endl;
                exit(1);
            }
        }
    }
    file_read.close();
}

void create_matrix() {

    switch(cfg->selection){
        case 0:
            for (int i = 0; i < cfg->n_row / cfg->mpi_threads +2; i++) {
                for (int j = 0; j < cfg->n_col; j++) {
                    if(i%2 != 0 && j%2 != 0){
                        readM[v(i,j)] = PREY;
                    }
                }
            }

            if(my_rank == cfg->mpi_root){
                for (int i = cfg->n_row / (cfg->mpi_threads*2); i < cfg->n_row/ cfg->mpi_threads; i++) {
                    for (int j = (cfg->n_col / cfg->mpi_threads) - 20; j < (cfg->n_col / cfg->mpi_threads) + 20; j++) {
                       if(i%3 != 0 && j%4 != 0){
                            readM[v(i,j)] = PREDATOR;
                        }
                    }
                }
            }
            break;
        case 1:
            create_matrix_in_file();
            read_matrix_from_file();
            break;
        case 2:
            read_matrix_from_file();
            break;
    }
}

void swap_row(int rank, MPI_Datatype row_type){

    MPI_Request request;
    MPI_Status status;

    MPI_Cart_shift(cart_comm, 0, 1, &upper_rank, &lower_rank); // direzione: 0, spostamento: 1

    const int dest = (rank == 0) ? cfg->mpi_threads-1 : 0;

    MPI_Isend(&readM[v(1, 0)], 1, row_type, upper_rank, 20, cart_comm, &request);
    MPI_Isend(&readM[v(cfg->n_row/cfg->mpi_threads, 0)], 1, row_type, lower_rank, 1, cart_comm, &request);

    MPI_Recv(&readM[v(cfg->n_row/cfg->mpi_threads+1,0)], 1, row_type, lower_rank, 20, cart_comm, &status);
    MPI_Recv(&readM[v(0, 0)], 1, row_type, upper_rank, 1, cart_comm, &status);
}

void neighbourt(int x, int y) {
    int count_prey = 0;
    int count_predator = 0;

    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            if (di == 0 && dj == 0)
                continue;

            int neighbor_x = x + di;
            int neighbor_y = y + dj;

            if (neighbor_x < 0){
                if(readM[v(cfg->n_row-1, neighbor_y)] == PREY){
                    count_prey++;
                }else if(readM[v(cfg->n_row-1, neighbor_y)] == PREDATOR){
                    count_predator++;
                }
            }
            if(neighbor_x >= cfg->n_row){
                if(readM[v(0, neighbor_y)] == PREY){
                    count_prey++;
                }else if(readM[v(0, neighbor_y)] == PREDATOR){
                    count_predator++;
                }
            }

            if(neighbor_y < 0){
                if(readM[v(neighbor_x, cfg->n_col-1)] == PREY){
                    count_prey++;
                } else if(readM[v(neighbor_x, cfg->n_col-1)] == PREDATOR){
                    count_predator++;
                }
            }
            else if(neighbor_y > cfg->n_col-1){
                if(readM[v(neighbor_x, 0)] == PREY){
                    count_prey++;
                } else if(readM[v(neighbor_x, 0)] == PREDATOR){
                    count_predator++;
                }
            }
            else{
                if (readM[v(neighbor_x, neighbor_y)] == PREY){
                    count_prey++;
                }else if (readM[v(neighbor_x, neighbor_y)] == PREDATOR){
                    count_predator++;
                }
            }

            int current_state = readM[v(x, y)];


            if (current_state == EMPTY){
                if(count_prey >= 2 && rand() % 100 < 25){
                    writeM[v(x,y)] = PREY;
                }
                else if(count_predator >= 2 && count_prey >= 1 ){
                    writeM[v(x,y)] = PREDATOR;
                }
                else{
                    writeM[v(x,y)] = EMPTY;
                }
            }
            else if(current_state == PREY){
                if(count_predator > 0){
                    writeM[v(x,y)] = EMPTY;
                }
                else {
                    writeM[v(x,y)] = PREY;
                }
            }
            else if(current_state == PREDATOR){
                if(count_prey == 0 || count_predator >= 3){
                    writeM[v(x,y)] = EMPTY;
                }
                else{
                    writeM[v(x,y)] = PREDATOR;
                }
            }
        }

    }
}

void* game(void* arg){

    int* p = (int*)arg;
    int rank_posix = *p;


    const int start_position = rank_posix * (cfg->n_col / cfg->posix_threads);
    int end_position = start_position + (cfg->n_col / cfg->posix_threads);

    // Se il numero non di colonne non è divisibile per il numero di posix_thread, aggiungi i restanti
    if(rank_posix == cfg->posix_threads-1){
        end_position += cfg->n_row % cfg->posix_threads;
    }

    for (int i = 1; i < cfg->n_row / cfg->mpi_threads+1; ++i) {
        for (int j = start_position; j < end_position; ++j) {

            neighbourt(i, j);

        }
    }

    pthread_barrier_wait(&barrier);
    
    delete p;
    return NULL;
}



void posix_threads_handler(){

    pthread_barrier_init(&barrier, NULL, cfg->posix_threads);

    //Sezione posix
    pthread_t tid[cfg->posix_threads]; // tid = thread id

    for (int i = 0; i < cfg->posix_threads; i++) {
        int* p = new int;
        *p = i;
        pthread_create(&tid[i], NULL, &game, p);
    }

    for (int i = 0; i < cfg->posix_threads; i++) {
        pthread_join(tid[i], NULL);
    }
    pthread_barrier_destroy(&barrier);
}


void write_header(std::ofstream& file_write) {

    std::string separator = ",  ";
    file_write << "total_time" << separator
        << "communication_time" << separator
        << "generation_time" << separator
        << "draw_time" << separator
        << "start_time" << separator
        << "end_time"<< separator
        << "matrix_size"
        << "\n";
}

void write_result(std::ofstream& file_write) {
    
    std::string separator = ",  ";
    file_write << total_time << separator
        << communication_time << separator
        << generation_time << separator
        << draw_time << separator
        << start_time << separator
        << end_time << separator
        << cfg->n_row << "x" << cfg->n_col
        <<"\n";

}

void write_times_on_file(){

    ofstream file;
    string file_path = "./result.txt";
    if(std::filesystem::exists(file_path)) {
        file.open(file_path, std::ios::app);
    }
    else {
        std::cout << "Creating results file" << std::endl;
        file.open(file_path);
        write_header(file);
    }
    write_result(file);
    file.close();

    std::cout << "Results written inside " << file_path << std::endl;

}

void finalize(){

    
    MPI_Type_free(&row_type);
    MPI_Comm_free(&cart_comm);
    
    MPI_Finalize();

    delete[] readM;
    delete[] writeM;
    delete[] combined_matrix;
    delete cfg;

}


int main(int argc, char *argv[]) {
    const auto start = std::chrono::high_resolution_clock::now();
    int media[cfg->mpi_threads];

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    // sezione topologia virtuale
    const int ndims = 1;
    const int dim_size[1] = {cfg->mpi_threads};
    const int periods[1] = {1};
    const int reorder = 0;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &cart_comm);

    MPI_Type_contiguous(cfg->n_col, MPI_INT, &row_type);
    MPI_Type_commit(&row_type);

    create_matrix();

    MPI_Barrier(MPI_COMM_WORLD);

    int i = 0;
    while(i < cfg->steps){
        start_time = MPI_Wtime();
        posix_threads_handler();
        swap_matrix();
        swap_row(my_rank, row_type);

        MPI_Gather(&readM[v(1,0)], (cfg->n_row/cfg->mpi_threads) * cfg->n_col, MPI_INT, &combined_matrix[v(0,0)], (cfg->n_row/cfg->mpi_threads) * cfg->n_col, MPI_INT, 0, cart_comm);

        if (my_rank == cfg->mpi_root){
            //print_matrix();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        i++;

        // calculate time
        end_time = MPI_Wtime();
        communication_time = end_time - start_time;
        total_time = MPI_Wtime() - start_time;
        sum_total_time = 0.0;
        MPI_Reduce(&total_time, &sum_total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        double avg_time = sum_total_time / size_;
        if (!flag){
            if (my_rank == cfg->mpi_root) {
                printf("Average communication time at step %d: %lf seconds\n", i, avg_time);
            }
            if(my_rank == cfg->mpi_root){
                printf("Communication time on step %d, is: %lf\n", i, communication_time);
            }
            flag=true;
        }
    }

    const auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = stop - start;
    if(my_rank == cfg->mpi_root){
        cout << "Time taken for the computation of " << cfg->steps <<" generations, with matrix size equal to " <<cfg->n_row <<"x"<< cfg->n_row<< ": " <<elapsed_time.count()<<" seconds"<<endl;
    }

    if(my_rank == cfg->mpi_root){
        write_times_on_file();
    }

    finalize();

    return 0;
}
