// c++ -o serial-project-no-allegro serial-project-no-allegro.cpp && ./serial-project-no-allegro

// or 

// c++ -o serial-project-no-allegro serial-project-no-allegro.cpp
// ./serial-project-no-allegro

#include <filesystem>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include "../config.cfg"
#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>
using namespace std;

const Config* cfg = new Config();

bool flag = false;

int* writeM = new int[(cfg->n_row) * cfg->n_col];
int* readM = new int[(cfg->n_row) * cfg->n_col];

#define v(r,c) ((r)* cfg->n_col +(c)) //accesso alla cella

#define EMPTY 0
#define PREY 1
#define PREDATOR 2

// time
std::chrono::duration<double> elapsed;

void swap_matrix() {
    int * m_temp = readM;
    readM = writeM;
    writeM = m_temp;
}

void print_matrix() {
    for (int i = 0; i < cfg->n_row; i++) {
        for (int j = 0; j < cfg->n_col; j++) {
            printf("%d ", readM[v(i,j)]);
        }
        printf("\n");
    }
    printf("\n");
    system("clear");
}

void create_matrix_in_file() {
    ofstream file_write("input.txt", ios::trunc);

    if (!file_write.is_open()) {
        cerr << "Unable to open file in writing mode" << endl;
        exit(1);
    }

    for (int i = 0; i < cfg->n_row; i++) {
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

void read_matrix_from_file() {
    ifstream file_read("input.txt");

    if (!file_read.is_open()) {
        cerr << "Unable to open file" << endl;
        exit(1);
    }

    // Lettura dei valori dalla matrice
    for (int i = 0; i < cfg->n_row; i++) {
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
            for (int i = 0; i < cfg->n_row; i++) {
                for (int j = 0; j < cfg->n_col; j++) {
                    if(i%2 != 0 && j%2 != 0){
                        readM[v(i,j)] = PREY;
                    }else if(i%3 != 0 && j%4 != 0 && i%2 != 0){
                        readM[v(i,j)] = PREDATOR;
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
            } else if(neighbor_y > cfg->n_col-1){
                if(readM[v(neighbor_x, 0)] == PREY){
                    count_prey++;
                } else if(readM[v(neighbor_x, 0)] == PREDATOR){
                    count_predator++;
                }
            } else{
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
                } else if(count_predator >= 2 && count_prey >= 1 ){ //&& rand() % 100 < 60
                    writeM[v(x,y)] = PREDATOR;
                } else{
                    writeM[v(x,y)] = EMPTY;
                }
            } else if(current_state == PREY){
                if(count_predator > 0){
                    writeM[v(x,y)] = EMPTY;
                } else {
                    writeM[v(x,y)] = PREY;
                }
            } else if(current_state == PREDATOR){
                if(count_prey == 0 || count_predator >= 3){
                    writeM[v(x,y)] = EMPTY;
                } else{
                    writeM[v(x,y)] = PREDATOR;
                }
            }
        }
    }
}

void* game(){
    for (int i = 1; i < cfg->n_row; ++i) {
        for (int j = 0; j < cfg->n_row; ++j) {
            neighbourt(i, j);
        }
    }
    return NULL;
}

void write_header(std::ofstream& file_write) {

    std::string separator = ",  ";
    file_write
        << "elapsed" << separator
        << "mpi_threads" << separator
        << "posix_threads" << separator
        << "steps"<< separator
        << "matrix_size"
        << "\n";
}

void write_result(std::ofstream& file_write) {
    
    std::string separator = ",  ";
    file_write
        << elapsed.count() << separator
        << "1" << separator
        << "0" << separator
        << cfg->steps << separator 
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

    delete[] readM;
    delete[] writeM;
    delete cfg;
}

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    create_matrix();

    int i = 0;
    while(i < cfg->steps){
        game();
        swap_matrix();
        //print_matrix();
        i++;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    elapsed = (stop - start);
    cout << "Time taken for the computation of  " << cfg->steps << " generations, with matrix size equal to " << cfg->n_row << "x" << cfg->n_row << ": " << elapsed.count()<<" seconds"<<endl;

    write_times_on_file();

    finalize();

    return 0;
}

