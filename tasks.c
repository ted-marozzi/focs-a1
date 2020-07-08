/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 910193
 *   Name        : Edward Phillip Shugg Marozzi
 *
 ***************************************************************************/
/* Compiled with: gcc -Wall -std=c99 *.c -o flow -lm -g */
/* Executed on windows with: flow flow_data.csv 10 */

/*******************************************************************************
 * Standard includes
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include "tasks.h"

/*******************************************************************************
 * Global constants
*******************************************************************************/
#define SEEK_SET 0
#define OFFSET 0
/* Task one defines */
#define X_THRESHOLD 20
#define TASK_1_ROWS 4
#define HEADER_CHARS 12
#define TASK_1 "task1.csv"
#define ROW_1 0
#define ROW_2 1
#define ROW_3 2
#define ROW_4 3

/* Task two defines */
#define TASK_2 "task2.csv"
#define X_LOWER_LIMIT -15
#define X_UPPER_LIMIT 85
#define Y_LOWER_LIMIT -20
#define Y_UPPER_LIMIT 20

/* Task three defines */
#define TASK_3 "task3.csv"
#define TRUE 1
#define FALSE 0
/* Binary Search */
#define BST_SUCCESS 1
#define BST_FAILURE 0
#define BST_PREORDER 0
#define BST_INORDER 1
#define BST_POSTORDER 2

/* Task four defines */
#define THRESHOLD_0 0
#define THRESHOLD_1 5
#define THRESHOLD_2 10
#define THRESHOLD_3 15
#define THRESHOLD_4 20
#define THRESHOLD_5 25
#define NUM_THRESHOLDS 5


/*******************************************************************************
 * Structures and type definitions
*******************************************************************************/

/* Row of data struct */
typedef struct 
{ 
    double rho, u, v, x, y, S; 
} row_t;

/* Cell of data struct*/
typedef struct 
{ 
    /* Cell contains one row with */
    row_t average_row;
    /* Meta data */
    unsigned int num_cell_points;

} cell_t;

/* Grid Struct */
typedef struct
{
    /* Dimensions */
    int n, m;
    long unsigned int total_lines;

} grid_t;


/* Binary search tree list_node type */
typedef struct bst_flux_node bst_flux_node_t;

/* Binary Search Node */
struct bst_flux_node {
    void* data;
    bst_flux_node_t* left;
    bst_flux_node_t* right;
};

/* Binary search tree type */
typedef struct {
    int num_elements;
    bst_flux_node_t* root;
    void (*del)(void*);
   
} bst_flux_t;

/* Linked list list node type */
typedef struct list_node l_node_t;

struct list_node {
    void* data;
    l_node_t* next;
};

/* linked list type */
typedef struct {
    int num_elements;
    l_node_t* head;
    l_node_t* tail;
    void (*del)(void*);
} list_t;

/* Used for program execution timing */
struct timeval start;
struct timeval stop;

/*******************************************************************************
 * Function headers
*******************************************************************************/

/* Function prototypes */
FILE* read_in_file(const char* file);
int row_scanner(FILE* input_file, row_t* row);

/* Task one functions */
void parse_data_1(FILE* input_file, row_t* max_flux_arr);
void max_min_flux_calc(row_t* row, row_t* max_flux_arr);
void write_to_task_1(char* headers, row_t* max_flux_arr);

/* Task two functions */
void zero_cell_points(int resolution, cell_t cell_arr[][resolution]);
void parse_data_2(FILE* task_2_file, row_t row, int resolution, 
    cell_t cell_arr[][resolution], row_t* linear_row_arr);
int valid_x_point(int resolution, row_t* row, int x_point_bool, int* i);
int valid_y_point(int resolution, row_t* row, int y_point_bool, int* j);
void write_to_task_2(char* headers, row_t* linear_row_arr, FILE* task_2_file,
    int resolution);
void average_sum(int resolution, int x_point_bool, int y_point_bool, 
    cell_t cell_arr[][resolution], row_t row, int i, int j);
void average_divide(int resolution, int x_point_bool, int y_point_bool, 
    cell_t cell_arr[][10], row_t row, int i, int j);
void calc_s(int resolution, cell_t cell_arr[][resolution], int i, int j);
/* Sort functions */
int row_array_s_comp(const void* a, const void* b);
int row_array_flux_comp(const void* a, const void* b);
void merge_sort(void *base, size_t n_memb, size_t elem_size,
                int (*cmp)(const void *, const void *));
static void merge(void *out, const void *pa, size_t na,
                  const void *pb, size_t nb, size_t elem_size,
                  int (*cmp)(const void *, const void *));

/* Task three functions*/
int find_num_center(FILE* input_file, row_t row);
void create_center_arr(FILE* input_file, row_t row, row_t *center_points);
void create_flux_array(int num_center, row_t* center_points, double* u_flux);
double find_closest(double* u_flux, int num_center);
void linear_array_search(FILE* task3, int num_center, double* u_flux, 
    double closest);
/* Linked list */
list_t* list_new(void (*delfunc)(void*));
void list_push_back(list_t* row_list, void* d);
void* list_pop_front(list_t* row_list);
void list_process_all(l_node_t* p, void (*process)(void*));
void list_free(list_t* list);
void create_linked_list(int num_center, row_t* center_points, list_t* row_list);
void linear_list_search(FILE* task3, list_t* row_list, double closest);
/* Binary Tree */
int make_unique(double* array, int n);
void no_free(void* v);
void perfect_insert(bst_flux_t* bst_flux, double* array, int low, int high);
int bst_flux_insert(bst_flux_t* bst_flux, void* d);
void bst_flux_free(bst_flux_t* bst_flux);
void bst_flux_free_subtree(bst_flux_t* bst_flux, bst_flux_node_t* n);
bst_flux_t* bst_flux_new(void (*delfunc)(void*));
/* Binary Search */
int binary_search(double arr[], int l, int r, double x, FILE* task3);
/* Binary search on array */
bst_flux_node_t* bst_flux_find(bst_flux_t* bst_flux, double d, FILE* task3);

/* Task four functions */
void calc_n_m(FILE* input_file, grid_t* grid);


/*******************************************************************************
 * Task one: Calculate the maximum flux difference in the data set 
*******************************************************************************/
void maxfluxdiff(const char* flow_file)
{
    FILE *input_file;
    row_t *max_flux_arr = {0};
    char *headers;

    /* Dynamically assigns the memory for the first row */
    headers = (char*)malloc(HEADER_CHARS * sizeof(char));

    /* Creates file object */
    input_file = read_in_file(flow_file);
    
  
    /* Dynamically defines an array to store the maximum flux difference data */
    max_flux_arr = calloc(TASK_1_ROWS, sizeof(row_t));

    /* Return check for memory allocation */
    if(max_flux_arr == NULL || headers == NULL)
    {
        printf("Memory can't be allocated\n");
        exit(EXIT_SUCCESS);
    }
    
    /* Scans in the row of attibutes */
    fscanf(input_file,"%s", headers);

    /* Read data in while possible */
    parse_data_1(input_file, max_flux_arr);

    /* Close the file */
    fclose(input_file);
    
    /* Close */
    write_to_task_1(headers, max_flux_arr);

    /* Free the data */
    free(max_flux_arr);
    free(headers);
}

/* Reads in the file */
FILE* read_in_file(const char* file)
{
    
    FILE *input_file;
    input_file = fopen(file, "r");

    /* Couldn't find file */
    if(input_file == NULL)
    {
        printf("File couldn't be opened");
        exit(EXIT_SUCCESS);
    }
    /* Return the file object */
    return input_file;
}

/* Passes data for min max to be calcuated on */
void parse_data_1(FILE* input_file, row_t* max_flux_arr)
{
    /* Task one */
    row_t row = {0};
    int nums_scanned = 0;

    /* EOF is the end on file */
    while(nums_scanned != EOF)
    {
        /* Scans in the row */
        nums_scanned = row_scanner(input_file, &row);
        
        /* Calculates the min max calc */
        max_min_flux_calc(&row, max_flux_arr);
    }
    
}

/* Scans in rows */
int row_scanner(FILE* input_file, row_t* row)
{   
    int nums_scanned;
    /* Nums scanned is the number of valid reads */
    nums_scanned = fscanf(input_file,"%lf,%lf,%lf,%lf,%lf", 
                    &row->rho, &row->u, &row->v, &row->x, &row->y);

    return nums_scanned;
}

/* Extracts the min and max */
void max_min_flux_calc(row_t* row, row_t* max_flux_arr)
{
    double flux_u, flux_v;

    /* Checks for valid x */
    if(row->x > X_THRESHOLD)
    {
        /* Reduces the length of calls */
        flux_u = row->rho*row->u;
        flux_v = row->rho*row->v;

        /* Max u flux assignment */
        if(flux_u > max_flux_arr[ROW_1].rho*max_flux_arr[ROW_1].u)
        {
            max_flux_arr[ROW_1] = *row;
        }

        /* Min u flux assignment */
        if(flux_u <  max_flux_arr[ROW_2].rho*max_flux_arr[ROW_2].u || 
             max_flux_arr[ROW_2].rho*max_flux_arr[ROW_2].u == 0)
        {
            max_flux_arr[ROW_2] = *row;
        }

        /* Max v flux assignment */
        if(flux_v > max_flux_arr[ROW_3].rho*max_flux_arr[ROW_3].v)
        {
            max_flux_arr[ROW_3] = *row;
        }

        /* Min v flux assignment */
        if((flux_v < max_flux_arr[ROW_4].rho*max_flux_arr[ROW_4].v ) || 
            max_flux_arr[ROW_4].rho*max_flux_arr[ROW_4].v == 0)
        {
            max_flux_arr[ROW_4] = *row;
        }   
        
    }

}

/* Writes to task one */
void write_to_task_1(char* headers, row_t* max_flux_arr)
{
    FILE* task_1_file;
    task_1_file = fopen(TASK_1, "w");
    fprintf(task_1_file,"%s\n", headers);
    /* prints the 4 rows */
    for(int i = 0; i < TASK_1_ROWS; i++)
    {
        fprintf(task_1_file, "%lf,%lf,%lf,%lf,%lf\n", max_flux_arr[i].rho,  
            max_flux_arr[i].u, max_flux_arr[i].v, 
            max_flux_arr[i].x, max_flux_arr[i].y);

    }
    /* close the data off */
	fclose(task_1_file);
}

/*******************************************************************************
 * Task two: Calculate the the course grid
*******************************************************************************/

void coarsegrid(const char* flow_file, int resolution)
{
    FILE *task_2_file, *input_file;
    row_t row;

    //int i;
    char *headers;

    // Dynamically assigns the memory for the first row
    headers = (char*)malloc(HEADER_CHARS * sizeof(char));
    row_t* linear_row_arr = malloc(sizeof(row_t)*resolution*resolution);
    /* 2d cell allocation array */
    cell_t (*cell_arr)[resolution] = 
        malloc(sizeof(cell_t[resolution][resolution]));

    /* Initialises points to zero */
    zero_cell_points(resolution, cell_arr);

    input_file = read_in_file(flow_file);

    /* Scans in the row of attibutes */
    fscanf(input_file,"%s", headers);
    
    /* Pass data to calulate the grid */
    parse_data_2(input_file, row, resolution, cell_arr, linear_row_arr);

    /* Open file to write too */
    task_2_file = fopen("task2.csv", "w");
    
    /* Merge sort the array */
    merge_sort(linear_row_arr, resolution*resolution, sizeof(row_t), 
        row_array_s_comp);

    /* Write results to file */
    write_to_task_2(headers, linear_row_arr, task_2_file, resolution);

    /* Clean up */
    fclose(input_file);
    fclose(task_2_file);

    free(linear_row_arr);
    free(headers);
    free(cell_arr);
}

/* Initalise all cell points */
void zero_cell_points(int resolution, cell_t cell_arr[][resolution])
{
    int i, j;
    /* Loop through every point */
    for(i=0;i<resolution;i++)
    {
        for(j=0;j<resolution;j++)
        {
            cell_arr[i][j].num_cell_points = 0;
            cell_arr[i][j].average_row.rho = 0;
            cell_arr[i][j].average_row.u = 0;
            cell_arr[i][j].average_row.v = 0;
            cell_arr[i][j].average_row.x = 0;
            cell_arr[i][j].average_row.y = 0;
        }
    }
}

/* Gives data to average sum function and average divide function */
void parse_data_2(FILE* task_2_file, row_t row, int resolution, 
                    cell_t cell_arr[][resolution], row_t* linear_row_arr)
{
    int x_point_bool, y_point_bool, i, j, l = 0, nums_scanned = 0;

    
    while(nums_scanned != EOF)
    {   
        nums_scanned = row_scanner(task_2_file, &row);
        
        x_point_bool = FALSE;
        y_point_bool = FALSE;

        /* Finds the i and j co-ord of where the row exists */
        x_point_bool = valid_x_point(resolution, &row, x_point_bool, &i);
        y_point_bool = valid_y_point(resolution, &row, y_point_bool, &j);

        /* Adds the data point to the sum */
        average_sum(resolution, x_point_bool, y_point_bool, cell_arr, row, i, j);

    }

    
    for(j = 0; j < resolution; j++)
    {
        for(i =0 ; i< resolution; i++)
        {
            /* Divides the sum by the number of points */
            average_divide(resolution, x_point_bool, y_point_bool, cell_arr, 
                row, i, j);

            /* Calculates S */
            calc_s(resolution, cell_arr, i, j);

            /* Assigns the 2d array to a 1d array */
            linear_row_arr[l] = cell_arr[i][j].average_row;
            l++;
        }
    }
}

/* calculates S */
void calc_s(int resolution, cell_t cell_arr[][resolution], int i, int j)
{
    cell_arr[i][j].average_row.S = 100*sqrt(pow(cell_arr[i][j].average_row.u,2)
        +pow(cell_arr[i][j].average_row.v,2))/
        sqrt(pow(cell_arr[i][j].average_row.x,2)+
        pow(cell_arr[i][j].average_row.y,2));


}

/* Used in calculating the average cell points */
void average_divide(int resolution, int x_point_bool, int y_point_bool, 
    cell_t cell_arr[][resolution], row_t row, int i, int j)
{
    cell_arr[i][j].average_row.rho /= cell_arr[i][j].num_cell_points;
    cell_arr[i][j].average_row.u /= cell_arr[i][j].num_cell_points;
    cell_arr[i][j].average_row.v /= cell_arr[i][j].num_cell_points;
    cell_arr[i][j].average_row.x /= cell_arr[i][j].num_cell_points;
    cell_arr[i][j].average_row.y /= cell_arr[i][j].num_cell_points;

}

/* Used in calculating the average cell points */
void average_sum(int resolution, int x_point_bool, int y_point_bool, 
    cell_t cell_arr[][resolution], 
row_t row, int i, int j)
{
    /* Checks for valid point */
    if(x_point_bool && y_point_bool)
    {
        cell_arr[i][j].num_cell_points++;
        cell_arr[i][j].average_row.rho += row.rho;
        cell_arr[i][j].average_row.u += row.u;
        cell_arr[i][j].average_row.v += row.v;
        cell_arr[i][j].average_row.x += row.x;
        cell_arr[i][j].average_row.y += row.y;

    }
}

/* i will be set to the location of [x][y] of where the row is located */
int valid_x_point(int resolution, row_t* row, int x_point_bool, int* i)
{

    double x_range = abs(X_UPPER_LIMIT) + abs(X_LOWER_LIMIT);
    double x_step = x_range/resolution;
    /* Go through X domain */
    for(*i = 0; *i<resolution; (*i)++)
    {
        if(row->x> X_LOWER_LIMIT + (*i)*x_step && 
            row->x< X_LOWER_LIMIT+(*(i)+1)*x_step)
        {
            x_point_bool = TRUE;
            break;
        }
    }
    /* Checks if x point was valid in the domain */
    return x_point_bool;
}

/* j wll be set to the location of [x][y] of where the row is located */
int valid_y_point(int resolution, row_t* row, int y_point_bool, int *j)
{

    double y_range = abs(Y_UPPER_LIMIT) + abs(Y_LOWER_LIMIT);
    double y_step = y_range/resolution;

    for((*j) = 0; *j < resolution; (*j)++)
    {
        /* Go through Y domain */
        if(row->y> Y_LOWER_LIMIT + (*j)*y_step && 
            row->y< Y_LOWER_LIMIT+((*j)+1)*y_step)
        {
            
            y_point_bool = TRUE;
            break;
        }
    }
    /* Checks if y point was valid in the domain */
    return y_point_bool;
}

/* write out the results for task 2 */
void write_to_task_2(char* headers, row_t* linear_row_arr, FILE* task_2_file,
    int resolution)
{
    fprintf(task_2_file,"%s,S\n", headers);
    int i, j, k = resolution*resolution - 1;
    for(i = 0; i < resolution; i++)
    {
        for(j = 0; j< resolution; j++)
        {
            /* Write it out in reverse */
            fprintf(task_2_file, "%lf,%lf,%lf,%lf,%lf,%lf\n", 
                linear_row_arr[k].rho, linear_row_arr[k].u, 
                linear_row_arr[k].v, linear_row_arr[k].x, linear_row_arr[k].y,
                linear_row_arr[k].S);
            k--;
        }
    }
    
}

/* Compares the s data for merge sort */
int row_array_s_comp(const void* a, const void* b)
{
    row_t* ra = (row_t*)a;
    row_t* rb = (row_t*)b;
    /* With double data we can just subtract to get the right behaviour */
    return (ra->S > rb->S) - (ra->S < ra->S);
}

/* Compares the flux data for merge sort */
int row_array_flux_comp(const void* a, const void* b)
{
    row_t* ra = (row_t*)a;
    row_t* rb = (row_t*)b;
    return (ra->rho*ra->u  > rb->rho*rb->u) - (ra->rho*ra->u  < ra->rho*ra->u );
}

/* Subroutine to merge two input arrays into an output array. */
static void merge(void *out, const void *pa, size_t na,
                  const void *pb, size_t nb, size_t elem_size,
                  int (*cmp)(const void *, const void *))
{
    while (na != 0 || nb != 0) {
        if (na == 0 || (nb != 0 && cmp(pa, pb) > 0)) {
            memcpy(out, pb, elem_size);
            pb = (const char *)pb + elem_size;
            nb--;
        } else {
            memcpy(out, pa, elem_size);
            pa = (const char *)pa + elem_size;
            na--;
        }
        out = (char *)out + elem_size;
    }
}

/* Merge sort an array. Function modified from*/
void merge_sort(void *base, size_t n_memb, size_t elem_size,
                int (*cmp)(const void *, const void *))
{
    size_t n_bottom;
    size_t n_top;
    void *mid_p;
    void *bottom;
    void *top;

    if (n_memb <= 1) {
        /* Too small to sort. */
        return;
    }
    /* Sort the bottom half and the top half. */
    n_bottom = n_memb / 2;
    n_top = n_memb - n_bottom;
    mid_p = (char *)base + (n_bottom * elem_size);
    merge_sort(base, n_bottom, elem_size, cmp);
    merge_sort(mid_p, n_top, elem_size, cmp);
    /* Make temporary copies of the sorted bottom half and top half. */
    bottom = malloc(n_bottom * elem_size);
    top = malloc(n_top * elem_size);
    memcpy(bottom, base, n_bottom * elem_size);
    memcpy(top, mid_p, n_top * elem_size);
    /* Do a sorted merge of the copies into the original. */
    merge(base, bottom, n_bottom, top, n_top, elem_size, cmp);
    /* Free temporary copies. */
    free(bottom);
    free(top);
}

/*******************************************************************************
 * Task three: Search the data structures
*******************************************************************************/
void searching(const char* flow_file)
{
    FILE* input_file, *task3;
    row_t row = {0}, *center_points;
    list_t* row_list = list_new(free);
    bst_flux_t* bst_flux;
    char *headers = (char*)malloc(HEADER_CHARS * sizeof(char));
    int num_center = 0, n = 0;
    double* u_flux, closest;
    
    /* Reads in file */
    input_file = read_in_file(flow_file);
    fscanf(input_file,"%s", headers);

    /* Find the number of points on the center y */
    num_center = find_num_center(input_file, row);

    /* Reset to the start of the file */
    fseek(input_file, OFFSET, SEEK_SET);

    /* Scan in the input file */
    fscanf(input_file,"%s", headers);
    
    center_points = (row_t*)calloc(num_center,sizeof(row_t));
    u_flux = (double*)malloc(num_center*sizeof(double));

    /* Creates the array of center points */
    create_center_arr(input_file, row, center_points);

    /* Merge sorts them */
    merge_sort(center_points, num_center, sizeof(row_t), row_array_flux_comp);

    /* Creates linear flux array */
    create_flux_array(num_center, center_points, u_flux);

    /* create linked list */
    create_linked_list(num_center, center_points, row_list);

    /* create binary search tree */
    bst_flux = bst_flux_new(no_free);
    n = make_unique(u_flux, num_center);
    perfect_insert(bst_flux, u_flux, 0, n - 1);
    assert(bst_flux->num_elements == n);

    /* Open task 3 csv to write results to */
    task3 = fopen("task3.csv", "w");

    /* Find value closest to rho_40 */
    closest = find_closest(u_flux, num_center);

    /* Search the linear array and time it */
    gettimeofday(&start, NULL);

    linear_array_search(task3, num_center, u_flux, closest);

    fprintf(task3, "\n");
    gettimeofday(&stop, NULL);
    double elapsed_us = ((stop.tv_sec - start.tv_sec) * 1000000) + 
        (stop.tv_usec - start.tv_usec);
    printf("TASK 3 Array Linear Search:  %.2f microseconds\n", elapsed_us);
    
    /* Binary search the linear array and time it */
    gettimeofday(&start, NULL);
   
    binary_search(u_flux, 0, num_center-1, closest, task3);

    fprintf(task3, "\n");
    gettimeofday(&stop, NULL);
    elapsed_us = ((stop.tv_sec - start.tv_sec) * 1000000) + (stop.tv_usec - 
        start.tv_usec);
    printf("TASK 3 Array Binary Search:  %.2f microseconds\n", elapsed_us);

    /* Linear search the linked list and time it */
    gettimeofday(&start, NULL);

    linear_list_search(task3, row_list, closest);

    fprintf(task3, "\n");
    gettimeofday(&stop, NULL);
    elapsed_us = ((stop.tv_sec - start.tv_sec) * 1000000) + (stop.tv_usec - 
        start.tv_usec);
    printf("TASK 3 List Linear Search:  %.2f microseconds\n", elapsed_us);

    /* Search the Binary search the bst and time it */
    gettimeofday(&start, NULL);

    bst_flux_find(bst_flux, closest, task3);

    gettimeofday(&stop, NULL);
    elapsed_us = ((stop.tv_sec - start.tv_sec) * 1000000) + (stop.tv_usec - 
        start.tv_usec);
    printf("TASK 3 BST Search:  %.2f microseconds\n", elapsed_us);
    fprintf(task3, "\n");

    /* Clean up */
    bst_flux_free(bst_flux);
    free(u_flux);
    list_free(row_list);
    free(center_points);
    free(headers);
    fclose(input_file);
    fclose(task3);
}

/* Searches the inked list */
void linear_list_search(FILE* task3, list_t* row_list, double closest)
{
    l_node_t* cur;
    
    cur = row_list->head;
    row_t* list_row;
    double temp;

    /* while current node is a valid node */
    while(cur)
    {
        list_row = (row_t*)(cur->data);

        temp = (list_row->rho)*(list_row->u);

        if(temp == closest)
        {
            fprintf(task3, "%lf", temp);
            break;
        }

        fprintf(task3, "%lf,", temp);
        /* Go to next node */
        cur = cur->next;
    }
}

/* Linear search the array*/
void linear_array_search(FILE* task3, int num_center, double* u_flux, 
    double closest)
{ 
    for(int i=0;i<num_center;i++)
    {
        if(u_flux[i]==closest)
        {
            fprintf(task3, "%lf", u_flux[i]);
            break;
        } 
        fprintf(task3, "%lf,", u_flux[i]);
    }
}

/* Find the value closest to 40 percent of max rho */
double find_closest(double* u_flux, int num_center)
{
    double rho_40 = u_flux[num_center-1]*0.4, min_flux_diff = rho_40, 
        prev_diff, next_diff, flux_diff, closest;

    /* Loop through and find the value closest to rho_40 */
    for(int i = 0; i < num_center; i++)
    {
        
        flux_diff = fabs(rho_40 - u_flux[i]);
        if(i>0)
        {
            prev_diff = fabs(rho_40 - u_flux[i-1]);
        } 
        
        if(i<num_center-1)
        {
            next_diff = fabs(rho_40 - u_flux[i+1]);
        }

        if(flux_diff < min_flux_diff)
        {
            min_flux_diff = flux_diff;
            if(next_diff > min_flux_diff && prev_diff > min_flux_diff)
            {
                closest = u_flux[i];
                break;
            }

        }
    
    }
    return closest;
}

/* Create the linear flux array */
void create_flux_array(int num_center, row_t* center_points, double* u_flux)
{
    int i;
    for(i=0;i<num_center;i++)
    {
        u_flux[i] = (center_points[i].rho)*(center_points[i].u);
    }
}

/* Creates the linked list */
void create_linked_list(int num_center, row_t* center_points, list_t* row_list)
{
    int i;
    for(i=0;i<num_center;i++)
    {
        row_t* row_node = (row_t*)malloc(sizeof(row_t));

        if(row_node == NULL)
        {
            printf("Memory not allocated");
            exit(EXIT_SUCCESS);
        }

        row_node -> rho = center_points[i].rho;
        row_node -> u = center_points[i].u;
        row_node -> v = center_points[i].v;
        row_node -> x = center_points[i].x;
        row_node -> y = center_points[i].y;
        list_push_back(row_list, row_node);

    }
}

/* Creates the array of center points */
void create_center_arr(FILE* input_file, row_t row, row_t *center_points)
{
    int nums_scanned = 0, i = 0;

    while (nums_scanned != EOF)
    {
        nums_scanned = row_scanner(input_file, &row);
        if(row.y == 0)
        {
            center_points[i] = row;
            i++;
        }
    }
}

/* Count the number of center points */
int find_num_center(FILE* input_file, row_t row)
{
    int nums_scanned = 0, num_center = 0;
    while (nums_scanned != EOF)
    {     
        nums_scanned = row_scanner(input_file, &row);
        if(row.y == 0)
        {
            num_center++;
        }
    }
    return num_center;
}

int binary_search(double arr[], int l, int r, double x, FILE* task3) 
{   

    if (r >= l) { 
        int mid = l + (r - l) / 2;
        /* If the element is present at the middle itself */
        if (arr[mid] == x)
        {
            fprintf(task3, "%lf", arr[mid]);
            
            return mid; 
        }

        fprintf(task3, "%lf,", arr[mid]);
        
        /* If element is smaller than mid, then it can only be present in left 
        subarray */
        if (arr[mid] > x) 
            return binary_search(arr, l, mid - 1, x, task3); 
  
        /* Else the element can only be present in right subarray */
        return binary_search(arr, mid + 1, r, x, task3); 
    } 
  
    /* We reach here when element is not present in array */
    return -1; 
}

/*******************************************************************************
 *  This code below is in part provided by the University of Melbourne,
 *  ENGR30003 by mpetri, see repo here.
 *  https://github.com/mpetri/ENGR3K3
*******************************************************************************/
bst_flux_node_t* bst_flux_find(bst_flux_t* bst_flux, double d, FILE* task3) {
	assert(bst_flux != NULL);
	
	bst_flux_node_t* tmp = bst_flux->root;
    double flux_value, *flux_value_ptr = &flux_value;

	while(tmp) {
        flux_value_ptr = (double*)(tmp->data);
  

		if(*flux_value_ptr - d > 0) 
        { // element is smaller
			tmp = tmp->left;
		} else if( *flux_value_ptr - d < 0) { // element is bigger
			tmp = tmp->right;
		} else if(*flux_value_ptr - d == 0){
			fprintf(task3, "%lf", *flux_value_ptr);
         
			break;
		}
        fprintf(task3, "%lf,", *flux_value_ptr);
        
        
	}
	return tmp;
}

/* Create a new empty bst_flux structure */
bst_flux_t* bst_flux_new(void (*delfunc)(void*))
{
    bst_flux_t* bst_flux;
    bst_flux = (bst_flux_t*)malloc(sizeof(bst_flux_t));
    assert(bst_flux != NULL);
    bst_flux->root = NULL;
    bst_flux->num_elements = 0;
    bst_flux->del = delfunc;
    
    return bst_flux;
}

/* Free all memory assocated with a subtree */
void bst_flux_free_subtree(bst_flux_t* bst_flux, bst_flux_node_t* n)
{
    assert(bst_flux != NULL);
    if (n) {
        bst_flux_free_subtree(bst_flux, n->left);
        bst_flux_free_subtree(bst_flux, n->right);
        bst_flux->del(n->data);
        free(n);
        bst_flux->num_elements--;
    }
}

/* free all memory associated with a bst_flux */
void bst_flux_free(bst_flux_t* bst_flux)
{
    assert(bst_flux != NULL);
    bst_flux_free_subtree(bst_flux, bst_flux->root);
    free(bst_flux);
}

/* insert a new element into the bst_flux */

int bst_flux_insert(bst_flux_t* bst_flux, void* d)
{
    assert(bst_flux != NULL);
    assert(d != NULL);
    bst_flux_node_t* parent = NULL;
    bst_flux_node_t* tmp = bst_flux->root;
    while (tmp) {
        parent = tmp;
        if ((tmp->data - d) > 0) { // element is smaller
            tmp = tmp->left;
        }
        else if ((tmp->data - d) < 0) { // element is bigger
            tmp = tmp->right;
        }
        else {
            /* ALREADY EXISTS! -> do nothing and return ERROR */
            return BST_FAILURE;
        }
    }

    /* insert as child of parent */
    bst_flux_node_t* new_node = 
        (bst_flux_node_t*)malloc(sizeof(bst_flux_node_t));
    assert(new_node != NULL);
    new_node->data = d;
    new_node->left = NULL;
    new_node->right = NULL;
    if (parent != NULL) {
        if ((parent->data - d) > 0) { // element is smaller
            assert(parent->left == NULL);
            parent->left = new_node;
        }
        else {
            assert(parent->right == NULL);
            parent->right = new_node;
        }
    }
    else {
        assert(bst_flux->root == NULL);
        bst_flux->root = new_node;
    }
    bst_flux->num_elements++;

    return BST_SUCCESS;
}

void perfect_insert(bst_flux_t* bst_flux, double* array, int low, int high)
{
    if (low <= high) {
    	// Choose root from array and insert
    	// Recursively do the same on left and right (1)
        int mid = low + (high - low) / 2;
        double* ptr = array + mid;
        bst_flux_insert(bst_flux, ptr);
        perfect_insert(bst_flux, array, low, mid - 1);
        perfect_insert(bst_flux, array, mid + 1, high);
    }
}

void no_free(void* v)
{
}

int make_unique(double* array, int n)
{
    int dest = 0;
    int itr = 1;
    while (itr != n) {
        if (array[dest] != array[itr]) {
            array[++dest] = array[itr];
        }
        itr++;
    }
    return dest+1;
}

list_t* list_new(void (*delfunc)(void*))
{
    list_t* list;
    list = (list_t*)malloc(sizeof(list_t));
    assert(list != NULL);
    list->head = NULL;
    list->tail = NULL;
    list->num_elements = 0;
    list->del = delfunc;
    return list;
}

/* process all elements in the list */
void list_process_all(l_node_t* p, void (*process)(void*))
{
    while (p) {
        process(p->data);
        p = p->next;
    }
}

/* remove list_node at the front of the list */
void* list_pop_front(list_t* list)
{
    assert(list != NULL);
    assert(list->num_elements > 0);
    l_node_t* old;
    assert(list->head != NULL);
    old = list->head;
    list->head = list->head->next;
    void* d = old->data;
    free(old);
    list->num_elements--;
    if (list->num_elements == 0) {
        list->head = NULL;
        list->tail = NULL;
    }
    return d;
}

/* add list_node add the back of the list */
void list_push_back(list_t* list, void* d)
{
    assert(list != NULL);
    l_node_t* new = (l_node_t*)malloc(sizeof(l_node_t));
    assert(new != NULL);
    new->data = d;
    new->next = NULL;
    if (list->tail)
        list->tail->next = new;
    list->tail = new;
    if (list->head == NULL)
        list->head = new;
    list->num_elements++;
}

/* free all memory associated with a list */
void list_free(list_t* list)
{
    assert(list != NULL);
    while (list->num_elements) {
        void* elem = list_pop_front(list);
        list->del(elem);
    }
    free(list);
}

/*******************************************************************************
 *  This code above is in part provided by the University of Melbourne,
 *  ENGR30003 by mpetri, see repo here.
 *  https://github.com/mpetri/ENGR3K3
*******************************************************************************/

/*******************************************************************************
 * Task four: Computing the vorticity
*******************************************************************************/
void vortcalc(const char* flow_file)
{
    FILE* input_file;
    int i, j;
    long unsigned int num_thres_1 = 0, num_thres_2 = 0, num_thres_3 = 0, 
        num_thres_4 = 0, num_thres_5 = 0;
    grid_t grid;
    row_t row = {0};
    char *headers = (char*)malloc(HEADER_CHARS * sizeof(char));

    input_file = read_in_file(flow_file);
    fscanf(input_file,"%s", headers);

    /* Calculates n and m */
    calc_n_m(input_file, &grid);

    fseek(input_file, OFFSET, SEEK_SET);

    fscanf(input_file,"%s", headers);

    /* Assings memory to the 2D array */
    double (*vort)[grid.m + 1] = malloc(sizeof(double[grid.n + 1][grid.m + 1]));
    double (*U)[grid.m + 1] = malloc(sizeof(double[grid.n + 1][grid.m + 1]));
    double (*V)[grid.m + 1] = malloc(sizeof(double[grid.n + 1][grid.m + 1]));
    double (*X)[grid.m + 1] = malloc(sizeof(double[grid.n + 1][grid.m + 1]));
    double (*Y)[grid.m + 1] = malloc(sizeof(double[grid.n + 1][grid.m + 1]));

    /* Setting the values */
    for(j = 0; j < grid.m; j++)
    {
        for(i = 0; i<grid.n; i++)
        {
            row_scanner(input_file, &row);
            
            U[i][j] = row.u;
            V[i][j] = row.v;
            X[i][j] = row.x;
            Y[i][j] = row.y;
        }
    }
    
    /* Computing the vort for the whole grid */
    for(j = 0; j < grid.m; j++)
    {
        for(i = 0; i<grid.n; i++)
        {            
            if(i == grid.n - 1 && j != grid.m -1) // Edge case 1
            {  
                vort[i][j] = (V[i][j] - V[i-1][j])/(X[i][j] - X[i-1][j]) - 
                    (U[i][j] - U[i][j])/(Y[i][j+1] - Y[i][j]);
            }
            else if (j == grid.m -1 && i != grid.n - 1 ) // Edge case 2
            {
                
                vort[i][j] = (V[i+1][j] - V[i][j])/(X[i+1][j] - X[i][j]) - 
                    (U[i][j] - U[i][j-1])/(Y[i][j] - Y[i][j - 1]);
            }
            else if(j == grid.m - 1 && i == grid.n - 1) // Edge case 3
            {
                vort[i][j] = (V[i][j] - V[i-1][j])/(X[i][j] - X[i-1][j]) - 
                    (U[i][j] - U[i][j-1])/(Y[i][j] - Y[i][j - 1]);
            }
            else // Normal case
            {
                vort[i][j] = (V[i+1][j] - V[i][j])/(X[i+1][j] - X[i][j]) - 
                    (U[i][j+1] - U[i][j])/(Y[i][j+1] - Y[i][j]);
            }
        }
    }

    /* Assign threshold counts */
    for(j = 0; j < grid.m; j++)
    {
        for(i = 0; i<grid.n; i++)
        {
            if(fabs(vort[i][j]) >= THRESHOLD_0 && 
                fabs(vort[i][j]) < THRESHOLD_1)
            {
                num_thres_1++;
            } 
            else if(fabs(vort[i][j]) >= THRESHOLD_1 && 
                fabs(vort[i][j]) < THRESHOLD_2)
            {
                num_thres_2++;
            }
            else if(fabs(vort[i][j]) >= THRESHOLD_2 && 
                fabs(vort[i][j]) < THRESHOLD_3)
            {
                num_thres_3++;
            }
            else if(fabs(vort[i][j]) >= THRESHOLD_3 && 
                fabs(vort[i][j]) < THRESHOLD_4)
            {
                num_thres_4++;
            }
            else if(fabs(vort[i][j]) >= THRESHOLD_4 && 
                fabs(vort[i][j]) < THRESHOLD_5)
            {
                num_thres_5++;
             
            }
        }
    }

    /* Print to task 4 */
    FILE* task_4_file = fopen("task4.csv", "w");
    
    fprintf(task_4_file, "threshold ,points\n");
    fprintf(task_4_file, "%d ,%lu\n", THRESHOLD_1, num_thres_1);
    fprintf(task_4_file, "%d ,%lu\n", THRESHOLD_2, num_thres_2);
    fprintf(task_4_file, "%d ,%lu\n", THRESHOLD_3, num_thres_3);
    fprintf(task_4_file, "%d ,%lu\n", THRESHOLD_4, num_thres_4);
    fprintf(task_4_file, "%d ,%lu\n", THRESHOLD_5, num_thres_5);

    /* Clean up */
    fclose(input_file);
    fclose(task_4_file);
    free(vort);
    free(U);
    free(V);
    free(X);
    free(Y);
    free(headers);
}

void calc_n_m(FILE* input_file, grid_t *grid)
{
    /* Finds the grid dimentions of the data */
    row_t row = {0}, prev_row = {0};
    int num_scanned = 0;
    int n = 0;
    long unsigned int total = 0;

    num_scanned = row_scanner(input_file, &row);
    prev_row = row;

    /* Finds x dimention */
    while(row.y == prev_row.y && num_scanned != EOF)
    {
        prev_row = row;
        num_scanned = row_scanner(input_file, &row);
        n++;
    }

    grid->n = n;
    total = grid -> n;
    /* Finds total */
    while(num_scanned != EOF)
    {
        num_scanned = row_scanner(input_file, &row);
        total++;
    }
    
    /* Calculates m */
    grid->total_lines = total;
    grid->m = total/grid->n;
    assert(grid->total_lines = grid->m*grid->n);
}