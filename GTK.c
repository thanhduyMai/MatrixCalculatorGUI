#include <gtk/gtk.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

typedef struct {
    int rows;
    int cols;
    double **data;
} Matrix;

typedef struct {
    GtkWidget *parent_window;
    Matrix matrix;
    Matrix new_matrix;
    Matrix result_matrix;
    int matrix_loaded;
    int new_matrix_loaded;
    int result_matrix_loaded;
    Matrix add_result;        // Add result matrix
    Matrix subtract_result;   // Subtract result matrix
    Matrix multiply_result;   // Multiply result matrix
    Matrix scalar_result;     // Scalar multiply result matrix
     // Results storage
    double determinant_result;
    Matrix inverse_result;
    int rank_result;
    double *eigenvalues_result;
    Matrix eigenvectors_result;
    int add_result_loaded;        // Flag for add result
    int subtract_result_loaded;   // Flag for subtract result
    int multiply_result_loaded;   // Flag for multiply result
    int scalar_result_loaded;     // Flag for scalar multiply result
    int results_calculated;
    int eigenvectors_result_loaded;
} AppData;

// Function declarations
Matrix add_matrices(Matrix A, Matrix B);
Matrix subtract_matrices(Matrix A, Matrix B);
Matrix multiply_matrices(Matrix A, Matrix B);
Matrix scalar_multiply_matrix(Matrix A, double scalar);
double determinant(Matrix A);
Matrix inverse_matrix(Matrix A);
int matrix_rank(Matrix A);
void show_matrix_dialog(GtkWidget *parent_window, Matrix matrix);
void on_load_matrix(GtkWidget *widget, gpointer data);
void on_load_new_matrix(GtkWidget *widget, gpointer data);
void on_generate_random_matrix(GtkWidget *widget, gpointer data);
void on_add_matrices(GtkWidget *widget, gpointer data);
void on_subtract_matrices(GtkWidget *widget, gpointer data);
void on_multiply_matrices(GtkWidget *widget, gpointer data);
void on_scalar_multiply_matrix(GtkWidget *widget, gpointer data);
void on_determinant(GtkWidget *widget, gpointer data);
void on_inverse(GtkWidget *widget, gpointer data);
void on_rank(GtkWidget *widget, gpointer data);
void eigenvalues_and_eigenvectors(Matrix A, double *eigenvalues, Matrix *eigenvectors);
void on_eigenvalues_and_eigenvectors(GtkWidget *widget, gpointer data);
void on_save_matrix(GtkWidget *widget, gpointer data);
void free_matrix(Matrix matrix);

int main(int argc, char *argv[]) {
    gtk_init(&argc, &argv);
    GtkWidget *window;
    GtkWidget *vbox;
    GtkWidget *hbox;
    GtkWidget *load_button;
    GtkWidget *load_new_button;
    GtkWidget *save_button;
    GtkWidget *generate_button;
    GtkWidget *add_button;
    GtkWidget *subtract_button;
    GtkWidget *multiply_button;
    GtkWidget *scalar_multiply_button;
    GtkWidget *determinant_button;
    GtkWidget *inverse_button;
    GtkWidget *rank_button;
    GtkWidget *eigen_button;
    
    AppData app_data;
    app_data.matrix_loaded = FALSE;
    app_data.new_matrix_loaded = FALSE;
    app_data.result_matrix_loaded = FALSE;
    app_data.results_calculated = FALSE;
    app_data.eigenvalues_result = NULL;

    window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), "Matrix Calculator-Mai Thành Duy-20227225");
    gtk_container_set_border_width(GTK_CONTAINER(window), 10);
    gtk_window_set_default_size(GTK_WINDOW(window), 400, 300);
    g_signal_connect(window, "destroy", G_CALLBACK(gtk_main_quit), NULL);
    
    app_data.parent_window = window;

    vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 20);
    gtk_container_add(GTK_CONTAINER(window), vbox);
    
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    load_button = gtk_button_new_with_label("Load Matrix from File");
    g_signal_connect(load_button, "clicked", G_CALLBACK(on_load_matrix), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), load_button, TRUE, TRUE, 0);

    load_new_button = gtk_button_new_with_label("Load Another Matrix");
    g_signal_connect(load_new_button, "clicked", G_CALLBACK(on_load_new_matrix), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), load_new_button, TRUE, TRUE, 0);

    generate_button = gtk_button_new_with_label("Generate Random Matrix");
    g_signal_connect(generate_button, "clicked", G_CALLBACK(on_generate_random_matrix), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), generate_button, TRUE, TRUE, 0);

    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    add_button = gtk_button_new_with_label("Add Matrices");
    g_signal_connect(add_button, "clicked", G_CALLBACK(on_add_matrices), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), add_button, TRUE, TRUE, 0);
    
    subtract_button = gtk_button_new_with_label("Subtract Matrices");
    g_signal_connect(subtract_button, "clicked", G_CALLBACK(on_subtract_matrices), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), subtract_button, TRUE, TRUE, 0);

    multiply_button = gtk_button_new_with_label("Multiply Matrices");
    g_signal_connect(multiply_button, "clicked", G_CALLBACK(on_multiply_matrices), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), multiply_button, TRUE, TRUE, 0);

    scalar_multiply_button = gtk_button_new_with_label("Multiply by Scalar");
    g_signal_connect(scalar_multiply_button, "clicked", G_CALLBACK(on_scalar_multiply_matrix), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), scalar_multiply_button, TRUE, TRUE, 0);

    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    determinant_button = gtk_button_new_with_label("Determinant");
    g_signal_connect(determinant_button, "clicked", G_CALLBACK(on_determinant), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), determinant_button, TRUE, TRUE, 0);

    inverse_button = gtk_button_new_with_label("Inverse");
    g_signal_connect(inverse_button, "clicked", G_CALLBACK(on_inverse), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), inverse_button, TRUE, TRUE, 0);

    rank_button = gtk_button_new_with_label("Rank");
    g_signal_connect(rank_button, "clicked", G_CALLBACK(on_rank), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), rank_button, TRUE, TRUE, 0);

    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    eigen_button = gtk_button_new_with_label("Eigenvalues and Eigenvectors");
    g_signal_connect(eigen_button, "clicked", G_CALLBACK(on_eigenvalues_and_eigenvectors), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), eigen_button, TRUE, TRUE, 0);
    
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    save_button = gtk_button_new_with_label("Save Results to File");
    g_signal_connect(save_button, "clicked", G_CALLBACK(on_save_matrix), &app_data);
    gtk_box_pack_start(GTK_BOX(hbox), save_button, TRUE, TRUE, 0);

    gtk_widget_show_all(window);
    gtk_main();

    if (app_data.matrix_loaded) {
        free_matrix(app_data.matrix);
    }

    if (app_data.new_matrix_loaded) {
        free_matrix(app_data.new_matrix);
    }

    return 0;
}
// Function of Calculating Matrix

Matrix add_matrices(Matrix A, Matrix B) {
    if (A.rows != B.rows || A.cols != B.cols) {
        printf("Matrices must have the same dimensions for addition\n");
        exit(1);
    }

    Matrix result;
    result.rows = A.rows;
    result.cols = A.cols;
    result.data = (double **)malloc(result.rows * sizeof(double *));
    for (int i = 0; i < result.rows; ++i) {
        result.data[i] = (double *)malloc(result.cols * sizeof(double));
        for (int j = 0; j < result.cols; ++j) {
            result.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }

    return result;
}

Matrix subtract_matrices(Matrix A, Matrix B) {
    if (A.rows != B.rows || A.cols != B.cols) {
        printf("Matrices must have the same dimensions for subtraction\n");
        exit(1);
    }

    Matrix result;
    result.rows = A.rows;
    result.cols = A.cols;
    result.data = (double **)malloc(result.rows * sizeof(double *));
    for (int i = 0; i < result.rows; ++i) {
        result.data[i] = (double *)malloc(result.cols * sizeof(double));
        for (int j = 0; j < result.cols; ++j) {
            result.data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }

    return result;
}

Matrix multiply_matrices(Matrix A, Matrix B) {
    if (A.cols != B.rows) {
        printf("Number of columns of the first matrix must equal number of rows of the second matrix for multiplication\n");
        exit(1);
    }

    Matrix result;
    result.rows = A.rows;
    result.cols = B.cols;
    result.data = (double **)malloc(result.rows * sizeof(double *));
    for (int i = 0; i < result.rows; ++i) {
        result.data[i] = (double *)malloc(result.cols * sizeof(double));
        for (int j = 0; j < result.cols; ++j) {
            result.data[i][j] = 0;
            for (int k = 0; k < A.cols; ++k) {
                result.data[i][j] += A.data[i][k] * B.data[k][j];
            }
        }
    }

    return result;
}

Matrix scalar_multiply_matrix(Matrix A, double scalar) {
    Matrix result;
    result.rows = A.rows;
    result.cols = A.cols;
    result.data = (double **)malloc(result.rows * sizeof(double *));
    for (int i = 0; i < result.rows; ++i) {
        result.data[i] = (double *)malloc(result.cols * sizeof(double));
        for (int j = 0; j < result.cols; ++j) {
            result.data[i][j] = A.data[i][j] * scalar;
        }
    }

    return result;
}

double determinant(Matrix A) {
    if (A.rows != A.cols) {
        printf("Matrix must be square to calculate determinant\n");
        exit(1);
    }

    int n = A.rows;
    double det = 1;
    Matrix B;
    B.rows = n;
    B.cols = n;
    B.data = (double **)malloc(n * sizeof(double *));
    if (B.data == NULL) {
        perror("Failed to allocate memory");
        exit(1);
    }

    for (int i = 0; i < n; ++i) {
        B.data[i] = (double *)malloc(n * sizeof(double));
        if (B.data[i] == NULL) {
            perror("Failed to allocate memory");
            exit(1);
        }
        for (int j = 0; j < n; ++j) {
            B.data[i][j] = A.data[i][j];
        }
    }

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(B.data[j][i]) > fabs(B.data[pivot][i])) {
                pivot = j;
            }
        }

        if (i != pivot) {
            double *temp = B.data[i];
            B.data[i] = B.data[pivot];
            B.data[pivot] = temp;
            det *= -1;
        }

        if (B.data[i][i] == 0) {
            for (int i = 0; i < n; ++i) {
                free(B.data[i]);
            }
            free(B.data);
            return 0;
        }

        for (int j = i + 1; j < n; ++j) {
            double factor = B.data[j][i] / B.data[i][i];
            for (int k = i; k < n; ++k) {
                B.data[j][k] -= factor * B.data[i][k];
            }
        }

        det *= B.data[i][i];
    }

    for (int i = 0; i < n; ++i) {
        free(B.data[i]);
    }
    free(B.data);

    return det;
}


Matrix inverse_matrix(Matrix A) {
    if (A.rows != A.cols) {
        printf("Matrix is not square\n");
        exit(1);
    }
    int n = A.rows;
    Matrix inverse;
    inverse.rows = n;
    inverse.cols = n;
    inverse.data = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; ++i) {
        inverse.data[i] = (double *)malloc(n * sizeof(double));
        for (int j = 0; j < n; ++j) {
            inverse.data[i][j] = (i == j) ? 1 : 0;
        }
    }

    double **matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; ++i) {
        matrix[i] = (double *)malloc(n * sizeof(double));
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = A.data[i][j];
        }
    }

    for (int i = 0; i < n; ++i) {
        int maxIdx = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(matrix[j][i]) > fabs(matrix[maxIdx][i])) {
                maxIdx = j;
            }
        }

        if (i != maxIdx) {
            double *temp = matrix[i];
            matrix[i] = matrix[maxIdx];
            matrix[maxIdx] = temp;

            temp = inverse.data[i];
            inverse.data[i] = inverse.data[maxIdx];
            inverse.data[maxIdx] = temp;
        }

        if (matrix[i][i] == 0) {
            printf("Matrix is singular and cannot be inverted\n");
            exit(1);
        }

        double pivot = matrix[i][i];
        for (int j = 0; j < n; ++j) {
            matrix[i][j] /= pivot;
            inverse.data[i][j] /= pivot;
        }

        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double factor = matrix[j][i];
                for (int k = 0; k < n; ++k) {
                    matrix[j][k] -= factor * matrix[i][k];
                    inverse.data[j][k] -= factor * inverse.data[i][k];
                }
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        free(matrix[i]);
    }
    free(matrix);

    return inverse;
}

int matrix_rank(Matrix A) {
    int m = A.rows;
    int n = A.cols;
    int rank = n;

    double **matrix = (double **)malloc(m * sizeof(double *));
    for (int i = 0; i < m; ++i) {
        matrix[i] = (double *)malloc(n * sizeof(double));
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = A.data[i][j];
        }
    }

    for (int row = 0; row < rank; row++) {
        if (matrix[row][row]) {
            for (int col = 0; col < m; col++) {
                if (col != row) {
                    double mult = matrix[col][row] / matrix[row][row];
                    for (int i = 0; i < rank; i++) {
                        matrix[col][i] -= mult * matrix[row][i];
                    }
                }
            }
        } else {
            int reduce = 1;
            for (int i = row + 1; i < m; i++) {
                if (matrix[i][row]) {
                    double *temp = matrix[row];
                    matrix[row] = matrix[i];
                    matrix[i] = temp;
                    reduce = 0;
                    break;
                }
            }

            if (reduce) {
                rank--;
                for (int i = 0; i < m; i++) {
                    matrix[i][row] = matrix[i][rank];
                }
            }

            row--;
        }
    }

    for (int i = 0; i < m; ++i) {
        free(matrix[i]);
    }
    free(matrix);

    return rank;
}

void eigenvalues_and_eigenvectors(Matrix A, double *eigenvalues, Matrix *eigenvectors) {
    if (A.rows != A.cols) {
        printf("Matrix must be square to calculate eigenvalues and eigenvectors\n");
        exit(1);
    }

    int n = A.rows;
    int max_iterations = 1000;
    double tolerance = 1e-10;

    eigenvectors->rows = n;
    eigenvectors->cols = n;
    eigenvectors->data = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; ++i) {
        eigenvectors->data[i] = (double *)malloc(n * sizeof(double));
    }

    for (int k = 0; k < n; ++k) {
        double *b = (double *)malloc(n * sizeof(double));
        for (int i = 0; i < n; ++i) {
            b[i] = (double)rand() / RAND_MAX;
        }

        double eigenvalue = 0;
        for (int iter = 0; iter < max_iterations; ++iter) {
            double *Ab = (double *)malloc(n * sizeof(double));
            for (int i = 0; i < n; ++i) {
                Ab[i] = 0;
                for (int j = 0; j < n; ++j) {
                    Ab[i] += A.data[i][j] * b[j];
                }
            }

            double norm = 0;
            for (int i = 0; i < n; ++i) {
                norm += Ab[i] * Ab[i];
            }
            norm = sqrt(norm);

            for (int i = 0; i < n; ++i) {
                b[i] = Ab[i] / norm;
            }

            double new_eigenvalue = 0;
            for (int i = 0; i < n; ++i) {
                new_eigenvalue += b[i] * Ab[i];
            }

            free(Ab);

            if (fabs(new_eigenvalue - eigenvalue) < tolerance) {
                break;
            }

            eigenvalue = new_eigenvalue;
        }

        eigenvalues[k] = eigenvalue;
        for (int i = 0; i < n; ++i) {
            eigenvectors->data[i][k] = b[i];
        }

        free(b);
    }
}
// Function of handle events
void show_matrix_dialog(GtkWidget *parent_window, Matrix matrix) {
    GtkWidget *dialog = gtk_dialog_new_with_buttons("Matrix", GTK_WINDOW(parent_window),
        GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT, "_OK", GTK_RESPONSE_OK, NULL);
    GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    GtkWidget *table = gtk_grid_new();

    for (int i = 0; i < matrix.rows; ++i) {
        for (int j = 0; j < matrix.cols; ++j) {
            char buffer[32];
            snprintf(buffer, sizeof(buffer), "%.2f", matrix.data[i][j]);
            GtkWidget *label = gtk_label_new(buffer);
            gtk_grid_attach(GTK_GRID(table), label, j, i, 1, 1);
        }
    }

    gtk_container_add(GTK_CONTAINER(content_area), table);
    gtk_widget_show_all(dialog);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
}

void on_load_matrix(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File", GTK_WINDOW(app_data->parent_window),
        GTK_FILE_CHOOSER_ACTION_OPEN, "_Cancel", GTK_RESPONSE_CANCEL, "_Open", GTK_RESPONSE_ACCEPT, NULL);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        FILE *file = fopen(filename, "r");
        if (file) {
            fscanf(file, "%d %d", &app_data->matrix.rows, &app_data->matrix.cols);
            app_data->matrix.data = (double **)malloc(app_data->matrix.rows * sizeof(double *));
            for (int i = 0; i < app_data->matrix.rows; ++i) {
                app_data->matrix.data[i] = (double *)malloc(app_data->matrix.cols * sizeof(double));
                for (int j = 0; j < app_data->matrix.cols; ++j) {
                    fscanf(file, "%lf", &app_data->matrix.data[i][j]);
                }
            }
            fclose(file);
            app_data->matrix_loaded = TRUE;
        }
        g_free(filename);
    }

    gtk_widget_destroy(dialog);
}

void on_load_new_matrix(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File", GTK_WINDOW(app_data->parent_window),
        GTK_FILE_CHOOSER_ACTION_OPEN, "_Cancel", GTK_RESPONSE_CANCEL, "_Open", GTK_RESPONSE_ACCEPT, NULL);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        FILE *file = fopen(filename, "r");
        if (file) {
            fscanf(file, "%d %d", &app_data->new_matrix.rows, &app_data->new_matrix.cols);
            app_data->new_matrix.data = (double **)malloc(app_data->new_matrix.rows * sizeof(double *));
            for (int i = 0; i < app_data->new_matrix.rows; ++i) {
                app_data->new_matrix.data[i] = (double *)malloc(app_data->new_matrix.cols * sizeof(double));
                for (int j = 0; j < app_data->new_matrix.cols; ++j) {
                    fscanf(file, "%lf", &app_data->new_matrix.data[i][j]);
                }
            }
            fclose(file);
            app_data->new_matrix_loaded = TRUE;
        }
        g_free(filename);
    }

    gtk_widget_destroy(dialog);
}

void on_generate_random_matrix(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    GtkWidget *dialog = gtk_dialog_new_with_buttons("Generate Random Matrix", GTK_WINDOW(app_data->parent_window),
        GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT, "_OK", GTK_RESPONSE_OK, "_Cancel", GTK_RESPONSE_CANCEL, NULL);
    GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    GtkWidget *grid = gtk_grid_new();

    GtkWidget *rows_label = gtk_label_new("Rows:");
    GtkWidget *rows_entry = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(rows_entry), "3");
    GtkWidget *cols_label = gtk_label_new("Columns:");
    GtkWidget *cols_entry = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(cols_entry), "3");

    gtk_grid_attach(GTK_GRID(grid), rows_label, 0, 0, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), rows_entry, 1, 0, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), cols_label, 0, 1, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), cols_entry, 1, 1, 1, 1);

    gtk_container_add(GTK_CONTAINER(content_area), grid);
    gtk_widget_show_all(dialog);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK) {
        int rows = atoi(gtk_entry_get_text(GTK_ENTRY(rows_entry)));
        int cols = atoi(gtk_entry_get_text(GTK_ENTRY(cols_entry)));

        app_data->matrix.rows = rows;
        app_data->matrix.cols = cols;
        app_data->matrix.data = (double **)malloc(rows * sizeof(double *));
        srand(time(NULL));
        for (int i = 0; i < rows; ++i) {
            app_data->matrix.data[i] = (double *)malloc(cols * sizeof(double));
            for (int j = 0; j < cols; ++j) {
                app_data->matrix.data[i][j] = (double)(rand() % 100) / 10.0;
            }
        }

        app_data->matrix_loaded = TRUE;
        show_matrix_dialog(app_data->parent_window, app_data->matrix);
    }

    gtk_widget_destroy(dialog);
}

void on_add_matrices(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded || !app_data->new_matrix_loaded) {
        printf("Both matrices must be loaded to perform addition\n");
        return;
    }

    Matrix result = add_matrices(app_data->matrix, app_data->new_matrix);
    show_matrix_dialog(app_data->parent_window, result);

    if (app_data->matrix_loaded && app_data->new_matrix_loaded) {
        Matrix result = add_matrices(app_data->matrix, app_data->new_matrix);
         app_data->add_result = add_matrices(app_data->matrix, app_data->new_matrix);
         app_data->add_result_loaded = TRUE;
    }
    
}

void on_subtract_matrices(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded || !app_data->new_matrix_loaded) {
        printf("Both matrices must be loaded to perform subtraction\n");
        return;
    }

    Matrix result = subtract_matrices(app_data->matrix, app_data->new_matrix);
    show_matrix_dialog(app_data->parent_window, result);
    
    if (app_data->matrix_loaded && app_data->new_matrix_loaded) {
        Matrix result = subtract_matrices(app_data->matrix, app_data->new_matrix);
        app_data->subtract_result = subtract_matrices(app_data->matrix, app_data->new_matrix);
        app_data->subtract_result_loaded = TRUE;
    }
}

void on_multiply_matrices(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded || !app_data->new_matrix_loaded) {
        printf("Both matrices must be loaded to perform multiplication\n");
        return;
    }

    Matrix result = multiply_matrices(app_data->matrix, app_data->new_matrix);
    show_matrix_dialog(app_data->parent_window, result);
    if (app_data->matrix_loaded && app_data->new_matrix_loaded) {
        Matrix result = multiply_matrices(app_data->matrix, app_data->new_matrix);
         app_data->multiply_result = multiply_matrices(app_data->matrix, app_data->new_matrix);
         app_data->multiply_result_loaded = TRUE;
    }
}

void on_scalar_multiply_matrix(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded) {
        printf("Matrix must be loaded to perform scalar multiplication\n");
        return;
    }

    GtkWidget *dialog = gtk_dialog_new_with_buttons("Multiply by Scalar", GTK_WINDOW(app_data->parent_window),
        GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT, "_OK", GTK_RESPONSE_OK, "_Cancel", GTK_RESPONSE_CANCEL, NULL);
    GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    GtkWidget *grid = gtk_grid_new();

    GtkWidget *scalar_label = gtk_label_new("Scalar:");
    GtkWidget *scalar_entry = gtk_entry_new();

    gtk_grid_attach(GTK_GRID(grid), scalar_label, 0, 0, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), scalar_entry, 1, 0, 1, 1);

    gtk_container_add(GTK_CONTAINER(content_area), grid);
    gtk_widget_show_all(dialog);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK) {
        double scalar = atof(gtk_entry_get_text(GTK_ENTRY(scalar_entry)));

        Matrix result = scalar_multiply_matrix(app_data->matrix, scalar);
        show_matrix_dialog(app_data->parent_window, result);
        
        // Store the result in app_data
        app_data->scalar_result = scalar_multiply_matrix(app_data->matrix, scalar);
        app_data->scalar_result_loaded = TRUE;
    }

    gtk_widget_destroy(dialog);

}

void on_determinant(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded) {
        printf("Matrix must be loaded to calculate determinant\n");
        return;
    }

    // Debugging information
    if (!GTK_IS_WINDOW(app_data->parent_window)) {
        printf("Invalid parent window\n");
        g_print("Parent window pointer: %p\n", app_data->parent_window);
        g_print("GTK_IS_WINDOW check: %d\n", GTK_IS_WINDOW(app_data->parent_window));
        return;
    }

    double det = determinant(app_data->matrix);
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "Determinant: %.2f", det);
    app_data->determinant_result = det;
    app_data->results_calculated = TRUE;

    GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(app_data->parent_window),
                                               GTK_DIALOG_MODAL,
                                               GTK_MESSAGE_INFO,
                                               GTK_BUTTONS_OK,
                                               "%s", buffer);

    if (dialog == NULL) {
        printf("Failed to create message dialog\n");
        return;
    }
    
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
}

void on_inverse(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded) {
        printf("Matrix must be loaded to calculate inverse\n");
        return;
    }

    Matrix result = inverse_matrix(app_data->matrix);
    show_matrix_dialog(app_data->parent_window, result);
    app_data->inverse_result = result;
    app_data->results_calculated = TRUE;
}

void on_rank(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded) {
        printf("Matrix must be loaded to calculate rank\n");
        return;
    }

    int rank = matrix_rank(app_data->matrix);
    app_data->rank_result = rank;
    app_data->results_calculated = TRUE;
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "Rank: %d", rank);
    GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(app_data->parent_window), GTK_DIALOG_MODAL, GTK_MESSAGE_INFO, GTK_BUTTONS_OK, "%s", buffer);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
}

void on_eigenvalues_and_eigenvectors(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded) {
        printf("Matrix must be loaded to calculate eigenvalues and eigenvectors\n");
        return;
    }

    int n = app_data->matrix.rows;
    double *eigenvalues = (double *)malloc(n * sizeof(double));
    Matrix eigenvectors;
    eigenvalues_and_eigenvectors(app_data->matrix, eigenvalues, &eigenvectors);
    app_data->eigenvalues_result = eigenvalues;
    app_data->eigenvectors_result.rows = n;
    app_data->eigenvectors_result.cols = n;
    app_data->eigenvectors_result.data = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; ++i) {
        app_data->eigenvectors_result.data[i] = (double *)malloc(n * sizeof(double));
    }
    app_data->results_calculated = TRUE;
    
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "Eigenvalues: ");
    for (int i = 0; i < n; ++i) {
        char val_str[32];
        snprintf(val_str, sizeof(val_str), "%.2f ", eigenvalues[i]);
        strcat(buffer, val_str);
    }
    GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(app_data->parent_window), GTK_DIALOG_MODAL, GTK_MESSAGE_INFO, GTK_BUTTONS_OK, "%s", buffer);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
    
    show_matrix_dialog(app_data->parent_window, eigenvectors);
    app_data->eigenvectors_result = eigenvectors;
    app_data->eigenvectors_result_loaded = TRUE;
}
void on_save_matrix(GtkWidget *widget, gpointer data) {
    AppData *app_data = (AppData *)data;

    if (!app_data->matrix_loaded) {
        printf("No matrix loaded to save\n");
        return;
    }

    GtkWidget *dialog = gtk_file_chooser_dialog_new("Save File", GTK_WINDOW(app_data->parent_window),
        GTK_FILE_CHOOSER_ACTION_SAVE, "_Cancel", GTK_RESPONSE_CANCEL, "_Save", GTK_RESPONSE_ACCEPT, NULL);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        FILE *file = fopen(filename, "w");
        if (file) {
            fprintf(file, "Matrix (%d x %d):\n", app_data->matrix.rows, app_data->matrix.cols);
            for (int i = 0; i < app_data->matrix.rows; ++i) {
                for (int j = 0; j < app_data->matrix.cols; ++j) {
                    fprintf(file, "%.2f ", app_data->matrix.data[i][j]);
                }
                fprintf(file, "\n");
            }
            
             // Save new matrix if loaded
            if (app_data->new_matrix_loaded) {
                Matrix to_save = app_data->new_matrix;
                fprintf(file, " Matrix which user's load to calculate (%d x %d):\n", to_save.rows, to_save.cols);
                for (int i = 0; i < to_save.rows; ++i) {
                    for (int j = 0; j < to_save.cols; ++j) {
                        fprintf(file, "%.2f ", to_save.data[i][j]);
                    }
                    fprintf(file, "\n");
                }
            }
            if (app_data->add_result_loaded) {
                Matrix to_save = app_data->add_result;
                fprintf(file, "Addition Result:\n");
                fprintf(file, "%d %d\n", to_save.rows, to_save.cols);
                for (int i = 0; i < to_save.rows; ++i) {
                    for (int j = 0; j < to_save.cols; ++j) {
                        fprintf(file, "%.2f ", to_save.data[i][j]);
                    }
                    fprintf(file, "\n");
                }
            }
            if (app_data->subtract_result_loaded) {
                Matrix to_save = app_data->subtract_result;
                fprintf(file, "Subtraction Result:\n");
                fprintf(file, "%d %d\n", to_save.rows, to_save.cols);
                for (int i = 0; i < to_save.rows; ++i) {
                    for (int j = 0; j < to_save.cols; ++j) {
                        fprintf(file, "%.2f ", to_save.data[i][j]);
                    }
                    fprintf(file, "\n");
                }
            } 
            if (app_data->multiply_result_loaded) {
                Matrix to_save = app_data->multiply_result;
                fprintf(file, "Multiplication Result:\n");
                fprintf(file, "%d %d\n", to_save.rows, to_save.cols);
                for (int i = 0; i < to_save.rows; ++i) {
                    for (int j = 0; j < to_save.cols; ++j) {
                        fprintf(file, "%.2f ", to_save.data[i][j]);
                    }
                    fprintf(file, "\n");
                }
            } 
            if (app_data->scalar_result_loaded) {
                Matrix to_save = app_data->scalar_result;
                fprintf(file, "Scalar Multiplication Result:\n");
                fprintf(file, "%d %d\n", to_save.rows, to_save.cols);
                for (int i = 0; i < to_save.rows; ++i) {
                    for (int j = 0; j < to_save.cols; ++j) {
                        fprintf(file, "%.2f ", to_save.data[i][j]);
                    }
                    fprintf(file, "\n");
                }
            } 

            if (app_data->results_calculated) {
                // Save determinant
                fprintf(file, "Determinant: %lf\n", app_data->determinant_result);

                // Save inverse
                fprintf(file, "Inverse Matrix:\n");
                for (int i = 0; i < app_data->inverse_result.rows; ++i) {
                    for (int j = 0; j < app_data->inverse_result.cols; ++j) {
                        fprintf(file, "%.2f ", app_data->inverse_result.data[i][j]);
                    }
                    fprintf(file, "\n");
                }

                // Save rank
                fprintf(file, "Rank: %d\n", app_data->rank_result);

                // Save eigenvalues and eigenvectors
                fprintf(file, "Eigenvalues:\n");
                for (int i = 0; i < app_data->matrix.rows; ++i) {
                    fprintf(file, "%.2f ", app_data->eigenvalues_result[i]);
                }
                fprintf(file, "\n");

                fprintf(file, "Eigenvectors:\n");
                for (int i = 0; i < app_data->eigenvectors_result.rows; ++i) {
                    for (int j = 0; j < app_data->eigenvectors_result.cols; ++j) {
                        fprintf(file, "%.2f ", app_data->eigenvectors_result.data[i][j]);
                    }
                    fprintf(file, "\n");
                }
            }
           
            fclose(file);
        }
        g_free(filename);
    }

    gtk_widget_destroy(dialog);
}

void free_matrix(Matrix matrix) {
    for (int i = 0; i < matrix.rows; ++i) {
        free(matrix.data[i]);
    }
    free(matrix.data);
}