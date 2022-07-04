
// This program estimates a battery system's SOC (state of charge) with respect to time
// based off simulated values for the circuit's V (voltage) and I (current).
// The program must read data from SOC_data_sim.csv
// A csv file with all the estimated values named kalman_SOC_estimate.csv is generated.
// In the generated csv, you can compare the actual SOC values to the estimated SOC values.

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stdbool.h>

#pragma warning (disable:4996)


// Function declarations
double VOC(double SOC, double Voc0);
double hk(double SOC_val, double I, double Voc0, double R0);
void fk(double**** xhatk_1, double I, double dt, double Cbat, double Cc, double Rc, double**** fk_);
void mat_multiply(double**** Aprime, double**** Pk_1, double**** sub_mat_1, int num_cols_of_second_matrix);
void transpose(double*** Aprime, double*** Aprime_transpose);
void EKF(double*** xhatk_1, double*** Pk_1, double I, double Ik_1, double V, double Voc0,
	double Rk, double*** Aprime, double** Cprime, double*** Eprime, double*** Fprime,
	double*** fk_, double dt, double Cbat, double Cc, double Rc, double*** Qk1,
	double*** Aprime_transpose, double*** Eprime_transpose,
	double*** sub_mat_1, double*** sub_mat_2, double*** sub_mat_3, double*** xhat,
	double*** P, double*** Cprime_transpose, double*** Lk, double R0,
	double*** xhatCorrected, double*** PCorrected, FILE** fp, double* actualSOC);
double** memory_allocate(int num_columns);
void insertFirst(double t_val, double SOC_actual, double SOC_estimated);
struct node* insertAtEnd(double t_val, double SOC_actual, double SOC_estimated);
int print_list_into_csv(FILE** fp2);



//------------------------------- SMALLER CALCULATION FUNCTIONS --------------------------------

// Functions we need to calculate
double VOC(double SOC, double Voc0) {
	double VOC_ = 0.007 * SOC + Voc0;
	return VOC_;
}

double hk(double SOC_val, double I, double Voc0, double R0) {
	//double SOC_val = (*xhat)[0][0];
	double hk_ = VOC(SOC_val, Voc0) - (R0 * I);		// hk_ is calculated voltage?
	return hk_;
}

void fk(double**** xhatk_1, double I, double dt, double Cbat, double Cc, double Rc, double**** fk_) {
	double my_mat[2][1] = { {(-1) * I / Cbat}, {(I / Cc) - ((**xhatk_1)[1][0] / (Cc * Rc))} };
	for (int i = 0; i < 2; i++) {
		my_mat[i][0] *= dt;
		(**fk_)[i][0] = (**xhatk_1)[i][0] + my_mat[i][0];
	}

	// xhat = previous xhat + change in xhat in time (dt)
	// xhat is a derivative
}


//------------------------------------ MATRIX FUNCTIONS -------------------------------------

// Function to allocate memory for a 2x1 or 2x2 matrix
double** memory_allocate(int num_columns) {
	int i, j;
	// memory allocation
	double** xhat;
	xhat = (double**)malloc((2) * sizeof(double*));
	for (i = 0; i < 2; i++) {
		xhat[i] = (double*)malloc((num_columns) * sizeof(double));
	}
	// Initialize values to 0.0
	for (i = 0; i < 2; i++) {
		for (j = 0; j < num_columns; j++) {
			xhat[i][j] = 0.0;
		}
	}
	return xhat;
}

// Function to multiply 2 matrixes
void mat_multiply(double**** Aprime, double**** Pk_1, double**** sub_mat_1, int num_cols_of_second_matrix) {
	int row, column, index;
	double temp_val = 0.0;
	for (row = 0; row < 2; row++) {
		for (column = 0; column < num_cols_of_second_matrix; column++) {
			for (index = 0; index < 2; index++) {
				temp_val += ((**Aprime)[row][index]) * ((**Pk_1)[index][column]);
			}
			(**sub_mat_1)[row][column] = temp_val;
			temp_val = 0;
		}
	}
}

// Function to transpose a 2x2 matrix
void transpose(double*** Aprime, double*** Aprime_transpose) {
	(*Aprime_transpose)[0][1] = (*Aprime)[1][0];
	(*Aprime_transpose)[0][0] = (*Aprime)[0][0];
	(*Aprime_transpose)[1][0] = (*Aprime)[0][1];
	(*Aprime_transpose)[1][1] = (*Aprime)[1][1];
}



//------------------------------------ EKF CALCULATOR FUNCTION -------------------------------------

// EKF Calculator function (to compute values of xhatCorrected and PCorrected continuously)
void EKF(double*** xhatk_1, double*** Pk_1, double I, double Ik_1, double V, double Voc0,
	double Rk, double*** Aprime, double** Cprime, double*** Eprime, double*** Fprime,
	double*** fk_, double dt, double Cbat, double Cc, double Rc, double*** Qk1,
	double*** Aprime_transpose, double*** Eprime_transpose,
	double*** sub_mat_1, double*** sub_mat_2, double*** sub_mat_3, double*** xhat,
	double*** P, double*** Cprime_transpose, double*** Lk, double R0,
	double*** xhatCorrected, double*** PCorrected, FILE** fp, double* actualSOC) {

	int row, col, index;

	//----------------- CALCULATIONS FOR xhat -----------------
	fk(&xhatk_1, I, dt, Cbat, Cc, Rc, &fk_);
	for (row = 0; row < 2; row++) {
		(*xhat)[row][0] = (*fk_)[row][0];
	}

	//----------------- CALCULATIONS FOR P -----------------
	// Multiplying Aprime * Pk_1 * Aprime_transpose
	mat_multiply(&Aprime, &Pk_1, &sub_mat_1, 2);
	mat_multiply(&sub_mat_1, &Aprime_transpose, &sub_mat_2, 2);

	// Multiplying Eprime * Qk1 * Eprime_transpose
	mat_multiply(&Eprime, &Qk1, &sub_mat_1, 2);
	mat_multiply(&sub_mat_1, &Eprime_transpose, &sub_mat_3, 2);

	// P = (Aprime * Pk_1 * Aprime_transpose) + (Eprime * Qk1 * Eprime_transpose)
	for (row = 0; row < 2; row++) {
		for (col = 0; col < 2; col++) {
			(*P)[row][col] = (*sub_mat_2)[row][col] + (*sub_mat_3)[row][col];
		}
	}

	//----------------- CALCULATIONS FOR Lk -----------------
	// error_estimate = (Cprime * P * Cprime_tranpose)
	mat_multiply(&P, &Cprime_transpose, &sub_mat_1, 1);
	double error_estimate = 0.0, combined_error = 0.0;
	for (row = 0; row < 2; row++) {
		error_estimate += (*sub_mat_1)[row][0] * (*Cprime)[row];
	}
	// combined_error = [(Cprime * P * Cprime_tranpose) + Rk]^(-1)
	combined_error = error_estimate + Rk;
	combined_error = pow(combined_error, (-1));

	// Lk = combined_error * (P * Cprime_transpose)
	for (row = 0; row < 2; row++) {
		(*Lk)[row][0] = (*sub_mat_1)[row][0] * combined_error;
	}

	//--------------- CALCULATIONS FOR xhatCorrected ---------------
	// temp_val = (SOC + Vc) - hk(100*SOC, I, Voc0, R0)
	double temp_val = 0.0;
	temp_val = (V + (*xhat)[1][0]) - hk(100 * (*xhat)[0][0], I, Voc0, R0);
	
	// xhatCorrected = (Lk * temp_val) * xhat
	for (row = 0; row < 2; row++) {
		(*sub_mat_1)[row][0] = (*Lk)[row][0] * temp_val;
		(*xhatCorrected)[row][0] = (*xhat)[row][0] + (*sub_mat_1)[row][0];
	}

	//----------------- CALCULATIONS FOR PCorrected -----------------
	// FINDING: Lk * Cprime
	for (row = 0; row < 2; row++) {
		for (index = 0; index < 2; index++) {
			(*sub_mat_1)[row][index] = (*Lk)[row][0] * (*Cprime)[index];
		}
	}
	// FINDING: (Lk * Cprime) * P
	mat_multiply(&sub_mat_1, &P, &sub_mat_2, 2);

	// PCorrected = P - [(Lk * Cprime) * P]
	for (row = 0; row < 2; row++) {
		for (col = 0; col < 2; col++) {
			(*PCorrected)[row][col] = (*P)[row][col] - (*sub_mat_2)[row][col];
		}
	}

}


//---------------------------------- LINKED LIST FUNCTIONS -----------------------------------

// Linked list struct
struct node {
	double t;
	double SOC_actual;
	double SOC_estimated;
	struct node* next;
};

// WE MUST HAVE A POINTER TO HEAD OF LIST AND POINTER TO CURRENT ELEM
struct node* head = NULL;
struct node* current = NULL;


// INSERT LINK AT THE FIRST LOCATION
void insertFirst(double t_val, double SOC_actual, double SOC_estimated) {
	//create a link
	struct node* link = (struct node*)malloc(sizeof(struct node));

	link->t = t_val;
	link->SOC_actual = SOC_actual;
	link->SOC_estimated = SOC_estimated;

	current = head = link;
	link->next = NULL;
}


// Insert link at end
struct node* insertAtEnd(double t_val, double SOC_actual, double SOC_estimated) {
	//create a link
	struct node* link = (struct node*)malloc(sizeof(struct node));
	link->t = t_val;
	link->SOC_actual = SOC_actual;
	link->SOC_estimated = SOC_estimated;
	link->next = NULL;

	current->next = link;
	current = link;
	return link;
}

// Print the linked list into csv file
int print_list_into_csv(FILE** fp2) {
	struct node* ptr = head;
	while (ptr != NULL) {
		fprintf(*fp2, "%.2Lf,%.6Lf,%.6Lf\n", ptr->t, ptr->SOC_actual, ptr->SOC_estimated);
		ptr = ptr->next;
	}
	return 0;
}



//--------------------------------------- DRIVER CODE ----------------------------------------

int main() {

	// Constants
	double dt = 0.01;		// Sampling Period
	int totalTime = 100;		// total #seconds
	double R0 = 0.01;
	double Rc = 0.015;
	int Cc = 2400;
	int Cbat = 18000;
	double Voc0 = 3.435;
	double alp = 0.007;

	// Opening the csv files
	FILE* fp;
	fp = fopen("SOC_data_sim.csv", "r");
	if (fp == NULL) {
		printf("SOC_data_sim.csv can't be opened");
		return 1;
	}
	FILE* fp2;
	fp2 = fopen("kalman_SOC_estimate.csv", "w");
	if (fp2 == NULL) {
		printf("kalman_SOC_estimate.csv can't be opened");
		return 1;
	}


	// Initializing probability matricies
	double** Aprime = memory_allocate(2);
	double** Aprime_transpose = memory_allocate(2);
	Aprime[0][0] = 1.0;
	Aprime[1][1] = exp(-dt / (Cc * Rc));
	transpose(&Aprime, &Aprime_transpose);
	
	double** Eprime = memory_allocate(2);
	double** Eprime_transpose = memory_allocate(2);
	Eprime[0][0] = 1.0;
	Eprime[1][1] = 1.0;
	transpose(&Eprime, &Eprime_transpose);

	double** Fprime = memory_allocate(2);
	Fprime[0][0] = 1.0;
	Fprime[1][1] = 1.0;

	double* Cprime;					// VOC(SOC)
	Cprime = (double*)malloc(2 * sizeof(double));
	Cprime[0] = 0.007;
	Cprime[1] = -1.0;

	double** Cprime_transpose;
	Cprime_transpose = memory_allocate(1);
	(Cprime_transpose)[0][0] = Cprime[0];
	(Cprime_transpose)[1][0] = Cprime[1];

	// Coefficients for probability
	double Rk = pow(10, -4);


	// --------------------------------------------------------------------------
	// At time k
	// --------------------------------------------------------------------------

	// Variables Qualitatively
	// (Constant)VOC: Open circuit voltage
	// (Time variant) Vc : Voltage across capacitor in battery circuit model
	// (Time variant) V : Voltage we measure(output from battery circuit model)
	// (Constant)Cc : Capacitor's capacitance
	// (Constant)Rc : Resistor in parallel with capacitor
	// (Constant)R0 : Ohmic resistance
	// (Constant)Cbat : Capacity of battery in AmpHours ?


	//// Initializing the matrixes we'll need to use
	double** xhat = memory_allocate(1);	//xhat is a 2-by-1 matrix, top value is: SOC_estimated, bottom is Vc_estimated
	double** Qk1 = memory_allocate(2);
	Qk1[0][0] = 2.5 * pow(10, -7);
	double** xhatk_1 = memory_allocate(1);
	double** fk_ = memory_allocate(1);
	double** xhatCorrected = memory_allocate(1);
	double** PCorrected = memory_allocate(2);
	double** Lk = memory_allocate(1);
	double** P = memory_allocate(2);
	double** Pk_1 = memory_allocate(2);
	P[0][0] = Rk;
	P[1][1] = Rk;
	Pk_1[0][0] = Rk;
	Pk_1[1][1] = Rk;
	double** sub_mat_1 = memory_allocate(2);
	double** sub_mat_2 = memory_allocate(2);
	double** sub_mat_3 = memory_allocate(2);

	// Variables subject to constant change
	int i, j;
	double t = 0.0, actualSOC = 1.0, Vc = 0.0;
	double I = 0.0;
	double Ik_1 = 0.0;
	double V = 0.0;		// current V
	(xhatk_1)[0][0] = actualSOC;
	(xhatk_1)[1][0] = Vc;
	//printf("xhatCorrected:\t\tPCorrected:\n");

	// Linked List Control Variables
	bool start_new_list = true;
	bool done_csv_printing = false;
	int while_loop_iter_count = 0;
	int iteration_count_50 = 50;
	fprintf(fp2, "Time(s),Actual_State_Of_Charge,Estimated_State_Of_Charge\n");



	//THE BIG FOR LOOP TO CONTINUOUSLY CALCULATE xhatCorrected & PCorrected
	while (!feof(fp)) {
		// read values for t, SOC, Vc, V, I from csv file
		fscanf(fp, "%lf,%lf,%lf,%lf,%lf", &t, &actualSOC, &Vc, &V, &I);

		// skip to second line of csv after first is read
		if (t == 0.0) {
			t += dt;
			continue;
		}

		// FUNCTION TO CALCULATE xhatCorrected & PCorrected
		EKF(&xhatk_1, &Pk_1, I, Ik_1, V, Voc0, Rk, &Aprime, &Cprime, &Eprime,
			&Fprime, &fk_, dt, Cbat, Cc, Rc, &Qk1,
			&Aprime_transpose, &Eprime_transpose, &sub_mat_1, &sub_mat_2,
			&sub_mat_3, &xhat, &P, &Cprime_transpose, &Lk, R0, &xhatCorrected,
			&PCorrected, &fp, &actualSOC);


		// Setting the xhatk_1 and Pk_1 values
		for (i = 0; i < 2; i++) {
			// set values of xhatk_1 to xhatCorrected
			xhatk_1[i][0] = xhatCorrected[i][0];
			// For each column of second matrix
			for (j = 0; j < 2; j++) {
				// set values of Pk_1 to PCorrected
				Pk_1[i][j] = PCorrected[i][j];
			}
		}

		//// PRINTING MATRIXES TO CONSOLE
		//printf("[");
		//// For each row of matrix
		//for (i = 0; i < 2; i++) {
		//	// set values of xhatk_1 to xhatCorrected
		//	xhatk_1[i][0] = xhatCorrected[i][0];
		//	printf("%.6lf", xhatCorrected[i][0]);
		//	if (i == 1)
		//		printf("]\t\t");
		//	else
		//		printf("\t\t[");

		//	// For each column of second matrix
		//	for (j = 0; j < 2; j++) {
		//		// set values of Pk_1 to PCorrected
		//		Pk_1[i][j] = PCorrected[i][j];
		//		printf("%.6lf", PCorrected[i][j]);
		//		if (j == 0)
		//			printf("  ");
		//	}
		//	if (i == 1)
		//		printf("]");
		//	printf("\n");
		//}
		//printf("\n");


		// Inserting the (t,SOC_actual,SOC_estimated) into linked list
		if (start_new_list == true) {
			insertFirst(t, actualSOC, xhatCorrected[0][0]);
			start_new_list = false;
		}
		else {
			insertAtEnd(t, actualSOC, xhatCorrected[0][0]);
		}

		// Inserting the linked list into the csv
		if (while_loop_iter_count % iteration_count_50 == 0) {
			print_list_into_csv(&fp2);
			if (t == totalTime) {
				done_csv_printing = true;
			}
			start_new_list = true;
			current = head;
		}

		while_loop_iter_count += 1;
	}

	// print last few elements into csv after while loop is over
	if (done_csv_printing == false)
		print_list_into_csv(&fp2);

	printf("DONE generating kalman_SOC_estimate.csv\n");
	// Close the csv's
	fclose(fp);
	fp = NULL;
	fclose(fp2);
	fp2 = NULL;
	return 0;
}

