#include <glpk.h>
#include <support.c>
#include "tinyexpr.c"

#define MAX_LINE 2048
#define MAX_PARS 100
#define MAX_FPARS 100
#define MAX_EQNS 100
#define MAX_OPTIONS 100

typedef struct {
    char name[50];
    double value;
} Parameter;
Parameter parameters[MAX_PARS], free_parameters[MAX_FPARS];

#define opt_count 5
Parameter options[opt_count]={{"t0\0",0.0},{"t1\0",10000},{"dt\0",0.1}, {"tshow\0",1.0}, {"LP_tolerance", 1e-18}};
double t0 = 0.0, t1 = 10000.0, dt = 0.1, tshow=1.0, LP_tolerance=1e-18;

typedef struct {
    char name[50];
    char expression[MAX_LINE];
} Intermediate;
Intermediate intermediates[MAX_PARS];
int intermediate_count = 0;
char line[MAX_LINE];
int par_count = 0, fpar_count = 0, eqn_count = 0, var_count = 0;
FILE *file, *fr;
char eqns[MAX_EQNS][MAX_LINE];
char vars_eqns[MAX_EQNS][MAX_LINE];
char vars[MAX_EQNS][MAX_LINE];
char init_conditions[MAX_EQNS][MAX_LINE];
double dinit_cond[MAX_EQNS];
te_variable pars[MAX_PARS];
te_variable pars_temp[MAX_PARS];
double **Amatrix, *beqn;



void removeSpaces(char *str) {
    int i = 0, j = 0;
    // Traverse the string
    while (str[i]) {
        // If the current character is not a space, copy it to the output
        if (str[i] != ' ') {
            str[j++] = str[i];
        }
        i++;
    }
    // Null-terminate the modified string
    str[j] = '\0';
}

// Function to substitute variables with expressions
int substituteExpression(char *equation) {
    char temp[MAX_LINE];
    strcpy(temp, equation);
	int isrepl=0;
    for (int i = 0; i < intermediate_count; i++) {
        char *pos = temp;
        while ((pos = strstr(pos, intermediates[i].name))) {
            // Check if it's an exact match (not part of a larger name)
            char before = (pos == temp) ? '\0' : *(pos - 1);
            char after = *(pos + strlen(intermediates[i].name));

            if ((before == '\0' || before == '+' || before == '-' || before == '*' || before == '/' || before == '(') &&
                (after == '\0' || after == '+' || after == '-' || after == '*' || after == '/' || after == ')')) {
                // Replace variable with its expression
                char before_str[MAX_LINE] = "";
                char after_str[MAX_LINE] = "";
				isrepl=1;
                // Copy parts of the equation before and after the match
                strncpy(before_str, temp, pos - temp);
                before_str[pos - temp] = '\0';
                strcpy(after_str, pos + strlen(intermediates[i].name));

                // Construct the new equation
                snprintf(temp, MAX_LINE, "%s(%s)%s", before_str, intermediates[i].expression, after_str);

                // Update position to continue searching after this substitution
                pos = temp + strlen(before_str) + strlen(intermediates[i].expression) + 2; // Account for added parentheses
            } else {
                // Skip to the next character to avoid infinite loop
                pos++;
            }
        }
    }
    // Copy back the substituted expression
    strcpy(equation, temp);
	if(isrepl==1) return 1; else return 0;
}


// Helper function to check if a substring is a valid scientific notation
int isScientificNotation(const char *str) {
    const char *p = str;
    // Optional leading sign
    if (*p == '+' || *p == '-') p++;
    // Must start with a digit or a dot
    if (!isdigit(*p) && *p != '.') return 0;
    // Digits before the dot
    while (isdigit(*p)) p++;
    // Optional dot and fractional part
    if (*p == '.') {
		p++;
        while (isdigit(*p)) {
        p++;
        }
    }
    // Optional exponent part
    if (*p == 'e' || *p == 'E') {
        p++;
        // Optional sign after exponent
        if (*p == '+' || *p == '-') {
            p++;
        }
        // At least one digit must follow the exponent
        if (!isdigit(*p)) {
            return 0;
        }
        while (isdigit(*p)) {
            p++;
        }
    }
    // The entire string must be consumed for a valid match
    return *p == '\0';
}

// Function to check if a string contains any operators (+, -, *, /)
int containsOperators(const char *str) {
    return (strpbrk(str, "+-*/") != NULL);
}

void print_data(void) {
	int i;
	 fprintf(fr, "INPUT DATA\n");
	// Print eqns
	  fprintf(fr,"\nEquations:\n");
	 for (int i = 0; i < eqn_count; i++) fprintf(fr,"d%s/dt=%s\n", vars[i], eqns[i]);
	 // Print eqns
	  fprintf(fr,"\nInitial Conditions:\n");
	 for (int i = 0; i < eqn_count; i++) fprintf(fr,"%s(0)=%s=%e\n", vars[i], init_conditions[i],dinit_cond[i] );
		double *temp;
	// Print the parameters
    fprintf(fr, "\nParameters:\n");
    for (int i = 0; i < par_count-var_count-fpar_count; i++) {
	   temp=pars[i].address;
       fprintf(fr, "%s = %e\n", pars[i].name, *temp);
    }
	 // Print the free parameters
    fprintf(fr,"\nFree Parameters:\n");
    for (int i = 0; i < fpar_count; i++) fprintf(fr,"%s\n", free_parameters[i].name);
	// Print the intermediates
	  fprintf(fr,"\nIntermediates:\n");
	 for (int i = 0; i < intermediate_count; i++) {
             fprintf(fr,"%s=%s\n", intermediates[i].name, intermediates[i].expression);
     }
	  // Print the options
    fprintf(fr,"\nOptions:\n");
    for (int i = 0; i < opt_count; i++) {
        fprintf(fr,"%s=%e\n", options[i].name,options[i].value );
    }
}

void print_data_console(void) {	
	int i;
	 printf( "INPUT DATA\n");
	// Print eqns
	  printf("\nEquations:\n");
	 for (int i = 0; i < eqn_count; i++) printf("d%s/dt=%s\n", vars[i], eqns[i]);
	 // Print eqns
	  printf("\nInitial Conditions:\n");
	 for (int i = 0; i < eqn_count; i++) printf("%s(0)=%s=%e\n", vars[i], init_conditions[i],dinit_cond[i] );
		double *temp;
	// Print the parameters
    printf("\nParameters:\n");
    for (int i = 0; i < par_count; i++) {
	   temp=pars[i].address;
       printf( "%s = %e\n", pars[i].name, *temp);
    }
	 // Print the free parameters
    printf("\nFree Parameters:\n");
    for (int i = 0; i < fpar_count; i++) printf("%s\n", free_parameters[i].name);
	// Print the intermediates
	  printf("\nIntermediates:\n");
	 for (int i = 0; i < intermediate_count; i++) printf("%s=%s\n", intermediates[i].name, intermediates[i].expression);
   // Print the options
    printf("\nOptions:\n");
    for (int i = 0; i < opt_count; i++) printf("%s=%e\n", options[i].name,options[i].value );
}

void evaluate_in_conds(void) {
	int err;
	for (int i = 0; i < var_count; i++) {
        te_expr *expr = te_compile(init_conditions[i], pars, par_count, &err);
		if (expr) {
			dinit_cond[i] = te_eval(expr);
			te_free(expr);
		} else {
			printf("In assess Initial Condition: Error in equation at position %d: %s\n", err, eqns[i]);
		}
    }
}

void assign_options(void) {
	for(int i=0;i<opt_count;i++) {
		if(strcmp("t0", options[i].name)==0) t0=options[i].value;
		if(strcmp("t1", options[i].name)==0) t1=options[i].value;
		if(strcmp("dt", options[i].name)==0) dt=options[i].value;
		if(strcmp("tshow", options[i].name)==0) tshow=options[i].value;
		if(strcmp("LP_tolerance", options[i].name)==0) tshow=options[i].value;
	}
}
int load_file(void) {
	int err;
	double result;
	char c;
    file = fopen("model.txt", "r");
    if (!file) {
        perror("Error opening file");
        return 0;
    }
    // Loop through each line in the file
    while (fgets(line, MAX_LINE, file)) {
        line[strcspn(line, "\n")] = 0; // Remove newline
        if (strlen(line) == 0 || line[0] == '#') continue; // Skip empty lines or comments
        removeSpaces(line); // Remove all spaces from the line
        
        char *equal_sign = strchr(line, '=');
        if (equal_sign) {
            *equal_sign = '\0';
            char *name = line;
            char *value_str = equal_sign + 1;
            // Check if the line specifies an initial condition
            if (strstr(name, "(0)")) {
                strncpy(vars[var_count], name, strchr(name, '(') - name); // Extract variable name before "(0)"
                strcpy(init_conditions[var_count], value_str); // Extract variable name before "(0)"
			  	var_count++;
			  } else if (strncmp(line, "opt", 3) == 0) {
				// Options
				char stemp[MAX_LINE];
				strncpy(stemp, line + 3, sizeof(stemp) - 3);
				for(int i=0;i<opt_count;i++) {
					if(strcmp(stemp, options[i].name)==0) {
						options[i].value = atof(value_str);
					}
				}
			}
            // Check if the line defines a differential equation
            else if (strstr(line, "d") && strstr(line, "/dt")) {
                char *variable_part = strtok(line, "/"); // Extract variable part (e.g., "dx/dt" -> "dx")
				variable_part++;
				 if (variable_part && value_str) {
                    // Extract variable name (skip the leading 'd' in "dx")
					strncpy(vars_eqns[eqn_count], variable_part, sizeof(vars_eqns[eqn_count]) - 1);
                    strncpy(eqns[eqn_count], value_str, sizeof(eqns[eqn_count]) - 1);
                    eqn_count++;
                }
            }
            // Check if the line defines an intermediate expression
            else if (containsOperators(value_str) && !isScientificNotation(value_str)) {
                strncpy(intermediates[intermediate_count].name, name, sizeof(intermediates[intermediate_count].name) - 1);
                strncpy(intermediates[intermediate_count].expression, value_str, sizeof(intermediates[intermediate_count].expression) - 1);
                intermediate_count++;
            }
            // Check if the line defines a constant parameter
            else if (isScientificNotation(value_str)) {
                strncpy(parameters[par_count].name, name, sizeof(parameters[par_count].name) - 1);
                parameters[par_count].value = atof(value_str);
                pars[par_count].name = parameters[par_count].name;
                pars[par_count].address = &parameters[par_count].value;
                par_count++;
            } 
            // Handle unexpected values or invalid lines
            else {
                printf("Unexpected value: %s\n", value_str);
            }
        } else if (strncmp(line, "free", 4) == 0) {
            // Free parameter
            strncpy(free_parameters[fpar_count].name, line + 4, sizeof(free_parameters[fpar_count].name) - 1);
            free_parameters[fpar_count].value = atof("NaN");
            fpar_count++;
            // A free parameter is also a parameter
            strncpy(parameters[par_count].name, line + 4, sizeof(parameters[par_count].name) - 1);
            parameters[par_count].value = atof("NaN");
            pars[par_count].name = parameters[par_count].name;
            pars[par_count].address = &parameters[par_count].value;
            par_count++;
        } 
    }
	 fclose(file);
    evaluate_in_conds();
	assign_options();
	for (int i = 0; i < var_count; i++){
		pars[par_count].name = vars[i];
        pars[par_count].address = &dinit_cond[i];
		double *temp2=pars[par_count].address;
		par_count+=1;
	}
    print_data();
	// Substitute intermediate expressions in the equations
    for (int i = 0; i < eqn_count; i++) {
        substituteExpression(eqns[i]);
		substituteExpression(eqns[i]);
    }
    return 1;
}


// Function to evaluate equations
double evaluate(const char *equation, double *y) {
	double temp;
	int err;
	for (int i = 0; i < par_count; i++) {
		for (int j = 0; j < eqn_count; j++) {
			if(strcmp(pars[i].name,vars[j])==0) {
				pars[i].address=&y[j];
			}
		}
	}
	te_expr *expr = te_compile(equation, pars, par_count, &err);
	if (expr) {
		temp = te_eval(expr);
		te_free(expr);
	} else {
		printf("In assess Initial Condition: Error in equation at position %d: %s\n", err, equation);
	}
    return temp; // Placeholder (you need to replace this with your logic)
}

// Runge-Kutta 4th Order Step
void runge_kutta_step(double t, double dt, double *y, int num_vars, double *k1, double *k2, double *k3, double *k4, char eqns[MAX_EQNS][MAX_LINE]) {
    double y_temp[MAX_EQNS];
    // Calculate k1
    for (int i = 0; i < num_vars; i++) k1[i] = dt * evaluate(eqns[i], y);
    // Calculate k2
    for (int i = 0; i < num_vars; i++) y_temp[i] = y[i] + 0.5 * k1[i];
    for (int i = 0; i < num_vars; i++) k2[i] = dt * evaluate(eqns[i], y_temp);
    // Calculate k3
    for (int i = 0; i < num_vars; i++) y_temp[i] = y[i] + 0.5 * k2[i];
    for (int i = 0; i < num_vars; i++) k3[i] = dt * evaluate(eqns[i], y_temp);
    // Calculate k4
    for (int i = 0; i < num_vars; i++) y_temp[i] = y[i] + k3[i];
    for (int i = 0; i < num_vars; i++) k4[i] = dt * evaluate(eqns[i], y_temp);
    // Update y values
    for (int i = 0; i < num_vars; i++) y[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
}

// Runge-Kutta Solver
void solve_odes_runge_kutta(double t0, double t1, double dt, double *y_initial, int num_vars, char eqns[MAX_EQNS][MAX_LINE]) {
    double t = t0, ttemp=0;
    double y[MAX_EQNS];
    double k1[MAX_EQNS], k2[MAX_EQNS], k3[MAX_EQNS], k4[MAX_EQNS];
    // Copy initial values
    for (int i = 0; i < num_vars; i++) y[i] = y_initial[i];
    // Open output file
    FILE *output_file = fopen("output.txt", "w");
    if (!output_file) {
        perror("Error opening output file");
        return;
    }
    // Write header
    fprintf(output_file, "time");
    for (int i = 0; i < num_vars; i++) {
        fprintf(output_file, "\t%s", vars[i]);
    }
    fprintf(output_file, "\n");
    // Time-stepping loop
    while (t <= t1) {
        // Write current values to file
		if(ttemp>tshow) {
			fprintf(output_file, "%.5f", t);
			for (int i = 0; i < num_vars; i++) {
				fprintf(output_file, "\t%e", y[i]);
			}
			fprintf(output_file, "\n");
			ttemp=0.0;
		}
        // Perform a single Runge-Kutta step
        runge_kutta_step(t, dt, y, num_vars, k1, k2, k3, k4, eqns);
        // Update time
        t += dt;
		ttemp += dt;
    }
    fclose(output_file);
}

int isfree(char *ss) {
	for(int i=0;i<fpar_count;i++) if(strcmp(free_parameters[i].name, ss)==0) return 1;
	return 0;
}

void assess_Amatrix(void) {
	int i,j, ieqn, ifpar, err;
	double temp,temp2, result;
	Amatrix=matrix(0,eqn_count-1,0,fpar_count-1);
	beqn=vector(0,eqn_count-1);
	// assess b
	// fill pars_temp
	for (int i = 0; i < par_count; i++) {
		pars_temp[i].name= pars[i].name;
		temp=0.0;
		if(isfree(pars[i].name)) pars_temp[i].address=&temp; else pars_temp[i].address=pars[i].address;
    }
   for (int i = 0; i < eqn_count; i++) {
        te_expr *expr = te_compile(eqns[i], pars_temp, par_count, &err);
        if (expr) {
           result = te_eval(expr);
           beqn[i]=-result;
           te_free(expr);
        } else {
            printf("In assess b: Error in equation at position %d: %s\n", err, eqns[i]);
        }
    }
	printf("\nb vector");
	for (int i = 0; i < eqn_count; i++) printf("\n%e", beqn[i]);
	printf("\n");
	//assess Amatrix
	for(ieqn=0;ieqn<eqn_count; ieqn++) {
		for(int ifpar = 0; ifpar < fpar_count; ifpar++) {
			for (int i = 0; i < par_count; i++) {
					pars_temp[i].name= pars[i].name;
					pars_temp[i].address=pars[i].address;
					temp=0.0;
					temp2=1.0;
					if(isfree(pars[i].name))pars_temp[i].address=&temp;
					if(strcmp(pars[i].name, free_parameters[ifpar].name)==0) pars_temp[i].address=&temp2; 
			}
			te_expr *expr = te_compile(eqns[ieqn], pars_temp, par_count, &err);
			if (expr) {
				result = te_eval(expr);
				Amatrix[ieqn][ifpar]=result+beqn[ieqn];
				te_free(expr);
			} else {
				printf("In assess Matrix: Error in equation at position %d: %s\n", err, eqns[i]);
			}
   		}
	}
	printf("\nAmatrix\n");
	for(ieqn=0;ieqn<eqn_count; ieqn++) {
		for(int ifpar = 0; ifpar < fpar_count; ifpar++) {
			printf("%e\t", Amatrix[ieqn][ifpar]);
		}
		printf("\n");
   	}
	printf("\n");
}

void assess_params(void) {
	int i,j;
    glp_prob *lp;
	int num_constraints=eqn_count;
	int num_variables=fpar_count;	
	int ia[MAX_EQNS], ja[MAX_PARS];  // Indices for constraints matrix (adjust size as needed)
    double ar[MAX_EQNS*MAX_PARS];        // Values for constraints matrix
    double min_values[100], max_values[100];
	int issolmin[100], issolmax[100];
	for(i=0;i<100;i++) {
		min_values[i]=0;
		max_values[i]=0;
		ia[i]=0;
		ja[i]=0;
		ar[i]=0;
		issolmin[i]=0;
		issolmax[i]=0;
	}
	lp = glp_create_prob();           // Create a problem instance
    glp_set_prob_name(lp, "linear_system");
    // Add constraints (rows)
    glp_add_rows(lp, num_constraints);
    for (int i = 1; i <= num_constraints; i++) {
        glp_set_row_bnds(lp, i, GLP_FX, beqn[i-1], beqn[i-1]); // Equality constraint
    }
    // Add variables (columns)
    glp_add_cols(lp, num_variables);
	for(i=0;i<num_variables;i++) glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
   	// Define constraint matrix
    int k = 1;
    for (int i = 0; i < num_constraints; i++) {
        for (int j = 0; j < num_variables; j++) {
        	ia[k] = i + 1;  // Row index
        	ja[k] = j + 1;  // Column index
      	    ar[k] = Amatrix[i][j];
      	    k++;
        }
    }
    glp_load_matrix(lp, k - 1, ia, ja, ar);
	// Enable automatic scaling
    glp_scale_prob(lp, GLP_SF_AUTO);
	// Set up simplex parameters with custom tolerances
    glp_smcp params;
    glp_init_smcp(&params);
    params.tol_bnd = LP_tolerance;  // Feasibility tolerance
    params.tol_bnd = LP_tolerance;  // Optimality tolerance
    // Solve for each variable's min and max
	for (int var = 1; var <= num_variables; var++) {
        glp_set_obj_dir(lp, GLP_MIN);	// Minimize x_var
        for (int j = 1; j <= num_variables; j++) 
			glp_set_obj_coef(lp, j, (j == var) ? 1.0 : 0.0);
		glp_simplex(lp, &params);
		issolmin[var - 1]=glp_get_status(lp); 
		min_values[var - 1] = glp_get_obj_val(lp);
		if(min_values[var - 1] < 0.0) min_values[var - 1] = 0.0;
		printf("min=%e\t%e\n",glp_get_obj_val(lp), min_values[var - 1]);
        
        glp_set_obj_dir(lp, GLP_MAX); // Maximize x_var
		glp_simplex(lp, &params);
		issolmax[var - 1]=glp_get_status(lp);
		max_values[var - 1] = glp_get_obj_val(lp);
		printf("max=%e\t%e\n",glp_get_obj_val(lp),max_values[var - 1]);
    }
    // Print results
	fprintf(fr, "\nOUTPUT DATA\nResults");
	for (int i = 0; i < num_variables; i++) {
		fprintf(fr, "\n%s: ", free_parameters[i].name);
        if(issolmin[i]==GLP_OPT) printf(" min = %e, ", min_values[i]);
		if(issolmin[i]==GLP_UNBND) printf(" min = 0.0 ");
		if(issolmin[i]==GLP_INFEAS || issolmin[i]==GLP_NOFEAS || issolmin[i]==GLP_UNDEF) printf(" min = no solution found ");
		if(issolmin[i]==GLP_OPT) fprintf(fr, " min = %e, ", min_values[i]);
		if(issolmin[i]==GLP_UNBND) fprintf(fr, " min = 0.0 ");
		if(issolmin[i]==GLP_INFEAS || issolmin[i]==GLP_NOFEAS || issolmin[i]==GLP_UNDEF) fprintf(fr, " min = no solution found ");
		if(issolmax[i]==GLP_OPT) printf("max = %e", max_values[i]); 
		if(issolmax[i]==GLP_UNBND) printf("max = +inf ");
		if(issolmax[i]==GLP_INFEAS || issolmax[i]==GLP_NOFEAS || issolmax[i]==GLP_UNDEF) printf(" min = no solution found ");
		if(issolmax[i]==GLP_OPT) fprintf(fr,"max = %e", max_values[i]); 
		if(issolmax[i]==GLP_UNBND) fprintf(fr,"max = +inf ");
		if(issolmax[i]==GLP_INFEAS || issolmax[i]==GLP_NOFEAS || issolmax[i]==GLP_UNDEF) fprintf(fr, " min = no solution found ");
	}
    glp_delete_prob(lp);  // Clean up
}


int main(void) {
	char c;
	fr=fopen("results.txt", "w");
	load_file();
	for(int i=0; i<eqn_count;i++) {
		if(strcmp(vars[i], vars_eqns[i])!=0) {
			printf("\nError in defining diff eqns and in conditions");
			scanf("\n", &c);
			return 0;
		}
	}
   double y_initial[MAX_EQNS]; // Initialize variables (replace with actual values)
	for (int j = 0; j < eqn_count; j++) {
		y_initial[j]=dinit_cond[j];
	}
    print_data_console();
	
	printf("opt_count=%d", opt_count);
    // Solve ODEs
    if(fpar_count==0) {
		solve_odes_runge_kutta(t0, t1, dt, y_initial, eqn_count, eqns);
	} else {
		assess_Amatrix();
		assess_params();
	}
	fclose(fr);
	
    return 1;
}
