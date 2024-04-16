#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lml.h"

typedef struct {
    int num_el;
    int num_no;
    int num_dof;
    
    double E;
    double I;
    double L;
    
    Matrix *force_vector;
    Matrix *fixed_vector;
    Matrix *disps_vector;

    Matrix *GSM, *LSM;
} Beam;

// calculates the local stiffness matrix of each element in a beam
// params: b - beam (see Beam struct)
void get_lsm(Beam *b) {
    double L = b->L;
    double arr[4][4] = {{   12,   6*L,  -12,   6*L },
                        {  6*L, 4*L*L, -6*L, 2*L*L },
                        {  -12,  -6*L,   12,  -6*L },
                        {  6*L, 2*L*L, -6*L, 4*L*L }};
    b->LSM = matrix_from_array(4, 4, arr);
    scale(b->LSM, b->E*b->I/pow(L, 3));
}

// assembles the LSMs of each element into global stiffness matrix
// params: b - beam (see Beam struct)
void get_gsm(Beam *b) {
    b->GSM = zeros(b->num_dof, b->num_dof);
    for (int i = 0; i < b->num_el; i++) {
        Matrix *sub = get_submatrix(b->GSM, 2*i, 2*i, 4, 4);
        sub = add(sub, b->LSM); 
        set_submatrix(b->GSM, 2*i, 2*i, sub);
    }
}

// shows the table of the beams constraints
// params: b - beam (see Beam struct)
void show_constraint_table(Beam *b) {
    printf("\nHere are your current constraints:\n\n");

    // print headings
    printf("|   Node   | Position |  Force   |  Moment  | Freedom  |   Disp.  |   Def.   |\n");
    for (int i = 0; i < 7; i++) { printf("|----------"); } printf("|");

    // print the valuse for each of the above headings
    for (int i = 0; i < b->num_no; i++) {
        printf("\n|");
        printf("%9d |", i);
        printf("%9.2f |", i*(b->L/b->num_el));
        printf("%9.2f |", b->force_vector->data[2*i    ][0] != 0 ? b->force_vector->data[2*i    ][0] / 1000 : 0);
        printf("%9.2f |", b->force_vector->data[2*i + 1][0] != 0 ? b->force_vector->data[2*i + 1][0] / 1000 : 0);

        // work out 'freedom of a node'
        if (b->fixed_vector->data[2*i][0] != 0) {
            if (b->fixed_vector->data[2*i + 1][0] != 0) { printf("    fixed |"); }
            else { printf("   pinned |"); }
        } else { printf("     free |"); }

        printf("%9.2f |", b->disps_vector->data[2*i    ][0] != 0 ? b->disps_vector->data[2*i    ][0] * 1000 : 0);
        printf("%9.2f |", b->disps_vector->data[2*i + 1][0] != 0 ? b->disps_vector->data[2*i + 1][0] * 1000 : 0);
    } printf("\n");
}

// allows the user to edit the beam's fixpoints and loads
// params: b - beam (see Beam struct)
void get_beam_info(Beam *b) {

    printf("\nLength of beam (m): "); scanf("%lf", &b->L);
    printf("Number of elements: "); scanf("%d", &b->num_el);
    printf("Youngs Modulus (GPa): "); scanf("%lf", &b->E);
    printf("Moment of inertia (cm4): "); scanf("%lf", &b->I);

    b->num_no = b->num_el + 1;
    b->num_dof = b->num_no * 2;
    b->E *= pow(10, 9); // gpa = 10^9 pa
    b->I *= pow(10, -8); // m4 = 10^-8 cm4

    int fix_count = 0;
    int load_count = 0;

    b->force_vector = zeros(b->num_dof, 1);
    b->disps_vector = zeros(b->num_dof, 1);
    b->fixed_vector = zeros(b->num_dof, 1);

    char adding_conditions = 'y';
    while (adding_conditions == 'y') {
        show_constraint_table(b);

        // get which constraint
        int constraint = 0;
        printf("\nWhat would you like to do?\n");
        printf("\n\t1. Add Force    2. Add Moment");
        printf("\n\t3. Add Pin      4. Add Fixpoint");
        printf("\n\t5. Finish and analyse\n");
        do {
            printf("\nEnter number 1-5: ");
            scanf("%d", &constraint); while (getchar() != '\n');
        } while (constraint < 1 || constraint > 5);

        if (constraint == 5){
            // break loop
            if (fix_count < 2) { printf("\nYou need to add more constraints."); }
            else if (load_count < 1) { printf("\nYou need to add more loads."); }
            else { adding_conditions = 'n'; }
        }
        else{
            // get which node
            int node = -1;
            printf("\nOn which node?\n");
            do {
                printf("\nEnter number 0-%d: ", b->num_no - 1); 
                scanf("%d", &node); while (getchar() != '\n');
            } while (node < 0 || node > b->num_no - 1);
        
            // add constraint
            switch(constraint){
                case 1:
                    printf("\nForce (kN): "); 
                    scanf("%lf", &b->force_vector->data[2*node][0]);
                    b->force_vector->data[2*node][0] *= 1000;
                    load_count++;
                    break;
                case 2:
                    printf("\nMoment (kNm): "); 
                    scanf("%lf", &b->force_vector->data[2*node + 1][0]);
                    b->force_vector->data[2*node+1][0] *= 1000;
                    load_count++;
                    break;
                case 3:
                    b->fixed_vector->data[2*node][0] = 1; 
                    fix_count++;
                    break;
                case 4:
                    b->fixed_vector->data[2*node][0] = 1; 
                    b->fixed_vector->data[2*node + 1][0] = 1; 
                    fix_count += 2;
                    break;
            }
        }
    }
}

int main() {
    Beam *b = malloc(sizeof(Beam));

    char running = 'y';
    while (running == 'y') {

        // get properties and boudary conditions from the user
        get_beam_info(b);

        printf("\nConditions applied.");
        printf("\nCalculating displacements...\n");

        // calculate stiffness matrices
        get_lsm(b);
        get_gsm(b);     

        // remove fixpoints from equation
        for (int i = b->num_dof - 1; i >= 0; i--) {
            if (b->fixed_vector->data[i][0] == 1) {
                remove_row(b->GSM, i);
                remove_col(b->GSM, i);
                remove_row(b->force_vector, i);
                remove_row(b->disps_vector, i);
            }
        }

        // solve for displacements
        b->disps_vector = solve(b->GSM, b->force_vector);

        // re-insert fixpoints for zero displacements
        Matrix *zero = zeros(1, 1);
        for (int i = 0; i < b->num_dof; i++) {
            if (b->fixed_vector->data[i][0] == 1) { 
                insert_row(b->disps_vector, i, zero); 
                insert_row(b->force_vector, i, zero); 
            }
        }
        release(zero);
        
        show_constraint_table(b);

        // free memory
        release(b->force_vector); release(b->fixed_vector); release(b->disps_vector);
        release(b->LSM); release(b->GSM);

        // loop condition
        printf("\nAnalyse another beam? (y/n): "); scanf(" %c", &running);
    }
}