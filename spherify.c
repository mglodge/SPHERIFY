/* 

19/07/23 - Code that takes DDSCAT-style shapes, doubles their resolution, and then "rounds" the edges. See
Aerosol are note Spherical Cows, Lodge (2023) for details on the algorithm and visual examples.

INSTRUCTIONS FOR USE

- Name the input file "shape.dat" and place in the same folder as this code
- Additionally place "STAG_spherify.py" in the same folder if you wish to visualise the input and output files immediately
- Compile and run the code!


#
# The code works on the mac it was designed on, but has not been tested on any other computers - please e-mail
# the following e-mail address if you have any questions or issues: matthew.g.lodge@gmail.com
#
# TROUBLESHOOTING
#
# This exact version of STAG_spherify.py is designed to produce two windows (side-by-side) on mac computers; however, this will not
# work if running on windows. Simply comment out any lines that produce errors as a workaround, e.g:
#
# - matplotlib.use("TkAgg") # use a backend to allow the current_fig_manager section to position the windows
# - plt.get_current_fig_manager().window.wm_geometry("+750+100")
#
# and then uncomment the lines under "windows version".
# 
# Matt Lodge 01/07/22

*/
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <string.h>


FILE* DDSCAT_infile;
FILE* DDSCAT_outfile;
FILE* original_grid_outfile;
FILE* new_grid_outfile;
int diagnostics, i, j, k, x, y ,z, original_N, original_lattice_dim, new_lattice_dim, dipole_count, JA, IX, IY, IZ, ICOMPX, ICOMPY, ICOMPZ, STAG_lattice_dim,edgecase;
double STAG_odd_even_offset, min[3], max[3], STAG_offset[3];
char buf[1000];
int*** original_grid;
int*** new_grid;
double** dipole_info;
double** STAG_dipole_positions;



int main()
{


    printf("\n\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n\n                                                  WELCOME TO SPHERIFY!                            \n\n");


    printf("\n                  ________________                     __________                          .  =  .                    ");
    printf("\n                 |                |                   |          |                       '         '                  ");
    printf("\n                 |                |                  (            )                    /             \\                ");
    printf("\n                 |                |                 |              |                  |               |               ");
    printf("\n                 |                |      --->      (                )      --->       |               |               ");
    printf("\n                 |                |                 |              |                   \\             /                ");
    printf("\n                 |                |                  (            )                      .         .                  ");
    printf("\n                 |________________|                   |__________|                         '  =  '                    ");



    printf("\n\n\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n ---------------------------------------------------------------------------------------------------------------------\n\n");

    diagnostics=0; //set = 1 to print diagnostic statements

    /* read in shape.dat file */

    if((DDSCAT_infile=fopen("shape.dat","r")) == NULL){
        printf("\n\nError- the shape file cannot be found!! \n\n\n");
        system("pause");
        return 1;
    }
    else{
        DDSCAT_infile=fopen("shape.dat","r");
    }

    printf("\n Shape data file opened successfully. Analysing data:\n");

    /* read the information from the header of the DDSCAT file and find where the data starts (these header sections can be flexible and have different numbers of lines before the data starts -- so we just find where line containing JA, IX, IY, IZ is and assume the data begins underneath)*/

    while (fgets(buf,1000, DDSCAT_infile)!=NULL){ //if something is read

        if(strstr(buf,"NAT") != NULL){ //check if "NAT" is in the string for this line
            if(1==sscanf(buf,"%*[^0123456789]%d", &original_N)){ //record the number of dipoles (the first instance where any of the numbers 0->9 appear in this line)
                printf("\n \t NAT found: %d dipoles are in the original data file.", original_N);
            }
        }

        if(strstr(buf,"JA") != NULL){
            if(strstr(buf,"IX") != NULL){       //check if the current line contains the strings "JA", "IX", "IY" and "IZ" (seperately in case the spacing is changed by the user).
                if(strstr(buf,"IY") != NULL){
                    if(strstr(buf,"IZ") != NULL){
                        break;		//then exit the loop -- data is about to start
                    }
                }
            }
        }
    }

    /* initialise matrix to store dipoles, now that we know how many there were originally */
    dipole_info=(double**)malloc(original_N*sizeof(double*));
    for (i=0; i<original_N; i++) {
        dipole_info[i]=(double*)malloc((6)*sizeof(double));   //matrix will record each of the values for "IX IY IZ ICOMPX ICOMPY ICOMPZ" in array elements [0]->[5]
    }

    /* begin recording values and saving only the required composition (default ==1, for soot) */

    k=0; //counts the number of dipoles with composition 1

    while(fscanf(DDSCAT_infile, " %d  %d  %d  %d  %d  %d  %d", &JA, &IX, &IY, &IZ, &ICOMPX, &ICOMPY, &ICOMPZ)>0){ //scan each line of data after this point and record the values as individual buffers JA, IX etc...

        dipole_info[k][0]=IX;
        dipole_info[k][1]=IY;
        dipole_info[k][2]=IZ;            //save info for all dipoles
        dipole_info[k][3]=ICOMPX;
        dipole_info[k][4]=ICOMPY;
        dipole_info[k][5]=ICOMPZ;
        k++; //we found a dipole - keep track of total number
    }

    if(k==0){
        printf("\n\nError- no dipoles were found!! \n\n\n");
        system("pause");
        return 1;
    }
    else{
        printf("\n\n %d dipoles successfully imported.", k);
    }

    fclose(DDSCAT_infile);

    /* convert dipole positions to STAG grid format -- all need to be > 0 (positive integers)*/

    printf(" \n Translating %d dipoles to positive values... ", original_N);
    /* make new matrix to store STAG version of dipoles for 3D visualisation */
    STAG_dipole_positions=(double**)malloc(original_N*sizeof(double*));
    for (i=0; i<original_N; i++) {
        STAG_dipole_positions[i]=(double*)malloc((3)*sizeof(double));  //makes matrix to store dipole positions, now that we know the number of dipoles within the object geometry and the number of spheres
    }

    /* Search to find the most negative x,y,z points in the dipole positions - in a moment, we will need to translate them again to make them all positive for viewing in S.T.A.G */
    min[0]=0;
    min[1]=0;
    min[2]=0;
    max[0]=0;
    max[1]=0;
    max[2]=0;

    for(i=0;i<original_N;i++){
        if(dipole_info[i][0]<min[0]){
            min[0]=dipole_info[i][0]; //find min x
        }
        if(dipole_info[i][1]<min[1]){
            min[1]=dipole_info[i][1]; //find min y
        }
        if(dipole_info[i][2]<min[2]){
            min[2]=dipole_info[i][2]; //find min z
        }

        if(dipole_info[i][0]>max[0]){
            max[0]=dipole_info[i][0]; //find max x
        }
        if(dipole_info[i][1]>max[1]){
            max[1]=dipole_info[i][1]; //find max y
        }
        if(dipole_info[i][2]>max[2]){
            max[2]=dipole_info[i][2]; //find max z
        }
    }

    STAG_lattice_dim=0;
    if((max[0]-min[0])>STAG_lattice_dim){
        STAG_lattice_dim = ceil(max[0]-min[0]); //if the x-coord has the biggest variation, use this as our lattice size in S.T.A.G
    }
    if((max[1]-min[1])>STAG_lattice_dim){
        STAG_lattice_dim = ceil(max[1]-min[1]); //if the y-coord has the biggest variation, use this as our lattice size in S.T.A.G
    }
    if((max[2]-min[2])>STAG_lattice_dim){
        STAG_lattice_dim = ceil(max[2]-min[2]); //if the z-coord has the biggest variation, use this as our lattice size in S.T.A.G
    }
    STAG_lattice_dim= STAG_lattice_dim + 1; //then add 1 to make the lattice large enough to hold the final values
    original_lattice_dim=STAG_lattice_dim; //save this value for the original grid

    //printf("\n\n\n\n min_x= %f \n min_y= %f \n min_z= %f \n max_x= %f \n max_y= %f \n max_z= %f \n\n STAG_lattice_dim = %d \n", min[0], min[1], min[2], max[0], max[1], max[2], STAG_lattice_dim);

    /* calculate how much the non-"controller" axes will be moved by in STAG (this basically finds the "middle" square of the lattice, which depends if the lattice is odd or even. See comments below for examples */
    if((STAG_lattice_dim%2)==0){ //if STAG_lattice is even
        STAG_odd_even_offset= STAG_lattice_dim/2.0;    //e.g. If the lattice is 16 x 16 x 16, the "middle" square is 8 (either 7 or 8 are the middle numbers between squares 0,1,2,3,4,5,6,7,!!8!!,9,10,11,12,13,14,15, and it doesn't matter which we pick -- you can't center a coordinate on an even lattice. Here we pick "8" and shift it slightly up, but it's barely noticeable for large N). We calculate it using 16/2 = 8.
    }
    else{ //if STAG_lattice is odd
        STAG_odd_even_offset= (STAG_lattice_dim-1)/2.0;   //e.g. If the lattice is 15 x 15 x 15, the "middle" square is "7" (because it's the middle number between squares 0,1,2,3,4,5,6,!!7!!,8,9,10,11,12,13,14). We find this by calculating (15-1)/2 =7.
    }

    /* find which axes is the "controller" of the grid size - this the the axis that is most negative, and thus needs the biggest shift to make it positive (whichever holds the smallest minimum value) */
    // These all work in the same way -- once we have the "controller", we base the shift for that axis on it's min[] value, and then simply "center" the other two axes using the offset calculated above (which finds the distance to the "middle" square of the odd or even STAG lattice)

    if(min[1]<min[0]){
        if(min[2]<min[1]){ //z axes is the controller
            if(ceil(min[0])==ceil(min[2])){
                STAG_offset[0]= -1.0*min[0]; //if the x-axis has the same minimum as the controller, also translate it by the same amount.
            }
            else{
                STAG_offset[0]= STAG_odd_even_offset - ceil((max[0]+min[0])/2.0); //otherwise center it in the lattice
            }

            if(ceil(min[1])==ceil(min[2])){
                STAG_offset[1]= -1.0*min[1]; //if the y-axis has the same minimum as the controller, also translate it by the same amount.
            }
            else{
                STAG_offset[1]= STAG_odd_even_offset - ceil((max[1]+min[1])/2.0); //otherwise center it in the lattice
            }

            STAG_offset[2]= -1.0*min[2]; //move the z-axis to get it's smallest value increased up to zero

            //printf("\n The z-axis is the controller (switch 1).");
            //printf("\n\n STAG_offset[0] = %f \n STAG_offset[1] = %f \n STAG_offset[2] = %f (controller) \n ", STAG_offset[0], STAG_offset[1], STAG_offset[2]);
        }
        else{ //y axis is the controller

            if(ceil(min[0])==ceil(min[1])){
                STAG_offset[0]= -1.0*min[0]; //if the x-axis has the same minimum as the controller, also translate it by the same amount.
            }
            else{
                STAG_offset[0]= STAG_odd_even_offset - ceil((max[0]+min[0])/2.0); //otherwise center it in the lattice
            }

            STAG_offset[1]= -1.0*min[1]; //move the y-axis to get it's smallest value increased up to zero

            if(ceil(min[2])==ceil(min[1])){
                STAG_offset[2]= -1.0*min[2]; //if the z-axis has the same minimum as the controller, also translate it by the same amount.
            }
            else{
                STAG_offset[2]= STAG_odd_even_offset - ceil((max[2]+min[2])/2.0); //otherwise center it in the lattice
            }

            //printf("\n The y-axis is the controller.");
            //printf("\n\n STAG_offset[0] = %f \n STAG_offset[1] = %f (controller) \n STAG_offset[2] = %f \n ", STAG_offset[0], STAG_offset[1], STAG_offset[2]);
        }
    }
    else if (min[2]<min[0]){ //z axis is the controller
        if(ceil(min[0])==ceil(min[2])){
            STAG_offset[0]= -1.0*min[0]; //if the x-axis has the same minimum as the controller, also translate it by the same amount.
        }
        else{
            STAG_offset[0]= STAG_odd_even_offset - ceil((max[0]+min[0])/2.0); //otherwise center it in the lattice
        }

        if(ceil(min[1])==ceil(min[2])){
            STAG_offset[1]= -1.0*min[1]; //if the y-axis has the same minimum as the controller, also translate it by the same amount.
        }
        else{
            STAG_offset[1]= STAG_odd_even_offset - ceil((max[1]+min[1])/2.0); //otherwise center it in the lattice
        }

        STAG_offset[2]= -1.0*min[2]; //move the z-axis to get it's smallest value increased up to zero

        //printf("\n The z-axis is the controller (switch 2).");
        //printf("\n\n STAG_offset[0] = %f \n STAG_offset[1] = %f \n STAG_offset[2] = %f (controller) \n ", STAG_offset[0], STAG_offset[1], STAG_offset[2]);
    }
    else{ //x axis is the controller
        STAG_offset[0]= -1.0*min[0]; //move the x-axis to get it's smallest value increased up to zero

        if(ceil(min[1])==ceil(min[0])){
            STAG_offset[1]= -1.0*min[1]; //if the y-axis has the same minimum as the controller, also translate it by the same amount.
        }
        else{
            STAG_offset[1]= STAG_odd_even_offset - ceil((max[1]+min[1])/2.0); //otherwise center it in the lattice
        }

        if(ceil(min[2])==ceil(min[0])){
            STAG_offset[2]= -1.0*min[2]; //if the z-axis has the same minimum as the controller, also translate it by the same amount.
        }
        else{
            STAG_offset[2]= STAG_odd_even_offset - ceil((max[2]+min[2])/2.0); //otherwise center it in the lattice
        }

        //printf("\n The x-axes is the controller.");
        //printf("\n\n STAG_offset[0] = %f (controller) \n STAG_offset[1] = %f \n STAG_offset[2] = %f \n ", STAG_offset[0], STAG_offset[1], STAG_offset[2]);
    }

    //printf(" \n\n STAG_lattice_dim = %d \n STAG_odd_even_offset = %f     (this is the value that the 'non controller' axes will shift by in STAG coords. It is the 'middle square' of the STAG lattice (depends if STAG lattice is odd or even).\n\n", STAG_lattice_dim, STAG_odd_even_offset);

    /* Adjust the dipole positions by the offsets found above and same them as positive values for STAG */

    //printf("\n Dipole position \t\t\t STAG position\n ");
    for(i=0;i<original_N;i++){
        STAG_dipole_positions[i][0]= dipole_info[i][0] + STAG_offset[0];
        STAG_dipole_positions[i][1]= dipole_info[i][1] + STAG_offset[1];
        STAG_dipole_positions[i][2]= dipole_info[i][2] + STAG_offset[2];
    }

    /* print diagnostics and dipole transformations - two different versions, either for imported coords or sphere-made coords

    //printf("\n\n\n                 SPHERE POSITION                   LOCATIONS                      TRANSLATED POSITION          FINAL STAG POSITION \n\n");
    //for(i=0;i<number_of_spheres;i++){
    //    for(j=0;j<sphere_N;j++){
    //        //printf(" %4d %4d     %5.2f %5.2f %5.2f       +       %5.2f %5.2f %5.2f     -->           %5.2f %5.2f %5.2f            %5.2f %5.2f %5.2f \n", i, j, sphere_dipole_positions[j][0], sphere_dipole_positions[j][1], sphere_dipole_positions[j][2], sphere_locations[i][0]+sphere_odd_even_offset, sphere_locations[i][1]+sphere_odd_even_offset, sphere_locations[i][2]+sphere_odd_even_offset, dipole_positions[i*sphere_N+j][0], dipole_positions[i*sphere_N+j][1], dipole_positions[i*sphere_N+j][2], STAG_dipole_positions[i*sphere_N+j][0], STAG_dipole_positions[i*sphere_N+j][1], STAG_dipole_positions[i*sphere_N+j][2]);
    //    }
    //}

    printf("\n\n\n             IMPORTED POSITION      ( + %.1f %.1f %.1f )       FINAL STAG POSITION \n\n", STAG_offset[0], STAG_offset[1], STAG_offset[2]);
    for(i=0;i<original_N;i++){
        printf(" %4d      %6.2f %6.2f %6.2f             -->              %6.2f %6.2f %6.2f            \n", i, dipole_info[i][0], dipole_info[i][1], dipole_info[i][2], STAG_dipole_positions[i][0], STAG_dipole_positions[i][1], STAG_dipole_positions[i][2]);
    }


    /* initialise original grid array */

    original_grid = (int***)malloc(original_lattice_dim*sizeof(int**));
    for (i=0;i<original_lattice_dim;i++) {
        original_grid[i] = (int**)malloc(original_lattice_dim*sizeof(int*));
        for (j=0;j<original_lattice_dim;j++) {
          original_grid[i][j] = (int*)malloc(original_lattice_dim*sizeof(int));
        }
    }

    /* initially, set values at all positions to 0 */

    for(x=0;x<original_lattice_dim;x++){
        for(y=0;y<original_lattice_dim;y++){
            for(z=0;z<original_lattice_dim;z++){
                original_grid[x][y][z]=0; //set values at all positions to 0
            }
        }
    }

    /* then go through list of dipoles, saving their x-y-z coords, and adjust the value or in the original_grid array to 1 if there is a dipole at this position */

    for(i=0;i<original_N;i++){
        x=STAG_dipole_positions[i][0];
        y=STAG_dipole_positions[i][1];
        z=STAG_dipole_positions[i][2];

        //printf("\n x=%d y= %d z= %d",x,y,z);
        original_grid[x][y][z]=1; //set the value at this position to 1
    }

   /*for(x=0;x<original_lattice_dim;x++){
        for(y=0;y<original_lattice_dim;y++){
            for(z=0;z<original_lattice_dim;z++){
                printf("\n ALTERED GRID: %d", original_grid[x][y][z]);
            }
        }
    }*/

    printf(" Translation complete. (%d x %d x %d) grid created. \n\n", original_lattice_dim, original_lattice_dim, original_lattice_dim);

    /* save original resolution output to STAG_spherify data file (to visualise it in 3D) */
    printf(" Exporting original data to S.T.A.G...");

    original_grid_outfile=fopen("original.txt","w"); //open file for saving dipole positions

    dipole_count=0;
    for(x=0;x<original_lattice_dim;x++){
        for(y=0;y<original_lattice_dim;y++){
            for(z=0;z<original_lattice_dim;z++){
                if(original_grid[x][y][z]==1){
                    fprintf(original_grid_outfile,"%d, %d, %d\n", x,y,z); //save x-y-z coords of any dipoles that have values > 0. Any dipoles will have values == 1 at this stage.
                    dipole_count++; //keep track of how many dipoles we are recording
                }
            }
        }
    }

    /* The final row is extra data needed for the python visulisation S.T.A.G program, NOT a dipole! */
    fprintf(original_grid_outfile,"%d, %d, %d\n", original_lattice_dim, dipole_count, original_lattice_dim); //the final row contains the number of dipoles, the grid size, and a random number just to keep the shape of three columns for python to read */
    fclose(original_grid_outfile);

    printf(" Export complete.\n");

    /* initialise new 3D grid at higher resolution */

    new_lattice_dim= 2*original_lattice_dim; // new grid resolution will be twice as large

    new_grid = (int***)malloc(new_lattice_dim*sizeof(int**));
    for (i=0;i<new_lattice_dim;i++) {
        new_grid[i] = (int**)malloc(new_lattice_dim*sizeof(int*));
        for (j=0;j<new_lattice_dim;j++) {
          new_grid[i][j] = (int*)malloc(new_lattice_dim*sizeof(int));
        }
    }

    //set all values == 0

    for(x=0;x<new_lattice_dim;x++){
        for(y=0;y<new_lattice_dim;y++){
            for(z=0;z<new_lattice_dim;z++){
                new_grid[x][y][z]=0;
            }
        }
    }

    printf("\n New high-resolution grid initialised (%d x %d x %d).\n\n", new_lattice_dim, new_lattice_dim, new_lattice_dim);

    /* search through y-z slices along the x-axis */

    for(x=0;x<original_lattice_dim;x++){
        for(y=0;y<original_lattice_dim;y++){
            for(z=0;z<original_lattice_dim;z++){

                // edge definitions
                // original_grid[x][y][z+1] == right
                // original_grid[x][y][z-1] == left
                // original_grid[x][y+1][z] == top
                // original_grid[x][y-1][z] == bottom

                edgecase=0;

                /* check if there are only two occupied edges, and define which type of shape they are in */

                //edge case 1: bottom-left

                if((z>0)&&(y>0)&&((original_grid[x][y][z-1]==1)&&(original_grid[x][y-1][z]==1))&&((z==(original_lattice_dim-1))||(original_grid[x][y][z+1]==0))&&((y==(original_lattice_dim-1))||(original_grid[x][y+1][z]==0))){ // LOGIC: Only check this edge if z and y are > 0 (if z and y are not in the bottom row/column). Check if bottom and left edges are occupied (==1). Finally, seperately check that the other two edges are either the boundaries of the grid (in the far-right column or top row of the grid for this case) or unoccupied (==0).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 1 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x][2*y][2*z]=3; //follow rule 1 for edge case 1
                        new_grid[2*x+1][2*y][2*z]=3; // repeat for new x-position (because we double the resolution, there are two x-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((y==(original_lattice_dim-1))||(z==0)||(original_grid[x][y+1][z-1]==0))&&((y==0)||(z==(original_lattice_dim-1))||(original_grid[x][y-1][z+1]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]++; //inner edge
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2x+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x+1][2*y][2*z]++; //inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                }

                //edge case 2: top-left

                else if((z>0)&&(y<(original_lattice_dim-1))&&((original_grid[x][y+1][z]==1)&&(original_grid[x][y][z-1]==1))&&((y==0)||(original_grid[x][y-1][z]==0))&&((z==(original_lattice_dim-1))||(original_grid[x][y][z+1]==0))){ // check if left and top edges are occupied (and check that the other two edges are unoccupied). Only check this edge if z>0 and y<new_lattice_dim (check z is not in first column and y is not in the top row).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 2 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x][2*y+1][2*z]=3; //follow rule 1 for edge case 2
                        new_grid[2*x+1][2*y+1][2*z]=3; // repeat for new x-position (because we double the resolution, there are two x-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((y==(original_lattice_dim-1))||(z==(original_lattice_dim-1))||(original_grid[x][y+1][z+1]==0))&&((y==0)||(z==0)||(original_grid[x][y-1][z-1]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z]++; //inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2x+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]++; //inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //edge case 3: top-right

                else if((z<(original_lattice_dim-1))&&(y<(original_lattice_dim-1))&&((original_grid[x][y][z+1]==1)&&(original_grid[x][y+1][z]==1))&&((z==0)||(original_grid[x][y][z-1]==0))&&((y==0)||(original_grid[x][y-1][z]==0))){ // check if top and right edges are occupied (and check that the other two edges are unoccupied). Only check this edge if z and y < new_lattice_dim (check z and y are not in the top row/column).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 3 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x][2*y+1][2*z+1]=3; //follow rule 1 for edge case 3
                        new_grid[2*x+1][2*y+1][2*z+1]=3; // repeat for new x-position (because we double the resolution, there are two x-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((y==(original_lattice_dim-1))||(z==0)||(original_grid[x][y+1][z-1]==0))&&((y==0)||(z==(original_lattice_dim-1))||(original_grid[x][y-1][z+1]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]++; //inner edge
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2x+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]++; //inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //edge case 4: bottom-right

                else if((z<(original_lattice_dim-1))&&(y>0)&&((original_grid[x][y-1][z]==1)&&(original_grid[x][y][z+1]==1))&&((y==(original_lattice_dim-1))||(original_grid[x][y+1][z]==0))&&((z==0)||(original_grid[x][y][z-1]==0))){ // check if right and bottom edges are occupied (and check that the other two edges are unoccupied). Only check this edge if y>0 and z<new_lattice_dim (check y is not in first row z is not in the final column). Finally, this code "((y==(N-1))||(original_grid[x][y+1][z]==0))" checks if the edges are the boundaries of the grid (and so are also unuccupied), or whether the next cell is there but unoccupied.

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 4 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x][2*y][2*z+1]=3; //follow rule 1 for edge case 4
                        new_grid[2*x+1][2*y][2*z+1]=3; // repeat for new x-position (because we double the resolution, there are two x-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((y==(original_lattice_dim-1))||(z==(original_lattice_dim-1))||(original_grid[x][y+1][z+1]==0))&&((y==0)||(z==0)||(original_grid[x][y-1][z-1]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y][2*z+1]++; //inner edge

                        //repeat for 2x+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]++; //inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //if none of the edge cases are satisfied (or if an edge was found, but the checks for RULE 2 failed, leaving edgecase==0 still), fill the new grid with 1's if occupied
                if((edgecase==0)&&(original_grid[x][y][z]==1)){
                    if(diagnostics==1){
                        printf("\n %d %d %d   No edge case found, filling with 1's.", x, y, z);
                    }
                    new_grid[2*x][2*y][2*z]++;
                    new_grid[2*x][2*y+1][2*z]++;
                    new_grid[2*x][2*y+1][2*z+1]++;
                    new_grid[2*x][2*y][2*z+1]++;

                    new_grid[2*x+1][2*y][2*z]++;       // 1 dipole in original grid -> 8 dipoles in new grid
                    new_grid[2*x+1][2*y+1][2*z]++;
                    new_grid[2*x+1][2*y+1][2*z+1]++;
                    new_grid[2*x+1][2*y][2*z+1]++;
                }
                //or 0's if unoccupied
                else if(edgecase==0){
                    if(diagnostics==1){
                        printf("\n %d %d %d   No edge case found, filling with 0's.", x, y, z);
                    }
                    new_grid[2*x][2*y][2*z]--;
                    new_grid[2*x][2*y+1][2*z]--;
                    new_grid[2*x][2*y+1][2*z+1]--;
                    new_grid[2*x][2*y][2*z+1]--;

                    new_grid[2*x+1][2*y][2*z]--;       // 1 dipole in original grid -> 8 dipoles in new grid
                    new_grid[2*x+1][2*y+1][2*z]--;
                    new_grid[2*x+1][2*y+1][2*z+1]--;
                    new_grid[2*x+1][2*y][2*z+1]--;
                }

            }
        }
    }

    /* diagnostics - print original grid to screen (only for small grid sizes that will display in the console) */

    if(diagnostics==1){
        printf ("\n\n --------------------- Original grid ---------------------\n");
        for(x=0;x<original_lattice_dim;x++){
            printf("\n  y-z slice at x = %d    \n", x);
            printf("\n        Coordinates                                 Values                                       xyz coords  \n\n");
            for(y=(original_lattice_dim-1);y>(-1);y--){
                for(z=0;z<original_lattice_dim;z++){
                    printf("   %d-%d   ", y,z); //print coord positions
                }
                printf("                  ");
                for(z=0;z<original_lattice_dim;z++){
                    printf("   %d   ", original_grid[x][y][z]); //print values at these coords
                }
                printf("                  ");
                for(z=0;z<original_lattice_dim;z++){
                    printf("   %d %d %d   ", x,y,z); //print values at these coords
                }
                printf("\n\n");
            }
            printf("\n\n\n\n");
        }

        printf ("\n\n --------------------- New grid ---------------------\n");
        for(x=0;x<new_lattice_dim;x++){
            printf("\n  y-z slice at x = %d    \n", x);
            printf("\n                       Coordinates                                                       Values                                                            xyz coords \n\n");
            for(y=(new_lattice_dim-1);y>(-1);y--){
                for(z=0;z<new_lattice_dim;z++){
                    printf("   %d-%d   ", y,z); //print coord positions
                }
                printf("                  ");
                for(z=0;z<new_lattice_dim;z++){
                    printf("   %d   ", new_grid[x][y][z]); //print values at these coords
                }
                printf("                  ");
                for(z=0;z<new_lattice_dim;z++){
                    printf("   %d %d %d   ", x,y,z); //print values at these coords
                }
                printf("\n\n");
            }
            printf("\n\n\n\n");
        }


        /* save high resolution output to STAG_spherify data file */

        new_grid_outfile=fopen("high_res.txt","w"); //open file for saving dipole positions

        dipole_count=0;
        for(x=0;x<new_lattice_dim;x++){
            for(y=0;y<new_lattice_dim;y++){
                for(z=0;z<new_lattice_dim;z++){
                    if(new_grid[x][y][z]>0){
                        fprintf(new_grid_outfile,"%d, %d, %d\n", x,y,z); //save x-y-z coords of any dipoles that have values > 0 (should still be dipoles, on average, after three "sweeps" in the x-,  y- and z- directions)
                        dipole_count++; //keep track of how many dipoles we are recording
                    }
                }
            }
        }

        /* The final row is extra data needed for the python visulisation S.T.A.G program, NOT a dipole! */
        fprintf(new_grid_outfile,"%d, %d, %d\n", new_lattice_dim, dipole_count, new_lattice_dim); //the final row contains the number of dipoles, the grid size, and a random number just to keep the shape of three columns for python to read */
        fclose(new_grid_outfile);

        /* run xSTAG_spherify on each file - views both the original and higher resolution images in 3D */

        system("xSTAG_spherify.bat"); //opens a batch file with a command to run STAG_spherify as a python script
    }



    /* search through z-x slices along the y-axis */

    for(y=0;y<original_lattice_dim;y++){
        for(z=0;z<original_lattice_dim;z++){
            for(x=0;x<original_lattice_dim;x++){

                // edge definitions
                // original_grid[x+1][y][z] == right
                // original_grid[x-1][y][z] == left
                // original_grid[x][y][z+1] == top
                // original_grid[x][y][z-1] == bottom

                edgecase=0;

                /* check if there are only two occupied edges, and define which type of shape they are in */

                //edge case 1: bottom-left

                if((z>0)&&(x>0)&&((original_grid[x][y][z-1]==1)&&(original_grid[x-1][y][z]==1))&&((x==(original_lattice_dim-1))||(original_grid[x+1][y][z]==0))&&((z==(original_lattice_dim-1))||(original_grid[x][y][z+1]==0))){ // LOGIC: Only check this edge if z and x are > 0 (if z and x are not in the bottom row/column). Check if bottom and left edges are occupied (==1). Finally, seperately check that the other two edges are either the boundaries of the grid (in the far-right column or top row of the grid for this case) or unoccupied (==0).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 1 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x][2*y][2*z]=3; //follow rule 1 for edge case 1
                        new_grid[2*x][2*y+1][2*z]=3; // repeat for new y-position (because we double the resolution, there are two y-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((z==(original_lattice_dim-1))||(x==0)||(original_grid[x-1][y][z+1]==0))&&((z==0)||(x==(original_lattice_dim-1))||(original_grid[x+1][y][z-1]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]++; //inner edge
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2y+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x][2*y+1][2*z]++; //inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                }

                //edge case 2: top-left

                else if((x>0)&&(z<(original_lattice_dim-1))&&((original_grid[x][y][z+1]==1)&&(original_grid[x-1][y][z]==1))&&((z==0)||(original_grid[x][y][z-1]==0))&&((x==(original_lattice_dim-1))||(original_grid[x+1][y][z]==0))){ // check if left and top edges are occupied (and check that the other two edges are unoccupied). Only check this edge if x>0 and z<new_lattice_dim (check x is not in first column and z is not in the top row).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 2 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x][2*y][2*z+1]=3; //follow rule 1 for edge case 2
                        new_grid[2*x][2*y+1][2*z+1]=3; // repeat for new y-position (because we double the resolution, there are two y-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((z==(original_lattice_dim-1))||(x==(original_lattice_dim-1))||(original_grid[x+1][y][z+1]==0))&&((z==0)||(x==0)||(original_grid[x-1][y][z-1]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y][2*z+1]++; //inner edge

                        //repeat for 2y+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]++; //inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //edge case 3: top-right

                else if((x<(original_lattice_dim-1))&&(z<(original_lattice_dim-1))&&((original_grid[x][y][z+1]==1)&&(original_grid[x+1][y][z]==1))&&((x==0)||(original_grid[x-1][y][z]==0))&&((z==0)||(original_grid[x][y][z-1]==0))){ // check if top and right edges are occupied (and check that the other two edges are unoccupied). Only check this edge if z and x < new_lattice_dim (check z and x are not in the top row/column).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 3 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x+1][2*y][2*z+1]=3; //follow rule 1 for edge case 3
                        new_grid[2*x+1][2*y+1][2*z+1]=3; // repeat for new y-position (because we double the resolution, there are two y-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((z==(original_lattice_dim-1))||(x==0)||(original_grid[x-1][y][z+1]==0))&&((z==0)||(x==(original_lattice_dim-1))||(original_grid[x+1][y][z-1]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]++; //inner edge
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2y+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]++; //inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //edge case 4: bottom-right

                else if((x<(original_lattice_dim-1))&&(z>0)&&((original_grid[x][y][z-1]==1)&&(original_grid[x+1][y][z]==1))&&((z==(original_lattice_dim-1))||(original_grid[x][y][z+1]==0))&&((x==0)||(original_grid[x-1][y][z]==0))){ // check if right and bottom edges are occupied (and check that the other two edges are unoccupied). Only check this edge if y>0 and z<new_lattice_dim (check y is not in first row z is not in the final column). Finally, this code "((y==(N-1))||(original_grid[x][y+1][z]==0))" checks if the edges are the boundaries of the grid (and so are also unuccupied), or whether the next cell is there but unoccupied.

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 4 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x+1][2*y][2*z]=3; //follow rule 1 for edge case 4
                        new_grid[2*x+1][2*y+1][2*z]=3; // repeat for new y-position (because we double the resolution, there are two y-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((z==(original_lattice_dim-1))||(x==(original_lattice_dim-1))||(original_grid[x+1][y][z+1]==0))&&((z==0)||(x==0)||(original_grid[x-1][y][z-1]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z]++; //inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2y+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]++; //inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //if none of the edge cases are satisfied fill new grid with 1's if occupied
                if((edgecase==0)&&(original_grid[x][y][z]==1)){
                    if(diagnostics==1){
                        printf("\n %d %d %d   No edge case found, filling with 1's.", x, y, z);
                    }
                    new_grid[2*x][2*y][2*z]++;
                    new_grid[2*x+1][2*y][2*z]++;
                    new_grid[2*x+1][2*y][2*z+1]++;
                    new_grid[2*x][2*y][2*z+1]++;

                    new_grid[2*x][2*y+1][2*z]++;       // 1 dipole in original grid -> 8 dipoles in new grid
                    new_grid[2*x+1][2*y+1][2*z]++;
                    new_grid[2*x+1][2*y+1][2*z+1]++;
                    new_grid[2*x][2*y+1][2*z+1]++;
                }
                //or 0's if unoccupied
                else if(edgecase==0){
                    if(diagnostics==1){
                        printf("\n %d %d %d   No edge case found, filling with 0's.", x, y, z);
                    }
                    new_grid[2*x][2*y][2*z]--;
                    new_grid[2*x+1][2*y][2*z]--;
                    new_grid[2*x+1][2*y][2*z+1]--;
                    new_grid[2*x][2*y][2*z+1]--;

                    new_grid[2*x][2*y+1][2*z]--;       // 1 dipole in original grid -> 8 dipoles in new grid
                    new_grid[2*x+1][2*y+1][2*z]--;
                    new_grid[2*x+1][2*y+1][2*z+1]--;
                    new_grid[2*x][2*y+1][2*z+1]--;
                }

            }
        }
    }

    if(diagnostics==1){
        printf ("\n\n --------------------- Original grid ---------------------\n");
        for(y=0;y<original_lattice_dim;y++){
            printf("\n  z-x slice at y = %d    \n", y);
            printf("\n        Coordinates                                 Values                                       xyz coords  \n\n");
            for(z=(original_lattice_dim-1);z>(-1);z--){
                for(x=0;x<original_lattice_dim;x++){
                    printf("   %d-%d   ", z,x); //print coord positions
                }
                printf("                  ");
                for(x=0;x<original_lattice_dim;x++){
                    printf("   %d   ", original_grid[x][y][z]); //print values at these coords
                }
                printf("                  ");
                for(x=0;x<original_lattice_dim;x++){
                    printf("   %d %d %d   ", x,y,z); //print values at these coords
                }
                printf("\n\n");
            }
            printf("\n\n\n\n");
        }

        printf ("\n\n --------------------- New grid ---------------------\n");
        for(y=0;y<new_lattice_dim;y++){
            printf("\n  z-x slice at y = %d    \n", y);
            printf("\n                       Coordinates                                                       Values                                                            xyz coords \n\n");
            for(z=(new_lattice_dim-1);z>(-1);z--){
                for(x=0;x<new_lattice_dim;x++){
                    printf("   %d-%d   ", z,x); //print coord positions
                }
                printf("                  ");
                for(x=0;x<new_lattice_dim;x++){
                    printf("   %d   ", new_grid[x][y][z]); //print values at these coords
                }
                printf("                  ");
                for(x=0;x<new_lattice_dim;x++){
                    printf("   %d %d %d   ", x,y,z); //print values at these coords
                }
                printf("\n\n");
            }
            printf("\n\n\n\n");
        }

        /* save high resolution output to STAG_spherify data file */

        new_grid_outfile=fopen("high_res.txt","w"); //open file for saving dipole positions

        dipole_count=0;
        for(x=0;x<new_lattice_dim;x++){
            for(y=0;y<new_lattice_dim;y++){
                for(z=0;z<new_lattice_dim;z++){
                    if(new_grid[x][y][z]>0){
                        fprintf(new_grid_outfile,"%d, %d, %d\n", x,y,z); //save x-y-z coords of any dipoles that have values > 0 (should still be dipoles, on average, after three "sweeps" in the x-,  y- and z- directions)
                        dipole_count++; //keep track of how many dipoles we are recording
                    }
                }
            }
        }

        /* The final row is extra data needed for the python visulisation S.T.A.G program, NOT a dipole! */
        fprintf(new_grid_outfile,"%d, %d, %d\n", new_lattice_dim, dipole_count, new_lattice_dim); //the final row contains the number of dipoles, the grid size, and a random number just to keep the shape of three columns for python to read */
        fclose(new_grid_outfile);

        /* run xSTAG_spherify on each file - views both the original and higher resolution images in 3D */

        system("xSTAG_spherify.bat"); //opens a batch file with a command to run STAG_spherify as a python script
    }





    /* search through x-y slices along the z-axis */

    for(z=0;z<original_lattice_dim;z++){
        for(x=0;x<original_lattice_dim;x++){
            for(y=0;y<original_lattice_dim;y++){

                // edge definitions
                // original_grid[x][y+1][z] == right
                // original_grid[x][y-1][z] == left
                // original_grid[x+1][y][z] == top
                // original_grid[x-1][y][z] == bottom

                edgecase=0;

                /* check if there are only two occupied edges, and define which type of shape they are in */

                //edge case 1: bottom-left

                if((x>0)&&(y>0)&&((original_grid[x-1][y][z]==1)&&(original_grid[x][y-1][z]==1))&&((x==(original_lattice_dim-1))||(original_grid[x+1][y][1]==0))&&((y==(original_lattice_dim-1))||(original_grid[x][y+1][z]==0))){ // LOGIC: Only check this edge if x and y are > 0 (if z and y are not in the bottom row/column). Check if bottom and left edges are occupied (==1). Finally, seperately check that the other two edges are either the boundaries of the grid (in the far-right column or top row of the grid for this case) or unoccupied (==0).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 1 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x][2*y][2*z]=3; //follow rule 1 for edge case 1
                        new_grid[2*x][2*y][2*z+1]=3; // repeat for new z-position (because we double the resolution, there are two z-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((x==(original_lattice_dim-1))||(y==0)||(original_grid[x+1][y-1][z]==0))&&((x==0)||(y==(original_lattice_dim-1))||(original_grid[x-1][y+1][z]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]++; //inner edge
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2z+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x][2*y][2*z+1]++; //inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                }

                //edge case 2: top-left

                else if((y>0)&&(x<(original_lattice_dim-1))&&((original_grid[x+1][y][z]==1)&&(original_grid[x][y-1][z]==1))&&((x==0)||(original_grid[x-1][y][z]==0))&&((y==(original_lattice_dim-1))||(original_grid[x][y+1][z]==0))){ // check if left and top edges are occupied (and check that the other two edges are unoccupied). Only check this edge if z>0 and y<new_lattice_dim (check z is not in first column and y is not in the top row).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 2 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x+1][2*y][2*z]=3; //follow rule 1 for edge case 2
                        new_grid[2*x+1][2*y][2*z+1]=3; // repeat for new z-position (because we double the resolution, there are two z-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((x==(original_lattice_dim-1))||(y==(original_lattice_dim-1))||(original_grid[x+1][y+1][z]==0))&&((x==0)||(y==0)||(original_grid[x-1][y-1][z]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z]++; //inner edge

                        //repeat for 2z+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]++; //inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //edge case 3: top-right

                else if((x<(original_lattice_dim-1))&&(y<(original_lattice_dim-1))&&((original_grid[x+1][y][z]==1)&&(original_grid[x][y+1][z]==1))&&((x==0)||(original_grid[x-1][y][z]==0))&&((y==0)||(original_grid[x][y-1][z]==0))){ // check if top and right edges are occupied (and check that the other two edges are unoccupied). Only check this edge if z and y < new_lattice_dim (check z and y are not in the top row/column).

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 3 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x+1][2*y+1][2*z]=3; //follow rule 1 for edge case 3
                        new_grid[2*x+1][2*y+1][2*z+1]=3; // repeat for new z-position (because we double the resolution, there are two z-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((x==(original_lattice_dim-1))||(y==0)||(original_grid[x+1][y-1][z]==0))&&((x==0)||(y==(original_lattice_dim-1))||(original_grid[x-1][y+1][z]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z]++; //inner edge
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2z+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]++; //inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //edge case 4: bottom-right

                else if((y<(original_lattice_dim-1))&&(x>0)&&((original_grid[x-1][y][z]==1)&&(original_grid[x][y+1][z]==1))&&((x==(original_lattice_dim-1))||(original_grid[x+1][y][z]==0))&&((y==0)||(original_grid[x][y-1][z]==0))){ // check if right and bottom edges are occupied (and check that the other two edges are unoccupied). Only check this edge if y>0 and z<new_lattice_dim (check y is not in first row z is not in the final column). Finally, this code "((y==(N-1))||(original_grid[x][y+1][z]==0))" checks if the edges are the boundaries of the grid (and so are also unuccupied), or whether the next cell is there but unoccupied.

                    if(diagnostics==1){
                        printf("\n %d %d %d   Edge case 4 found.", x, y, z);
                    }
                    //rule 1
                    if(original_grid[x][y][z]==0){ //if cell is unoccupied
                        new_grid[2*x][2*y+1][2*z]=3; //follow rule 1 for edge case 4
                        new_grid[2*x][2*y+1][2*z+1]=3; // repeat for new z-position (because we double the resolution, there are two z-positions to apply this rule to!)
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }

                    //rule 2
                    if((original_grid[x][y][z]==1)&&((x==(original_lattice_dim-1))||(y==(original_lattice_dim-1))||(original_grid[x+1][y+1][z]==0))&&((x==0)||(y==0)||(original_grid[x-1][y-1][z]==0))){ //if cell is occupied, and both cells along the diagonal are also unoccupied (if either is occupied, we ignore rule 2 -- see 08/08/22). The code snippets like: "y==(N-1))||(z==0)||" are to check if we are at the edge of the grid -- if so, we don't want to check the next cell because we will be going out of bounds of our matrix, so we exit the statement. But as a whole, "((y==(N-1))||(z==0)||(original_grid[x][y+1][z-1]==0))" just checks "Is the North-West cell empty?"
                        new_grid[2*x][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z]++; //inner edge
                        new_grid[2*x+1][2*y+1][2*z]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z]--; //follow rule 2 and set all new cells to zero except inner edge

                        //repeat for 2z+1 positions (doubled-resolution in x-y-z causes 1 cube -> 8 cubes)
                        new_grid[2*x][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x][2*y+1][2*z+1]++; //inner edge
                        new_grid[2*x+1][2*y+1][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        new_grid[2*x+1][2*y][2*z+1]--; //follow rule 2 and set all new cells to zero except inner edge
                        edgecase=1; //set switch to state that an "edge case has been assigned"
                    }
                }

                //if none of the edge cases are satisfied fill new grid with 1's if occupied
                if((edgecase==0)&&(original_grid[x][y][z]==1)){
                    if(diagnostics==1){
                        printf("\n %d %d %d   No edge case found, filling with 1's.", x, y, z);
                    }
                    new_grid[2*x][2*y][2*z]++;
                    new_grid[2*x][2*y+1][2*z]++;
                    new_grid[2*x+1][2*y+1][2*z]++;
                    new_grid[2*x+1][2*y][2*z]++;

                    new_grid[2*x][2*y][2*z+1]++;       // 1 dipole in original grid -> 8 dipoles in new grid
                    new_grid[2*x][2*y+1][2*z+1]++;
                    new_grid[2*x+1][2*y+1][2*z+1]++;
                    new_grid[2*x+1][2*y][2*z+1]++;
                }
                //or 0's if unoccupied
                else if(edgecase==0){
                    if(diagnostics==1){
                        printf("\n %d %d %d   No edge case found, filling with 0's.", x, y, z);
                    }
                    new_grid[2*x][2*y][2*z]--;
                    new_grid[2*x][2*y+1][2*z]--;
                    new_grid[2*x+1][2*y+1][2*z]--;
                    new_grid[2*x+1][2*y][2*z]--;

                    new_grid[2*x][2*y][2*z+1]--;       // 1 dipole in original grid -> 8 dipoles in new grid
                    new_grid[2*x][2*y+1][2*z+1]--;
                    new_grid[2*x+1][2*y+1][2*z+1]--;
                    new_grid[2*x+1][2*y][2*z+1]--;
                }

            }
        }
    }

    if(diagnostics==1){
        printf ("\n\n --------------------- Original grid ---------------------\n");
        for(z=0;z<original_lattice_dim;z++){
            printf("\n  x-y slice at z = %d    \n", z);
            printf("\n        Coordinates                                 Values                                       xyz coords  \n\n");
            for(x=(original_lattice_dim-1);x>(-1);x--){
                for(y=0;y<original_lattice_dim;y++){
                    printf("   %d-%d   ", x,y); //print coord positions
                }
                printf("                  ");
                for(y=0;y<original_lattice_dim;y++){
                    printf("   %d   ", original_grid[x][y][z]); //print values at these coords
                }
                printf("                  ");
                for(y=0;y<original_lattice_dim;y++){
                    printf("   %d %d %d   ", x,y,z); //print values at these coords
                }
                printf("\n\n");
            }
            printf("\n\n\n\n");
        }

        printf ("\n\n --------------------- New grid ---------------------\n");
        for(z=0;z<new_lattice_dim;z++){
            printf("\n  x-y slice at z = %d    \n", z);
            printf("\n                       Coordinates                                                       Values                                                            xyz coords \n\n");
            for(x=(new_lattice_dim-1);x>(-1);x--){
                for(y=0;y<new_lattice_dim;y++){
                    printf("   %d-%d   ", x,y); //print coord positions
                }
                printf("                  ");
                for(y=0;y<new_lattice_dim;y++){
                    printf("   %d   ", new_grid[x][y][z]); //print values at these coords
                }
                printf("                  ");
                for(y=0;y<new_lattice_dim;y++){
                    printf("   %d %d %d   ", x,y,z); //print values at these coords
                }
                printf("\n\n");
            }
            printf("\n\n\n\n");
        }
    }

    /* save high resolution output to STAG_spherify data file */

    printf(" Analysis complete. Exporting high-resolution data to S.T.A.G...");

    new_grid_outfile=fopen("high_res.txt","w"); //open file for saving dipole positions

    dipole_count=0;
    for(x=0;x<new_lattice_dim;x++){
        for(y=0;y<new_lattice_dim;y++){
            for(z=0;z<new_lattice_dim;z++){
                if(new_grid[x][y][z]>0){
                    fprintf(new_grid_outfile,"%d, %d, %d\n", x,y,z); //save x-y-z coords of any dipoles that have values > 0 (should still be dipoles, on average, after three "sweeps" in the x-,  y- and z- directions)
                    dipole_count++; //keep track of how many dipoles we are recording
                }
            }
        }
    }

    /* The final row is extra data needed for the python visulisation S.T.A.G program, NOT a dipole! */
    fprintf(new_grid_outfile,"%d, %d, %d\n", new_lattice_dim, dipole_count, new_lattice_dim); //the final row contains the number of dipoles, the grid size, and a random number just to keep the shape of three columns for python to read */
    fclose(new_grid_outfile);

    printf(" Exported %d dipoles.\n", dipole_count);

    /* SAVE DATA IN DDSCAT FORMAT AT HIGHER RESOLUTION */

    printf("\n Re-centering and exporting high-resolution data in DDSCAT format.\n");

    /* we have all the info we need from the first scan -- write the new file */

    DDSCAT_infile=fopen("shape.dat","r");
    DDSCAT_outfile=fopen("shape2.dat","w");

    /* duplicate the information from the header of the DDSCAT file */

    printf("\n \t Duplicating header section... ");

    while (fgets(buf,1000, DDSCAT_infile)!=NULL){

        if(strstr(buf,"NAT") != NULL){ //check if "NAT" is in the string for this line
            fprintf(DDSCAT_outfile,"   %d   = NAT\n", dipole_count); // is so, print a special line for the new dipole number (at higher resolution)

        }
        else if(strstr(buf,"JA") != NULL){
            if(strstr(buf,"IX") != NULL){       //check if the current line contains the strings "JA", "IX", "IY" and "IZ" (seperately in case the spacing is changed by the user).
                if(strstr(buf,"IY") != NULL){
                    if(strstr(buf,"IZ") != NULL){
                        fputs(buf, DDSCAT_outfile); //print this final line and then...
                        break;		// ...we have finished the header section -- exit the loop to start writing data
                    }
                }
            }
        }
        else{
            fputs(buf, DDSCAT_outfile); //otherwise, copy and paste the line exactly "as is" from the old file to the new one
        }
    }

    /* save high-res dipole data */
    printf("\n \t Saving high-resolution dipole data for %d dipoles... ", dipole_count);

    k=0;
    for(x=0;x<new_lattice_dim;x++){
        for(y=0;y<new_lattice_dim;y++){
            for(z=0;z<new_lattice_dim;z++){
                if(new_grid[x][y][z]>0){

                    k++; //keep track of how many dipoles we are recording

                    fprintf(DDSCAT_outfile,"%10d %10.0f %10.0f %10.0f %10d %10d %10d\n", k, x-2*STAG_offset[0], y-2*STAG_offset[1], z-2*STAG_offset[2], ICOMPX, ICOMPY, ICOMPZ); //reverse the offset (doubled, because the grid size is doubled) and save the dipoles in their original "centred" positions, but with the high res interpolations and at twice the resolution
                }
            }
        }
    }


    printf("Done! \n\n Exported data for %d dipoles. Spherify program compete! Enjoy your new smooth shapes.\n\n", k);

    fclose(DDSCAT_infile);
    fclose(DDSCAT_outfile);

    /* run xSTAG_spherify on each file - views both the original and higher resolution images in 3D */

    //system("xSTAG_spherify.bat"); //WINDOWS VERSION: opens a batch file with a command to run STAG_spherify as a python script
    system("python STAG_spherify.py"); // MAC version -- open file in python

    /* free memory for arrays */

    for (i=0; i<original_N; i++) {
        free((void*)STAG_dipole_positions[i]);
    }
    free((void*)STAG_dipole_positions);

    for (i=0; i<original_N; i++) {
        free((void*)dipole_info[i]);
    }
    free((void*)dipole_info);

    for(i=0;i<original_lattice_dim;i++){
        for(j=0;j<original_lattice_dim;j++){
            free((void*)original_grid[i][j]);
        }
        free((void*)original_grid[i]);
    }
    free((void*)original_grid);

    for(i=0;i<new_lattice_dim;i++){
        for(j=0;j<new_lattice_dim;j++){
            free((void*)new_grid[i][j]);
        }
        free((void*)new_grid[i]);
    }
    free((void*)new_grid);

    return 0;
}
