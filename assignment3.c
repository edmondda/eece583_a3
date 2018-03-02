#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <limits.h>

#include "graphics.h"

//static void delay (void);
//static void button_press (float x, float y);
//static void drawscreen (void);



struct exampleData {
   int   num_cell;
   int   num_con;
   int   num_rows;
   int   num_cols;
   int * num_terminals;
   int ** net_connect;
};

struct partitionStruct {
   //float   initTemp;
   //float   currTemp;
   //bool    betaDecay;
   //bool    tempBetaSwitch;
   //float   tempBeta;
   //float   tempBeta2;
   //int     tempIter ;
   int     maxSteps;
   //bool    stepAccept;
   //float   accRate;
   //int     accCnt;
   //int     bb_size;
   //int     bb_x_max;
   //int     bb_x_min;
   //int     bb_y_max;
   //int     bb_y_min;
   //float   stopCostStd;
   //int     numStopIter;
};

struct location {
   int   x;
   int   y;
};

struct partition_cell {
   int   partition;
   bool  locked;
   float   gain;
   struct location location;
};



void readfile (const char* filename, struct exampleData * data);
void printPartition ( struct partition_cell * partitionCfg, struct exampleData * cfgStruct);
int calcCutSet ( struct partition_cell * partitionCfg, struct exampleData * cfgStruct);
int KerrighanLinStep(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct, int cutSet);
void KerrighanLinSwap(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct);
void initSchedule(struct exampleData * cfgStruct, struct partitionStruct * schedStruct);
float calcStd(int * x, int num_samples);
void unlockNodes( struct partition_cell * partitionCfg, struct exampleData * cfgStruct);
void calcGains(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct);

/***** GUI functions **************/

char       gRunMode        = 's'     ; /* run mode, default is step     */
char       gFooterLabel[1024]  =""; /* global footer text label   */
char       gHeaderLabel[1024]  =""; /* global header text label   */
char       gFooterMessage[1024]=""; /* global footer text message */
float      gWorldX = 1000         ; /* world X dimension          */
float      gWorldY = 1000         ; /* world X dimension          */

  struct partition_cell * partitionCfg;
  struct exampleData * cfgStruct;

inline void runStep  (void (*drawScreen_ptr)(void)) {gRunMode='s' ;} /* run one step */
inline void runPass  (void (*drawScreen_ptr)(void)) {gRunMode='p' ;} /* run one pass */
inline void runAll   (void (*drawScreen_ptr)(void)) {gRunMode='a' ;} /* run all      */


/* redrawing routine for still pictures. Redraw if user changes the window           */
void fpDraw(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, float xWorld, float yWorld);
void drawScreen() { clearscreen(); fpDraw(partitionCfg, cfgStruct,gWorldX,gWorldY); } /* clear & redraw  */

/* called whenever event_loop gets a button press in the graphics area.              */
void buttonPress  (float x, float y) {                                     }

/* receives the current mouse position in the current world as in init_world         */
void mouseMove    (float x, float y) {                                                }

/* function to handle keyboard press event, the ASCII character is returned          */
void keyPress     (int i) {                                                           }

/* show global message, wait until 'Proceed' pressed then draw screen                */
void waitLoop (struct partition_cell * partitionCfg, struct exampleData * cfgStruct)  {
  update_message(gFooterMessage);
  drawScreen(partitionCfg, cfgStruct);
  event_loop(buttonPress,drawScreen);
}

/*******************/


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


int main(int argc, char **argv) {

  printf("You have entered file = %s\n",  argv[1]);
  const char* dir = "benchmarks/";
  const char* filename = argv[1];
//  struct exampleData cfgStruct;
  struct exampleData * cfgStructPtr;

  struct partitionStruct   schedCfg;
  struct partitionStruct * schedStructPtr;

  bool printPart = false;

  if (argc > 2) {
    printPart = true;
  }

  int i, j, k;
  int size_connect, net_cons;
  int cutSet, cutSet_new, delta_cutSet;
  int * cost_list;
  float cost_std = 0.0;
  float acc_rate;
  int numStopIter = 0;

  cfgStructPtr = cfgStruct;

  char* read_file = (char*)malloc(strlen(dir)+strlen(filename));
  strcpy(read_file, dir);
  strcat(read_file, filename);
  readfile(read_file, cfgStructPtr);

  printf("num_cells = %d \n" , cfgStruct->num_cell);
  printf("num_con = %d \n"   , cfgStruct->num_con);
  printf("num_rows = %d \n"  , cfgStruct->num_rows);
  printf("num_cols = %d \n"  , cfgStruct->num_cols);
  printf("num_nets = %d \n"  , cfgStruct->num_con);
  for(i = 0; i < cfgStruct->num_con; i++ ) {
     printf("net # %d :", i);
     for(j = 0; j < cfgStruct->num_terminals[i]; j++) {
       printf(" %d", cfgStruct->net_connect[i][j] );
     }
     printf("\n");
  }

   // check that the number of cells will fit within the design
  if (cfgStruct->num_cell > cfgStruct->num_cols*cfgStruct->num_rows) {
    printf("ERROR : %dx%d grid is too small to place %d cells\n", cfgStruct->num_cols,
                                                               cfgStruct->num_rows,
                                                               cfgStruct->num_cell );
    return(0);
  }

  bool gGUI = true;
  if (gGUI) {

    // initialize display
    init_graphics((char*)"Kernighan-Lin bi-partitioning algorithm");
    init_world (0.,0.,gWorldX,gWorldY);

    /* Create new buttons */
    create_button ((char*)"Window"    , (char*)"---1"      , NULL     ); /* Separator    */
    // DE _ FIXME
    //create_button ((char*)"---1"      , (char*)"Run Step"  , runStep  ); /* run one step */
    //create_button ((char*)"Run Step"  , (char*)"Run Pass"  , runPass  ); /* run one pass */
    //create_button ((char*)"Run Pass"  , (char*)"Run All"   , runAll   ); /* run all      */
  }


  // DE - output file for results
  FILE *fp;
  char* results_file = (char*)malloc(strlen(filename)+strlen(".out"));
  strcpy(results_file, filename);
  strcat(results_file, ".out");
  //printf( results_file );
  fp = fopen(results_file, "w+");

  // DE - partition cells randoming into two groups
  //      unevenness of 1 is permitted


  // is it better to have two list (one A, one B)?
  // - we could just go through A
  //    - go through each net connection list
  //
  // list of cell with A/B flag (similiar to previous implementation)
  //  - Go through the net lists... if all entries are from the same partition, add cost_new
  //  - otherwise

//  struct partition_cell * partitionCfg;

  // place cells randomly in partition A and B to start
  partitionCfg = (struct partition_cell *)malloc (sizeof (struct partition_cell) * cfgStruct->num_cell);
  for(i = 0; i < cfgStruct->num_cell; i++) {
    partitionCfg[i].partition = i % 2;
    partitionCfg[i].locked = false;
    partitionCfg[i].gain = 0.0;
  }

  //DE - show initial 'random' placement

  if( printPart ) {
    printf("Initial partition: \n");
    printPartition ( partitionCfg, cfgStructPtr);
  }

  if (gGUI) {
  /* update global message and wait for 'Proceed' to proceed */
  sprintf(gHeaderLabel,"Design: %s  |  Cells: %u  |  Nets: %u : %s",
          filename, cfgStructPtr->num_cell, cfgStructPtr->num_con);
  if (gRunMode=='f') { /* show initial floorplan for 'final' mode */
    sprintf(gFooterMessage,"Initial partitioning. Press 'Proceed' to find final solution");
    sprintf(gFooterLabel,"Initial partitioning");
    waitLoop(partitionCfg, cfgStructPtr);
    //if (gPostScript) postscript(drawScreen);
  }
}


  // Calculate initial cut-set size
  cutSet = calcCutSet ( partitionCfg, cfgStructPtr );
  printf("initial cut set is : %d\n", cutSet);

  // K-L loop

  // initialize run parameters
  schedStructPtr = &schedCfg;
  initSchedule(cfgStructPtr, schedStructPtr);

  // output csv header
  // fputs("temperature, cost, acc_rate, bb_size\n", fp);

  for(i = 0; i < schedStructPtr->maxSteps; i++ ) {

    // Perform iteration of K-L algorithm
    cutSet = KerrighanLinStep(partitionCfg, cfgStructPtr, schedStructPtr, cutSet);
    //cutSet = calcCutSet ( partitionCfg, cfgStructPtr );

    // report statistics of run every 20 iterations
    printf("iteration %d, cut set = %d\n", i, cutSet );
  }

  // close output results file
  //fclose(fp);

  if (printPart) {
    printf("End partition :\n");
    printPartition( partitionCfg, cfgStructPtr);
  }
  cutSet = calcCutSet( partitionCfg, cfgStructPtr );

  printf("end cut set is : %d\n", cutSet);

  /* finished! wait still until 'Exit" is pressed */
  if (gGUI)
    while(1) waitLoop(partitionCfg, cfgStructPtr);

  return (0);
}


void readfile (const char* filename, struct exampleData * data) {

 FILE *fp;
 struct exampleData dataRead;
 int num_con_cells, cell_i;
 int test;
 int i,j;

 printf("opening file : %s\n", filename);

 fp = fopen(filename, "r");

 fscanf(fp, "%d %d %d %d\n",  &dataRead.num_cell,
                              &dataRead.num_con,
                              &dataRead.num_rows,
                              &dataRead.num_cols);

 dataRead.net_connect = (int**)malloc (sizeof (int*) * dataRead.num_con);
 dataRead.num_terminals = (int*)malloc (sizeof (int*) * dataRead.num_con);
 for(i = 0; i < dataRead.num_con; i++) {
   fscanf(fp, "%d", &num_con_cells);
   dataRead.num_terminals[i] = num_con_cells;
   dataRead.net_connect[i] = (int*)malloc (sizeof (int) * num_con_cells);
   for(j = 0; j < num_con_cells; j++ ) {
       fscanf(fp, " %d",  &dataRead.net_connect[i][j]);
   }
   fscanf(fp, "\n");
 }
 fclose(fp);

 *data = dataRead;


}



// Prints the cells in each partition
void printPartition ( struct partition_cell * partitionCfg, struct exampleData * cfgStruct) {

  int i,j;
  int num_digits;
  int * listA;
  int asize = 0;
  int bsize = 0;
  int * listB;

  // Determine partitioning

  listA = (int *)malloc(cfgStruct->num_cell * sizeof(int));
  listB = (int *)malloc(cfgStruct->num_cell * sizeof(int));

  // DE - FIXME - STORE IN STRUCT SO THAT WE DON'T HAVE TO RECALCULATE
  for(i=0; i < cfgStruct->num_cell; i++) {
    if ( partitionCfg[i].partition == 0 ) {
       listA[asize] = i;
       asize++;
    }
    else {
      listB[bsize] = i;
      bsize++;
    }
  }

  num_digits = (int)log10( (double)cfgStruct->num_cell) + 1;

  printf("LIST  A :\n" );
  for(i = 0; i < asize; i++ ) {
    printf(" %*d", num_digits, listA[i] );
  }
  printf("\n");

  printf("LIST  B :\n" );
  for(i = 0; i < bsize; i++ ) {
    printf(" %*d", num_digits, listB[i] );
  }
  printf("\n");

  free(listA);
  free(listB);
}


int calcCutSet ( struct partition_cell * partitionCfg, struct exampleData * cfgStruct) {

  int i,j;
  int cutSet;
  int cell;
  int partition_temp;

  // for reach net, determine the bounding box of the nells in the placement
  // grid
  cutSet = 0;

  for(i = 0; i < cfgStruct->num_con; i++ ) {

    cell = cfgStruct->net_connect[i][0];
    partition_temp = partitionCfg[ cell ].partition;

    for(j = 1; j < cfgStruct->num_terminals[i]; j++) {
      cell = cfgStruct->net_connect[i][j];
      if(partitionCfg[ cell ].partition != partition_temp) {
        cutSet++;
        break;
      }
    }
  }

  return cutSet;
}


void calcGains(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct) {


  int i, j;
  int cell;
  int partition_temp;

  int * listA;
  int * listB;
  int asize = 0;
  int bsize = 0;

  for(i = 0; i < cfgStruct->num_cell; i++ ) {
    partitionCfg[i].gain= 0.0;
  }

  // Iterate through all nets
  for(i = 0; i < cfgStruct->num_con; i++ ) {

    // Determine multi-terminal partitioning for each net

    //listA = (int *)malloc( sizeof(int));
    //listB = (int *)malloc( sizeof(int));
    asize = 0;
    bsize = 0;

    //listA = (int *)malloc(cfgStruct->num_cell * sizeof(int));
    //listB = (int *)malloc(cfgStruct->num_cell * sizeof(int));

    // DE - FIXME - STORE IN STRUCT SO THAT WE DON'T HAVE TO RECALCULATE
    for(j = 0; j < cfgStruct->num_terminals[i]; j++) {
      cell = cfgStruct->net_connect[i][j];
      if ( partitionCfg[cell].partition == 0 ) {
         //listA[asize] = j;
         asize++;
      }
      else {
        //listB[bsize] = j;
        bsize++;
      }
    }

    // Start with the original K-L gain equation
    float attract_alpha = 1.0;
    float repulse_alpha = 1.0;
    float scalar;
    // DE - let's try to improve this result be weighting smaller cliques  more
    if (asize + bsize > 4) {
      scalar = pow( (float)(asize + bsize - 4), 0.5);
    } else {
      scalar = 1.0;
    }
      scalar = 1.0;

    // with the partition information for each net, iterate again to get calcGains
    for(j = 0; j < cfgStruct->num_terminals[i]; j++) {

      cell = cfgStruct->net_connect[i][j];

      // determine if there is a single cross
      // DE - FIXME - somehow when this is wrong, we get better results? What has this enabled in the algo?
      if ( partitionCfg[ cell ].partition == 0 && asize == 1 ||
           partitionCfg[ cell ].partition == 1 && bsize == 1    ){
         // increment the gain of the cell
         partitionCfg[ cell ].gain += 1.0;
      }
      //  determine if a there will be a cross if moved
      else if ( partitionCfg[ cell ].partition == 0 && asize == cfgStruct->num_terminals[ i ] ||
                partitionCfg[ cell ].partition == 1 && bsize == cfgStruct->num_terminals[ i ]    ){
         // decrement the gains of the cell
         partitionCfg[ cell ].gain -= 1.0;
      }
      // Additional gain equation for multi-terminal nets
      else {

        if ( partitionCfg[ cell ].partition == 0 ) { // partition a

          partitionCfg[ cell ].gain += scalar*( attract_alpha*( (float) bsize / (float) (asize + bsize -1) ) - (repulse_alpha)* ( (float) (asize-1) / (float) (asize + bsize -1 ) ) );
          //partitionCfg[ cell ].gain += scalar*( attract_alpha*pow( ( (float) bsize / (float) (asize + bsize -1) ), 2) - (repulse_alpha)*pow( ( (float) (asize-1) / (float) (asize + bsize -1 ) ), 2) );

          //if (asize > bsize) {
          //  partitionCfg[ cell ].gain += 0.5 - pow(asize+1, -2)*(1 - bsize / (float) (asize))/2;
          //} else {
          //  partitionCfg[ cell ].gain += 0.5 + pow(asize+1, -2)*(1 - asize / (float) (bsize))/2;
          //}

          //partitionCfg[ cell ].gain += 1 - pow((( (float) bsize / (float) (asize + bsize - 1) ) - ( (float) (asize-1) / (float) (asize + bsize - 1) )), 2 );
          // if (asize > bsize ) {
          //   partitionCfg[ cell ].gain -= ( 1.0 / (float)asize );
          // }
          // else {
          //   partitionCfg[ cell ].gain += ( 1.0 / (float)bsize );
          // }

        }
        else { // partition b

          partitionCfg[ cell ].gain += scalar*( attract_alpha*( (float) asize / (float) (asize + bsize -1) ) - (repulse_alpha)*( (float) (bsize-1) / (float) (asize + bsize -1 ) ) );
          //partitionCfg[ cell ].gain += scalar*( attract_alpha*pow( ( (float) asize / (float) (asize + bsize -1) ), 2) - (repulse_alpha)*pow( ( (float) (bsize-1) / (float) (asize + bsize -1 ) ), 2) );

          //if (bsize > asize) {
          //  partitionCfg[ cell ].gain -= 1 - pow((( (float) asize / (float) (asize + bsize - 1) ) - ( (float) (bsize-1) / (float) (asize + bsize - 1) )), 2 );
          //} else {
          //  partitionCfg[ cell ].gain += 1 - pow((( (float) asize / (float) (asize + bsize - 1) ) - ( (float) (bsize-1) / (float) (asize + bsize - 1) )), 2 );
          //}

          //if (bsize > asize) {
          //  partitionCfg[ cell ].gain += 0.5 - pow(asize+1, -2)*(1 - asize / (float) (bsize))/2;
          //} else {
          //  partitionCfg[ cell ].gain += 0.5 + pow(asize+1, -2)*(1 - bsize / (float) (asize))/2;
          //}

          //partitionCfg[ cell ].gain += 1 - pow((( (float) asize / (float) (asize + bsize - 1) ) - ( (float) (bsize-1) / (float) (asize + bsize - 1) )), 2 );
          // if (bsize > asize ) {
          //   partitionCfg[ cell ].gain -= ( 1.0 / (float)bsize );
          // }
          // else {
          //   partitionCfg[ cell ].gain += ( 1.0 / (float)asize );
          // }

        }

      }

    }

  }


  for(i=0; i < cfgStruct->num_cell; i++) {
    //printf("cell %d, gain = %0.4f\n", i, partitionCfg[ i ].gain);
  }


}




void KerrighanLinSwap(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct) {

  // choose	node	with	highest	gain	whose movement	would	not	cause	an imbalance
  //   move	nodeto	other	block	and	lock	it

  // Find cell with the highest gain that is unlocked and won't cause an imbalance

  int i,j;
  float max_gain;
  int max_idx;
  int * listA;
  int asize = 0;
  int * listB;
  int bsize = 0;

    // Determine partitioning

  listA = (int *)malloc(cfgStruct->num_cell * sizeof(int));
  listB = (int *)malloc(cfgStruct->num_cell * sizeof(int));

  // DE - FIXME - STORE IN STRUCT SO THAT WE DON'T HAVE TO RECALCULATE
  for(i=0; i < cfgStruct->num_cell; i++) {
    if ( partitionCfg[i].partition == 0 ) {
      listA[asize] = i;
      asize++;
    }
    else {
      listB[bsize] = i;
      bsize++;
    }
  }

  max_gain = (float) -cfgStruct->num_con;

  if (bsize < asize) { // pick from list A

    //printf(" A > B \n");

    // pick the cell in a with the largest gain that is unlocked
    for(i=0; i < asize; i++) {
      if( partitionCfg[ listA[i] ].gain > (float)max_gain &&
          partitionCfg[ listA[i] ].locked == false ) {
        max_idx = listA[i];
        max_gain = partitionCfg[ max_idx ].gain;
      }
    }

  }
  else if (bsize > asize) { // pick from list B

    //printf(" B > A \n");
    //printf(" bsize = %d\n", bsize);

    // pick the cell in b with the largest gain that is unlocked
    for(i=0; i < bsize; i++) {

      if( partitionCfg[ listB[i] ].gain > max_gain &&
          partitionCfg[ listB[i] ].locked == false) {
        max_idx = listB[i];
        max_gain = partitionCfg[ max_idx ].gain;
      }
    }

  }
  else { //bsize == a_size, pick from either
    // pick the cell in a or b with the largest gain that is unlocked

    //printf(" B == A \n");

    for(i=0; i < cfgStruct->num_cell; i++) {
      if( partitionCfg[ i ].gain > max_gain &&
          partitionCfg[ i ].locked == false) {
        max_idx = i;
        max_gain = partitionCfg[ max_idx ].gain;
      }
    }

  }

  // Swap the maximum gain cell
  //printf("swapping cell # %d\n", max_idx);
  partitionCfg[ max_idx ].partition = partitionCfg[ max_idx ].partition ^ 1;
  partitionCfg[ max_idx ].locked = true;

  // for(i=0; i < bsize; i++) {
  //   printf("list B[%d] = %d\n", i, listB[i] );
  // }
  // for(i=0; i < asize; i++) {
  //   printf("list A[%d] = %d\n", i, listA[i] );
  // }
  // printf("\n");

  free(listA);
  free(listB);

}


void unlockNodes( struct partition_cell * partitionCfg, struct exampleData * cfgStruct) {

  int i;

  for(i = 0; i < cfgStruct->num_cell; i++ ) {
    partitionCfg[i].locked = false;
    partitionCfg[i].gain= 0.0;
  }

}


int KerrighanLinStep(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct, int cutSet) {

  int lockCnt = 0;
  int cutSet_temp = 0;
  int cutSet_best;
  struct partition_cell * partitionCfg_best;

  // place cells randomly in partition A and B to start
  partitionCfg_best = (struct partition_cell *)malloc (sizeof (struct partition_cell) * cfgStruct->num_cell);

  cutSet_best = cutSet;
  memcpy(partitionCfg_best, partitionCfg, sizeof (struct partition_cell) * cfgStruct->num_cell);

  // unlock	all	nodes
  unlockNodes(partitionCfg, cfgStruct);

  // while some nodes are unlocked
  while ( lockCnt < cfgStruct->num_cell) {

    // calculate	all	gains
    calcGains(partitionCfg, cfgStruct, schedStruct);

    // choose	node	with	highest	gain	whose movement	would	not	cause	an imbalance
    //   move	node	to	other	block	and	lock	it
    KerrighanLinSwap(partitionCfg, cfgStruct, schedStruct);
    lockCnt++;

    // calculate and cutSet, store partition if this is the best cutSet result
    cutSet_temp = calcCutSet (partitionCfg, cfgStruct);
    //printf("Done swap, cut_set = %d \n", cutSet_temp);
    if (cutSet_temp <= cutSet_best) {
      cutSet_best = cutSet_temp;
      // copy the partition configuration for later
      memcpy(partitionCfg_best, partitionCfg, sizeof (struct partition_cell) * cfgStruct->num_cell);
    }

  }

  // choose	best	cut	seen	in	this	pass
  memcpy(partitionCfg, partitionCfg_best, sizeof (struct partition_cell) * cfgStruct->num_cell);
  return cutSet_best;

}

void initSchedule(struct exampleData * cfgStruct, struct partitionStruct * schedStruct) {


  struct partitionStruct partitionInit;

  //partitionInit.initTemp = 1000.0;
  //partitionInit.currTemp = 1000.0;
  //partitionInit.betaDecay = true;
  //partitionInit.tempBeta = 0.80;
  //partitionInit.tempBetaSwitch = false;
  //partitionInit.tempBeta2 = 0.98;
  //partitionInit.tempIter = (int)( 10.0 * pow( (double)cfgStruct->num_cell, (double)(4.0/3.0) ) );
  partitionInit.maxSteps = 8;
  //partitionInit.stepAccept = false;
  //partitionInit.accRate = 0.44;
  //partitionInit.accCnt = 0;
  //partitionInit.bb_size = MAX( cfgStruct->num_cols, cfgStruct->num_rows);
  //partitionInit.stopCostStd = 0.1;
  //partitionInit.numStopIter = 100;

  *schedStruct = partitionInit;

}

/* change 1D linear world index ix into 2D index in (nx,ny) size world                  */
struct location        index1Dto2D(unsigned int ix                ,unsigned int nx,unsigned int ny){
  struct location index2D;
  if (ix >= (nx*ny)) {
    index2D.x = UINT_MAX;
    index2D.y = UINT_MAX;
    printf("-W- Geometry: Out of range index passed to index1Dto2D\n");
    return index2D;
  }
  index2D.x = ix%nx;
  index2D.y = ix/nx;
  return index2D;
}



/* draw floorplan using EasyGl graphics module. World size is xWorld*yWorld */
void fpDraw(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, float xWorld, float yWorld) {

  unsigned int nx, ny; /* floorplan dimention               */
  unsigned int i, j, i0, i1        ; /* general counters                  */
  unsigned int celli,neti          ; /* cell/net counters                 */
  struct location     pnt                 ; /* point to hold cell location       */
  float        step                ; /* grid step                         */
  char         label[1024]         ; /* general label                     */
  float        scaleStep           ; /* a step for improvement scale      */
  float        improvementRate     ; /* crossing_nets / total_nets        */
  unsigned int prvCell , curCell   ; /* previous/current cell number      */
  unsigned int srcCell0, srcCell1  ; /* source cell for each group number */
  //net          curNet              ; /* current net                       */
  unsigned int palette=0           ; /* net colors                        */

  /* calculate smallest floorplan dimension */
  nx = ceil(sqrt(cfgStruct->num_cell+1));
  if (nx%2 == 1) nx++;
  ny = nx;

  /* draw header, including title */
  setcolor(LIGHTGREY);
  fillrect(0,0,1000,40);
  setcolor(BLACK);
  setfontsize(13);
  drawtext(500,20,"Kernighan-Lin based bi-partitioning algorithm",1000);

  /* draw second line header, including partitioning parameters */
  setcolor(LIGHTGREY);
  fillrect(0,60,880,100);
  setcolor(BLACK);
  setfontsize(10);
  drawtext(440,80,gHeaderLabel,880);

  /* draw improvement scale */
  //scaleStep=800/20;
  //improvementRate = (float)(fp->crossingNetsN)/(float)(fp->netsN);
  //setcolor(LIGHTGREY); fillrect(900,60 ,1000,940);
  //setcolor(WHITE    ); fillrect(950,120,980 ,920);
  //if (fp->isWorst) setcolor(RED); else setcolor(GREEN);
  //fillrect(950,920-(improvementRate/0.05)*scaleStep, 980 ,920);
  //setcolor(BLACK); drawrect(950,120,980 ,920);
  //drawtext(950,75,"Crossing",100);
  //drawtext(950,100,"Total Nets",100);
  //drawline(910,85,990,85);
  //for(i=1;i<=19;i++) {
  //  drawline(940,i*scaleStep+120,980,i*scaleStep+120);
  //  sprintf(label,"%.2f",1-i*.05);
  //  for(j=0; j<4; j++) label[j]=label[j+1]; /* shift label to remove leading zero */
  //  drawtext(920,i*scaleStep+120,label,30);
  //}

  /* draw footer */
  setcolor(LIGHTGREY); fillrect(0,960,1000,1000);
  setcolor(BLACK    ); drawtext(500,980,gFooterLabel,1000);

  /* draw grid */

  /* draw separating line */
  setlinestyle (DASHED); setcolor(RED);
  drawline(440,120,440,900);

  /* draw two groups */
  setlinestyle (SOLID);
  setcolor(LIGHTGREY);
  fillrect(20,140,390,880);
  fillrect(490,140,860,880);
  setcolor(RED);
  drawrect(20,140,390,880);
  drawrect(490,140,860,880);

  /* draw groups/seperator labels */
  //setcolor(BLACK);
  //sprintf(label,"Cells#:%u",fp->groupCellsN[0]);
  //drawtext(205,900,label,370);
  //printf(label,"Cells#:%u",fp->groupCellsN[1]);
  //drawtext(657,900,label,370);
  //sprintf(label,"nets#:%u",fp->crossingNetsN);
  //drawtext(440,920,label,200);


  step= (740./(float)(2*ny+2));

  /* draw cells */
  i0=0; i1=0;
  for(celli=0; celli<cfgStruct->num_cell; celli++) { /*  draw every cell */
    sprintf(label,"%u",celli);
    if (partitionCfg[celli].partition == 0) { /* if cell in group 0, draw in the left group */
      pnt = index1Dto2D(i0, nx/2, ny); /* set 2D location */
      partitionCfg[celli].location.x = 20 +1.5*step+2*pnt.x*step;
      partitionCfg[celli].location.y = 140+1.5*step+2*pnt.y*step;
      setcolor(BLACK);
      i0++;
    }
    else { /* if cell in group 1, draw in the right group */
      pnt = index1Dto2D(i1, nx/2, ny); /* set 2D location */
      partitionCfg[celli].location.x = 490+1.5*step+2*pnt.x*step;
      partitionCfg[celli].location.y = 140+1.5*step+2*pnt.y*step;
      setcolor(BLACK);
      i1++;
    }
    if (partitionCfg[celli].locked) { /* if locked, draw with yellow */
      setcolor(YELLOW);
      fillarc(partitionCfg[celli].location.x, partitionCfg[celli].location.y, step/2, 0., 360.);
    }
    /* draw a cell as a filled circle */
    setcolor(BLACK);
    drawarc (partitionCfg[celli].location.x, partitionCfg[celli].location.y, step/2, 0., 360.);
    drawtext(partitionCfg[celli].location.x, partitionCfg[celli].location.y, label , step    );
  }

  /* draw nets */
  for(neti=0;neti<(cfgStruct->num_con);neti++) {
    //curNet = cfgStruct->net_connect[neti];
    //prvCell = curNet.cells[0];
    srcCell0 = UINT_MAX;
    srcCell1 = UINT_MAX;

    /* find net cells in each groups, if possible */
    for(celli==0; celli<cfgStruct->num_terminals[neti]; celli++) {
      curCell = cfgStruct->net_connect[neti][celli];
        if (partitionCfg[curCell].partition == 0)
          srcCell0 = curCell;
        else
          srcCell1 = curCell;
    }

    setcolor((enum color_types)((palette++)%7+4)); /* cahnge color */
    /* draw line between source cell and all other net cells */
    for(celli==0; celli<cfgStruct->num_terminals[neti]; celli++) {
      curCell = cfgStruct->net_connect[neti][celli];
        if (partitionCfg[curCell].partition == 0) { /* if chosen cell in group 0 */
          drawline(partitionCfg[srcCell0].location.x,partitionCfg[srcCell0].location.y,
               partitionCfg[ curCell ].location.x,partitionCfg[curCell ].location.y);
        } else { /* if chosen cell in group 1 */
          drawline(partitionCfg[srcCell1].location.x, partitionCfg[srcCell1].location.y,
              partitionCfg[curCell ].location.x, partitionCfg[curCell ].location.y);
        }
    }
    /* draw line between the two groups */
    if ( (srcCell0 != UINT_MAX) && (srcCell1 != UINT_MAX) )
      drawline(partitionCfg[srcCell0].location.x,partitionCfg[srcCell0].location.y,
           partitionCfg[srcCell1].location.x,partitionCfg[srcCell1].location.y);
  }

} /* fpDraw */
