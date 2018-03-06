#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <limits.h>
#include <ctype.h>

#include "graphics.h"

struct exampleData {
   int   num_cell;
   int   num_con;
   int   num_rows;
   int   num_cols;
   int * num_terminals;
   int ** net_connect;
   char * filename;
};

struct partitionStruct {
   int     maxSteps;
   int     nIters;
};

struct location {
   int   x;
   int   y;
};

struct partition_cell {
   int     partition;
   bool    locked;
   float   gain;
   struct location location;
};



void readfile (const char* filename, struct exampleData * data);
void printPartition ( struct partition_cell * partitionCfg, struct exampleData * cfgStruct);
int calcCutSet ( struct partition_cell * partitionCfg, struct exampleData * cfgStruct);
int KerrighanLinStep(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct, int cutSet, int step);
void KerrighanLinSwap(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct);
void initSchedule(struct exampleData * cfgStruct, struct partitionStruct * schedStruct);
float calcStd(int * x, int num_samples);
void unlockNodes( struct partition_cell * partitionCfg, struct exampleData * cfgStruct);
void calcGains(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct);
unsigned int commandlineParse(int argc, char *argv[]  );

/***** GUI functionsand global variables  **************/

char       gRunMode        = 'p'     ; /* run mode, default is step     */
char       gFooterLabel[1024]  =""; /* global footer text label   */
char       gHeaderLabel[1024]  =""; /* global header text label   */
char       gFooterMessage[1024]=""; /* global footer text message */
float      gWorldX = 1000         ; /* world X dimension          */
float      gWorldY = 1000         ; /* world X dimension          */
bool       gGUI = false;
bool       gVerbose = false;

struct partition_cell * gPartitionCfg;
struct exampleData * gCfgStruct;

void runStep  (void (*drawScreen_ptr)(void)) {gRunMode='s' ;} /* run one step */
void runPass  (void (*drawScreen_ptr)(void)) {gRunMode='p' ;} /* run one pass */
void runAll   (void (*drawScreen_ptr)(void)) {gRunMode='a' ;} /* run all      */


/* redrawing routine for still pictures. Redraw if user changes the window           */
void fpDraw(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, float xWorld, float yWorld);
void drawScreen() { clearscreen(); fpDraw(gPartitionCfg, gCfgStruct,gWorldX,gWorldY); } /* clear & redraw  */

/* called whenever event_loop gets a button press in the graphics area.              */
void buttonPress  (float x, float y) {                                     }

/* receives the current mouse position in the current world as in init_world         */
void mouseMove    (float x, float y) {                                                }

/* function to handle keyboard press event, the ASCII character is returned          */
void keyPress     (int i) {                                                           }

/* show global message, wait until 'Proceed' pressed then draw screen                */
void waitLoop (struct partition_cell * partitionCfg, struct exampleData * cfgStruct)  {
  update_message(gFooterMessage);
  drawScreen();
  event_loop(buttonPress,drawScreen);
}

/*******************/

// Random number generation

unsigned int getUIntRand(unsigned int minRand, unsigned int maxRand						);
void uarrRandInit( unsigned int *uarr   , unsigned int  uarrSize   ,
                   unsigned int  minRand, unsigned int  maxRand   );
int uarrValueFound( unsigned int *uarr   , unsigned int uarrSize, unsigned int val  ,
                    unsigned int  fromInd, unsigned int toInd                      );


int main(int argc, char **argv) {

  int fileInd;

  fileInd = commandlineParse(argc, argv);

  printf("You have entered file = %s\n",  argv[fileInd]);
  const char* dir = "benchmarks/";
  const char* filename = argv[1];
  struct exampleData cfgStruct;
  struct exampleData * cfgStructPtr;

  struct partitionStruct   schedCfg;
  struct partitionStruct * schedStructPtr;

  bool printPart = false;
//  bool psEnable = false;

//  if (argc > 2 && (int)argv[2] == 1) {/
//    printPart = true;
//  }
//  if (argc > 3  && (int)argv[3] == 1) {
//    psEnable = true;
//  }

  int i, j, k;
  int minCutSet, cutSet, cutSet_new, delta_cutSet, cutSetSum;
  float cutSetAvg = 0.0;

  minCutSet = pow(2,16);

  cfgStructPtr = &cfgStruct;
  char* read_file = (char*)malloc(strlen(dir)+strlen(filename));
  strcpy(read_file, dir);
  strcat(read_file, filename);
  readfile(read_file, cfgStructPtr);
  gCfgStruct = cfgStructPtr;

  if (gVerbose) {
    printf("num_cells = %d \n" , cfgStructPtr->num_cell);
    printf("num_con = %d \n"   , cfgStructPtr->num_con);
    printf("num_rows = %d \n"  , cfgStructPtr->num_rows);
    printf("num_cols = %d \n"  , cfgStructPtr->num_cols);
    printf("num_nets = %d \n"  , cfgStructPtr->num_con);
    for(i = 0; i < cfgStructPtr->num_con; i++ ) {
       printf("net # %d :", i);
       for(j = 0; j < cfgStructPtr->num_terminals[i]; j++) {
         printf(" %d", cfgStructPtr->net_connect[i][j] );
       }
       printf("\n");
    }
  }

   // check that the number of cells will fit within the design
  if (cfgStructPtr->num_cell > cfgStructPtr->num_cols*cfgStructPtr->num_rows) {
    printf("ERROR : %dx%d grid is too small to place %d cells\n", cfgStructPtr->num_cols,
                                                               cfgStructPtr->num_rows,
                                                               cfgStructPtr->num_cell );
    return(0);
  }

  if (gGUI) {

    // initialize display
    init_graphics((char*)"Kernighan-Lin bi-partitioning algorithm");
    init_world (0.,0.,gWorldX,gWorldY);

    /* Create new buttons */
    /* Create new buttons */
    create_button ((char*)"Window"    , (char*)"Run Step"  , runStep  ); /* run one step */
    create_button ((char*)"Run Step"  , (char*)"Run Pass"  , runPass  ); /* run one pass */
    create_button ((char*)"Run Pass"  , (char*)"Run All"   , runAll   ); /* run all      */
  }


  unsigned int * partArr_inds;
  struct partition_cell * partitionCfg;

  // initialize run parameters
  schedStructPtr = &schedCfg;
  initSchedule(cfgStructPtr, schedStructPtr);

  // randomize seed based on time
  srandom(time(NULL));

  // place cells alternately in partition A and B to start
  for(k = 0; k < schedStructPtr->nIters; k++ ) {

    partArr_inds = (unsigned int *)malloc(sizeof(unsigned int) * cfgStructPtr->num_cell );
    uarrRandInit( partArr_inds , cfgStructPtr->num_cell, 0, cfgStructPtr->num_cell-1  );

    partitionCfg = (struct partition_cell *)malloc (sizeof (struct partition_cell) * cfgStructPtr->num_cell);
    for(i = 0; i < cfgStructPtr->num_cell; i++) {
      //partitionCfg[i].partition = i % 2;
      partitionCfg[partArr_inds[i]].partition = i % 2;
      partitionCfg[i].locked = false;
      partitionCfg[i].gain = 0.0;
    }

    gPartitionCfg = partitionCfg;

    // Show initial placement
    if( gVerbose ) {
      printf("Initial partition: \n");
      printPartition ( partitionCfg, cfgStructPtr);
    }

    // Calculate initial cut-set size
    cutSet = calcCutSet ( partitionCfg, cfgStructPtr );
    if (gVerbose) {
      printf("initial cut set is : %d\n", cutSet);
    }

    // display initial partitioning
    if (gGUI) {
      /* update global message and wait for 'Proceed' to proceed */
      sprintf(gHeaderLabel,"Design: %s  |  Cells: %u  |  Nets: %u | Cut Set : %d",
              cfgStructPtr->filename, cfgStructPtr->num_cell, cfgStructPtr->num_con, cutSet);
      if (gRunMode=='f' || gRunMode=='a') { /* show initial floorplan for 'final' mode */
        sprintf(gFooterMessage,"Initial partitioning. Press 'Proceed' to find final solution");
        sprintf(gFooterLabel,"Initial partitioning");
        waitLoop(partitionCfg, cfgStructPtr);
      }
    }

    // K-L loop


    // output csv header
    // fputs("\n", fp);

    for(i = 0; i < schedStructPtr->maxSteps; i++ ) {

      // Perform iteration of K-L algorithm
      cutSet = KerrighanLinStep(partitionCfg, cfgStructPtr, schedStructPtr, cutSet, i);

      // report cutset statistic every iteration
      if (gVerbose) {
        printf("iteration %d, cut set = %d\n", i, cutSet );
      }

    }

    // display final partition on stdout
    if (gVerbose) {
      printf("End partition :\n");
      printPartition( partitionCfg, cfgStructPtr);
    }
    cutSet = calcCutSet( partitionCfg, cfgStructPtr );
    cutSetSum += cutSet;
    if ( cutSet < minCutSet ) {
      minCutSet = cutSet;
    }

    if (gVerbose) {
      printf("end cut set is : %d\n", cutSet);
    }

    // display final partition on GUI

    /* finished! wait still until 'Exit" is pressed */
    // unlock	all	nodes
    unlockNodes(partitionCfg, cfgStructPtr);
    if (gGUI)
      while(1) {
        sprintf(gHeaderLabel,"Design: %s  |  Cells: %u  |  Nets: %u | Cut Set : %d",
                cfgStructPtr->filename, cfgStructPtr->num_cell, cfgStructPtr->num_con, cutSet);

        sprintf(gFooterMessage,"Final partitioning. Press 'Exit' to terminate");
        sprintf(gFooterLabel,"Partitioning after pass# %u",schedStructPtr->maxSteps);
        waitLoop(partitionCfg, cfgStructPtr);
      }

  }

  cutSetAvg = (float) cutSetSum / (float) schedStructPtr->nIters;
  printf("average cut set is : %f\n", cutSetAvg);
  printf("min cut set is : %d\n", minCutSet);

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

 int length = strlen(filename);
 dataRead.filename = (char *)malloc (sizeof (char*) * length);
 strncpy(dataRead.filename, filename, length);

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

// calculated the cut set of all hypergraphs
int calcCutSet ( struct partition_cell * partitionCfg, struct exampleData * cfgStruct) {

  int i,j;
  int cutSet;
  int cell;
  int partition_temp;

  // for reach net, determine the cut set and add to cutset counter

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

// Calculates gains of nodes in the hypergraph
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
    asize = 0;
    bsize = 0;

    for(j = 0; j < cfgStruct->num_terminals[i]; j++) {
      cell = cfgStruct->net_connect[i][j];
      if ( partitionCfg[cell].partition == 0 ) {
         asize++;
      }
      else {
        bsize++;
      }
    }

    // Start with the original K-L gain equation
    float attract_alpha = 1.0;
    float repulse_alpha = 1.0;

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

          partitionCfg[ cell ].gain += ( ( (float) bsize / (float) (asize + bsize -1) ) - ( (float) (asize-1) / (float) (asize + bsize -1 ) ) );

        }
        else { // partition b

          partitionCfg[ cell ].gain += ( ( (float) asize / (float) (asize + bsize -1) ) - ( (float) (bsize-1) / (float) (asize + bsize -1 ) ) );

        }

      }

    }

  }

}



// Indiviation Swap interation
void KerrighanLinSwap(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct) {


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

    for(i=0; i < cfgStruct->num_cell; i++) {
      if( partitionCfg[ i ].gain > max_gain &&
          partitionCfg[ i ].locked == false) {
        max_idx = i;
        max_gain = partitionCfg[ max_idx ].gain;
      }
    }

  }

  // Swap the maximum gain cell
  partitionCfg[ max_idx ].partition = partitionCfg[ max_idx ].partition ^ 1;
  partitionCfg[ max_idx ].locked = true;

  free(listA);
  free(listB);

  if (gGUI && (gRunMode=='s')) {
    waitLoop(partitionCfg, cfgStruct);
  }

}

// unlocks all nodes
void unlockNodes( struct partition_cell * partitionCfg, struct exampleData * cfgStruct) {

  int i;

  for(i = 0; i < cfgStruct->num_cell; i++ ) {
    partitionCfg[i].locked = false;
    partitionCfg[i].gain= 0.0;
  }

}

// Performs one step of the KL-algorithm
int KerrighanLinStep(struct partition_cell * partitionCfg, struct exampleData * cfgStruct, struct partitionStruct * schedStruct, int cutSet, int step) {

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

  if (gGUI && (gRunMode=='p')) {
    sprintf(gHeaderLabel,"Design: %s  |  Cells: %u  |  Nets: %u | Cut Set : %d",
            cfgStruct->filename, cfgStruct->num_cell, cfgStruct->num_con, cutSet);
    if (step == 0) { /* first pass */
      sprintf(gFooterMessage,"Initial partitioning. Press 'Proceed' to preform pass# 1");
      sprintf(gFooterLabel,"Initial partitioning");
    } else {
      sprintf(gFooterMessage,"Initial partitioning. Press 'Proceed' to preform pass# %u",step+1);
      sprintf(gFooterLabel,"partitioning after pass# %u",step);
    }
    waitLoop(partitionCfg, cfgStruct);
  }

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

  partitionInit.maxSteps = 6;
  partitionInit.nIters = 1;

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


  /* draw footer */
  setcolor(LIGHTGREY); fillrect(0,960,1000,1000);
  setcolor(BLACK    ); drawtext(500,980,gFooterLabel,1000);

  /* draw grid */

  /* draw separating line */
  setlinestyle (DASHED); setcolor(RED);
  drawline(440,120,440,900);

  /* draw two groups */
  setlinestyle (SOLID);
  setcolor(WHITE);
  fillrect(20,140,390,880);
  fillrect(490,140,860,880);
  setcolor(BLACK);
  drawrect(20,140,390,880);
  drawrect(490,140,860,880);


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
    for(celli=0; celli<cfgStruct->num_terminals[neti]; celli++) {
      curCell = cfgStruct->net_connect[neti][celli];
        if (partitionCfg[curCell].partition == 0)
          srcCell0 = curCell;
        else
          srcCell1 = curCell;
    }

    setcolor((enum color_types)((palette++)%7+4)); /* cahnge color */
    /* draw line between source cell and all other net cells */
    for(celli=0; celli<cfgStruct->num_terminals[neti]; celli++) {

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


/* commandline arguments parsing        */
/* results used to set global variables */
/* returns filename index in argv       */
unsigned int commandlineParse(int argc, char *argv[]  ){

  /* arguments parsing                                              */
  int argi;                             /* arguments index          */
  int fileNameArgInd=-1;                /* file name argument index */
  for(argi=1;argi<argc;argi++) {        /* check all argument       */
    if (argv[argi][0]=='-') {           /* switch argument          */
      switch (tolower(argv[argi][1])) { /* consider first letter    */

        case 'h': /* help */
          printf("\nKernighan-Lin based bi-partitioning algorithm for hypergraphs\n"     );
          printf("Usage:\n"                                                                       );
          printf("  assignment3 INFILE [OPTIONS]\n"                                              );
          printf("Options:\n"                                                                     );
          printf("  -help       (or -h): Print this message\n"                                    );
          printf("  -gui        (or -g): Enable graphical user interface (GUI) \n"                );
          printf("                       Automatically enabled if postscript is enabled\n"        );
          printf("  -verbose    (or -v): Enable verbose logging,\n"                               );
          printf("                       Automatically enabled if GUI is disabled\n"              );
          printf("  -runmode    (or -r): Run mode, followed by one of the following\n"            );
          printf("                       'step'  (or 's'): display results every one step\n"      );
          printf("                       'pass'  (or 'p'): display results every full pass\n"     );
          printf("                       'final' (or 'f'): display print final result only\n"     );
          printf("                       Default mode is 'step'\n"                                );
          printf("Input file syntax:\n"                                                           );
          printf("  <CELLS#> <NET#> <ROWS#> <COLUMNS#>\n"                                         );
          printf("  <#CELLS_CONNECTED_TO_NET_1> <LIST_OF_CELLS_CONNECTED_TO_NET_1>\n"             );
          printf("  <#CELLS_CONNECTED_TO_NET_2> <LIST_OF_CELLS_CONNECTED_TO_NET_2>\n"             );
          printf("  ...\n"                                                                        );
          printf("  <#CELLS_CONNECTED_TO_NET_n> <LIST_OF_CELLS_CONNECTED_TO_NET_n>\n"             );
          printf("Examples:\n"                                                                    );
          printf("  assignment3 cps.txt (using default options)\n"                               );
          printf("  assignment3 cps.txt -gui \n"      );
          printf("  assignment3 cps.txt -g -r f \n"              );
          exit(1);

        case 'v': /* verbose mode    */
          gVerbose=1; /* set verbose */
          break;

        case 'g': /* GUI mode */
          gGUI=1; /* set GUI  */
          break;

        case 'r': /* run mode      */
          argi++; /* next argument */
          gRunMode = tolower(argv[argi][0]);
          if ((gRunMode != 's')&&(gRunMode != 'p')&&(gRunMode != 'f')) {
            printf("-E- Commandline error: -runmode should be followed by: 'step', 'pass', or 'final'! Exiting...\n");
            exit(-1);
          }
          break;

        default : /* unknown argument */
          printf("-E- unknown argument %s\n",argv[argi]);
          exit(-1);

      } /* switch */
    } else fileNameArgInd = argi; /* file name argument index */
  } /* arguments loop (for) */

  /* check if infile is supplied */
  if (fileNameArgInd<0) {printf(" -E- infile should be supplied\n"); exit(-1);}

  /* return filename index in arguments array */
  return fileNameArgInd;
}


unsigned int getUIntRand(unsigned int minRand, unsigned int maxRand						){
 /* random()								gives a random number in [0,RAND_MAX]		*/
 /* random()%maxRand						gives a random number between [0,maxRand-1]	*/
 /* random()%(maxRand-minRan+1)			gives a random number in [0,maxRand-minRand]*/
 /* random()%(maxRand-minRan+1)+minRand	gives a random number in [minRand,maxRand]	*/
 return ( random() % (maxRand - minRand + 1) ) + minRand;
}

/* generate a random numbers array between minRand and maxRand with different numbers */
void uarrRandInit( unsigned int *uarr   , unsigned int  uarrSize   ,
                   unsigned int  minRand, unsigned int  maxRand   ){
  unsigned int randNum, i;
  for(i=0; i<uarrSize; i++) {
    randNum = getUIntRand(minRand,maxRand);
    while (uarrValueFound(uarr,uarrSize,randNum,0,i-1)) /* until not found previously */
      randNum = getUIntRand(minRand,maxRand);           /* generate new number        */
    uarr[i] = randNum;
  }
}

/* return 1 if 'value' is found in the array between fromInd to toInd, 0 otherwise */
int uarrValueFound( unsigned int *uarr   , unsigned int uarrSize, unsigned int val  ,
                    unsigned int  fromInd, unsigned int toInd                      ){
  unsigned int i;
  if ((fromInd>toInd)||(toInd>=uarrSize)) return 0;
  for(i=fromInd; i<=toInd; i++)
    if (uarr[i] == val) return 1;
  return 0;
}
