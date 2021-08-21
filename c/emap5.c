/***************************************************************************
 *                                                                         *
 *  The EMAP5 hybrid FEM/MoM modeling codes were created at the            *
 *  University of Missouri-Rolla Electromagnetic Compatibility Laboratory. *
 *  They are intended for use in teaching or research.  They may be freely *
 *  copied and distributed PROVIDED THE CODE IS NOT MODIFIED IN ANY MANNER.*
 *                                                                         *
 *  Suggested modifications or questions about these codes can be          *
 *  directed to:										   *
 *       Dr. Todd H. Hubing,                                               *
 *       Department of Electrical and Computer Engineering                 *
 *       University of Missouri-Rolla, 					         *
 *       Rolla, MO  65409.                                                 *
 *       E-mail:  hubing@ece.umr.edu					         *
 *                                                                         *
 *  Neither the authors nor the University of Missouri makes any warranty, *
 *  express or implied, or assumes any legal responsibility for the        *
 *  accuracy, completeness or usefulness of these codes or any information *
 *  distributed with these codes.                                          *
 *                                                                         *
 *   PROGRAM     :  emap5.c						               *
 *   VERSION     :  1.1							               *
 *   AUTHORS     :  Dr. Mohammad W. Ali, Yun Ji, Dr. Todd T. Hubing        *
 *   LAST UPDATE :  Dec 16,  1998 					               *
 *   DESCRIPTION :  A 3-D hybrid FEM/MoM Code for Analyzing                *
 *                  Time Varying Complex Electromagnetic Fields.           *
 *									                     *
 ***************************************************************************/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "complex.c"
#include "alloc.c"
#include "util.c"
										

/***************************  MARCRO DEFINITION ***************************/

/* TOLERANCE can be changed to 1E-3 or 1E-4  */
#define  TOLERANCE 1E-4

#ifndef M_PI
#define M_PI 3.1415926
#endif 

/****** Don't change the following four marcros ******/

/*  Use 7 point Gaussian Quadrature to evaluate elements of MOM matrices */
#define  TotQuadPoint 7

/* the maximum non-zero elements of B matrix */
#define  MaxBmatElementNum 5

/*  the maximum non-zero elements of A matrix */ 
#define  MaxAmatElementNum 30

/* the maximum forced edge number */
#define  MaxFORCDMatElementNum 15




/*************************** FUNCTION PROTOTYPE ****************************/

void ComputeBddMatrix(int, int, double *, double *, int *, int *,
                      int *, int *, double *, double [][3]);
void ComputeCMatrix(int, int,  double *, double *, int *, int *,
                    int *, int *, double **);
void ComputeCorrectionTerm();
void ComputeDMatrix(int, int, double *, double *, int *, int *,
                    int *, int *, double **, double [][3]);
void ComputeGaussQuadPoint(int, int *, double *);
void ComputeParameters();
void ComputeSingul(int *, double *, double, double*, double *, double *); 
void ComputeSourceVector( double, double *, complex *);
void ComputeNonSingul(int *, double, double *,complex *, 
                      complex *, complex *, complex *);
void ComputeNonSingul1(int *, double *, complex *, complex *, 
                       complex *, complex *);
void ComputeTetraHedronVolume(int,double **,int **,double *);
void ComputeTrngleCentroid(int, double **, int **, double **);
void ConjugateSolver(int,int *, int **, complex **, complex *,int, double);
void CreateRHSVector(int ,int, int, double *, complex *);
void FemMatrixCompute(int,int, int **, int **,double **, complex *, 
                      complex *, int **, complex **, int *, double *);
complex InnerProduct(complex *, complex *,int );
void InvertMatrix(complex **, int);
void matrixVectorproduct(char, int, complex *, complex *, 
                         int i, int **, complex **);
int  SearchNonZeroElement(int, int);
void SurfaceFieldCompute();
void TetFaceAreaNormal(int,double **,int **, double **,double *);
void PrintOutput();
void ReadInputFile();
double Sign(int);
double VectorNorm(complex *,int);



/****************************** GLOBAL VARIABLES ****************************/



/****************************************************************************
       the following varialbes are used to represent resistor edges 
*****************************************************************************/
/* to represent */
int ResistorEdgeStat[100];      /* the resistor edge index */
double ResistorEdgeValue[100];  /* the resistor value */
int ResistorEdgeNum;            /* total number of resistor edges */


/****************************************************************************
       the following varialbes are used to represent junctions 
*****************************************************************************/
int JunctionNum=0;     /* total number of junctions, initialize to zero */

/* Junction[][0] ---- the hybrid edge number 
   Junction[][1] ---- the other external metal edge with same coordinates 
                      with Junction[][0]. these two form a juction 
   Junction[][2] ---- the T++ or T-- triangle associate with the junction  */

int Junction[100][3];  /* maxium junction edge allowed is 100. */


char SourceType;   /* Source type. 
                      'V' means voltage source
                      'P' means incident plane wave source
                      'I' meams current source  */

/****************************************************************************
       the following varialbes are used to represent voltage sources 
 ****************************************************************************/
int VSourceNum,         /* voltage source edge number */
    VSourceEdge[20];    /* voltage source edge table. The maximum
                           voltage edges allowed is 20  */
double VSourceMag[20];  /* voltage edge magnitude */


/****************************************************************************
       the following varialbes are used to represent plane wave 
 ****************************************************************************/
double K_Theta, K_Phi,  /* the unit vector of K in the spherical  
                           coordiantes   */

       E_Theta, E_Phi,  /* the unit vecotr of E in the sperical  
                                    coordinates */
       PlaneWaveE_Mag;  /* the E field magnitude of the plane wave */


int 

    DielBoundEdgeNum,     /* total number of the diectric boundary edges */
    HybrdBoundEdgeNum,    /* total number of the hybrid boundary edges */
    TotEdgeNum,           /* total number of edges */
    TotTetElement,        /* total number of tetrahedrons */
    TotNodeNum,           /* total number of nodes */
    TotTrngleNum,         /* total number of triangles */
    TotBoundEdgeNum,      /* total number of boudary edgs */
    TotInnerEdgeNum,      /* total number of inner edgs */
    TotExtMetalEdgeNum,   /* total number of the exteral metal-boundary edges */

    TotISourceEdgeNum,
    TotMatrixEqnNum,      /* the size of final matrix equation after MoM is  
                             coupled to FEM matrix equation */

    **TrngleNode,      /* triangle node table */
    **TrngleEdge,      /* triangle edge table */

    *BoundEdgeStat,       /* boudary edge table */
    *InnerEdgeStat,       /* inner edge table */
    *ISourceEdgeStat,       /* the current source edge table */

    **PlusTrngleDetect,      /* plus triangle detected table */
    *PlusTrngleIndex,     /* plux triangle table */
    **MinusTrngleDetect,     /* minus triangle detected table */
    *MinusTrngleIndex,    /* minus triangle table */
    *VcRowNum ,
    **GBLmatColAddr,      /* the column table for each row in FEM matrix, which
                             is used by the row-indexed scheme */
    *GBLmatColIndex,      /* the number of non-zero element in each row in FEM
                             matrix, which is used by the row-indexed scheme */
    **GlobalEdgeEnds,  /* glabal edge table */

    **TetNode,   /* tetrahedron node table */
    **TetEdge;   /* tetrahedron edge table */




double

    OperateFreq,    /* the working frequency, the unit is Hz */
    FreeSpaceVel,   /* the free space light velocity: 3X10^8 m/s */
    AbsPermeable,   /* the absolute permiability: 4*Pi*10^-7  */
    AbsPermitt,     /* the absolute permittivity: 8.854*10^-12 */
    WaveNumber,     /* wavenumber, defined as wavelength/(2*pi) */
    WaveLength,     /* wavelength in meters */
    ImpeDance,      /* the impedance of free space, which equals to 377 ohms */
    *Wght,          /* the weighting coefficients of Gassian Qadrature */
    **Qpnt,         /* the position coefficients of Gassian Qadrature */
    *EdgeLength,    /* edge length table */
    *TetVolume,     /* tetrahedron volume table */
    **NodeCord,     /* global node table */
    **TrngleNormal, /* the normal vector of each triangle. TrngleNormal[i][0],
                       TrngleNormal[i][1] and TrngleNormal[i][2] are the x, y
                       and z component of the normal vector of triangle i+1. */

    **TrngleCenTroid,  /* the centroid point of each triangle.
                          TrngleCenTroid[i][0], TrngleCenTroid[i][1], 
                          and TrngleCenTroid[i][2] are the x, y and z 
                          coordinates of the centroid of triangle i+1 */


    *TrngleArea,    /* the triangle area table. TrngleArea[i] is the
                       area of triangle i+1 */
    **Bdd,          /* the surface term in the FEM equation. */
    *ISourceEdgeMag;   /* the current source value */


complex

    CFactor1,       /* one constant associated with MoM equation, see
                       ComputeParameters() */
    CFactor2,       /* one constant associated with MoM equation, see
                       ComputeParameters() */
    TParam[2],      /* two constants associated with FEM equation, see
                       ComputeParameters() */
    *Epsilon,       /* the complex permittivity table of each tetrahedron */
    *RHSVector,     /* the Right-Hand Side vector (RHS) of the 
                       final hybrid matrix equation */
    **Cdd, **Cdc, **Ccd, **Ccc,   /* four submatices of the MoM C matrix */
    *Einc,          /* the incident wave term generated by the testing
                       functions.  */  
    **Ddd, **Dcd,   /* the two submatrices of the MoM D matrix */
    **GBLmat,       /* the FEM matrix data stored in row-indexed scheme */


    /*********************************************************************
       the following six vectors are used in ComputerCorrectionTerm().
       They are used to couple the MoM equation to the FEM equation
     *********************************************************************/
    *GdVector, /* ComputerCorrectionTerm() CreateRHSVector */
    *VcVector, /* obtained in main, used in ComputerCorrectionTerm() */
    *VdVector, /* obtained in main, used in ComputerCorrectionTerm() */
    *UdVector, /* ComputerCorrectionTerm()*/
    *MdVector, /* ComputerCorrectionTerm(), SurfaceFieldCompute() */
    *McVector, /* ComputerCorrectionTerm(), SurfaceFieldCompute() */  

    *EdVector,   /* the magnetic current density on the dielectric surface */
    *JcVector,   /* the electric current density on the conductive surface */
    *JdVector;   /* the electric current density on the dielectric surface */


time_t start_time, mid_time, end_time;  /* to count CPU time used */
char  log_mesg[300];     /* store the log information generated in 
                            ConjugateSolver( ) */

FILE *InF,*OutF,*OutF1;  /* InF: file pointer of the input file 
                            OutF: file pointer of the output file
                            OutF1: file pointer of the log file  */



/**************************************************************************
 *		 	   MAIN PROGRAM 				  *
 **************************************************************************/
int main()
{
    int IEdge, JEdge, ObservTrngle, SourceTrngle,     
        ObservEdge, QuadPointM, QuadPointN;
    char filename[40];
    double **SrcPointArray, **ObsPointArray, ***RhoCtrd;
    char argv[40];
    strcpy( argv, "data/E40.in" );


    // if( argc!=2 )
    // {
    //     fprintf(stderr, "Usage:emap5 inputfile\n");
    //     exit(1);
    // }
 
    if( strlen(argv)<=3 ) 
    {
        fprintf(stderr, "please use *.in  as the input file\n");
        exit(1);
    }

    /* get the surfix */
    strcpy( filename, argv+strlen(argv)-3 );

    if( strcmp(filename, ".in")!=0 ) 
    {
        fprintf(stderr, "please use .in surfix for the input file\n");
        exit(1);
    }

    /* open the input file */
    if( (InF = fopen(argv,"r"))==NULL )
    {
        fprintf(stderr, "Sorry, can't open input file\n");
        exit(1);
    }
 

    /* open the log file */
    strncpy( filename, argv, strlen(argv)-3 );
    filename[ strlen(argv)-3 ] ='\0';    /* add numm character */

    strcat(filename, ".log");
 
    if( (OutF1=fopen(filename, "w"))==NULL )
    {
        fprintf(stderr, "Sorry, can't open log file\n");
        exit(1);
    } 

    ReadInputFile();

    if (SourceType=='V')
    {     
        int flag=0;

        for(IEdge=0;IEdge<VSourceNum;IEdge++)
        {
        
 	    for(JEdge=0;JEdge<TotBoundEdgeNum;JEdge++)  
            {
                if( BoundEdgeStat[JEdge]==VSourceEdge[IEdge] ) 
                {

                    /* vsource can not be defined on the hybrid boudary */
                   
                    if( JEdge<HybrdBoundEdgeNum ){
                        printf("\nvsource can not be defined on\
 the hybrid boudary\n Sorry, exit!!");
                        exit(1);
                    }

                    flag=1;
                    break;

                }  /*  if(BoundEdgeStat[JEdge]==  */
            }     
            
            /* Can't find the Vsource edge */
            if( flag==0 ) 
            {
                printf("Sorry: Check your Vsource edge\n");
                exit(1);
            } /* end of if(flag) */
    
        } /*for(IEdge) */

    }  /* if (!strcmp(  */
    

  
    /*******************************************************/
    TotExtMetalEdgeNum = TotBoundEdgeNum - HybrdBoundEdgeNum;
  
    /* revise here on Feb 2 */
    TotMatrixEqnNum = TotInnerEdgeNum + HybrdBoundEdgeNum ;
  
     
    /********************************************************************* 
     *    BOUNDARY ELEMENT PART OF THE HYBRID CODE                       *
     *********************************************************************/

    ObsPointArray = double_Matrix(TotQuadPoint,3);
    SrcPointArray = double_Matrix(TotQuadPoint,3);
    Qpnt          = double_Matrix(TotQuadPoint,3);
    Wght          = double_Vector(TotQuadPoint);
    RhoCtrd       = double_Matrix2(TotTrngleNum,3,3);

    Wght[0]=0.112500000000000; 
    Wght[2]=Wght[1]= Wght[3]=0.062969590272413; 
    Wght[6]=Wght[5]= Wght[4]=0.066197076394253; 

    Qpnt[0][1]=Qpnt[0][0]=0.333333333333333;  
    Qpnt[2][1]=Qpnt[1][0]=0.797426985353087; 
    Qpnt[3][0]=Qpnt[3][1]=Qpnt[2][0]=Qpnt[1][1]=0.101286507323456;
  
    Qpnt[6][1]=Qpnt[5][0]=Qpnt[4][1]=Qpnt[4][0]=0.470142064105115;
    Qpnt[6][0]=Qpnt[5][1]=0.059715871789770;
      
   
    for(QuadPointM=0;QuadPointM<TotQuadPoint;QuadPointM++)
        Qpnt[QuadPointM][2] = 1.0-Qpnt[QuadPointM][0]-Qpnt[QuadPointM][1];

    /*********** malloc menory *********************/

    Cdd = CMPLX_Matrix(HybrdBoundEdgeNum,HybrdBoundEdgeNum);
    Cdc = CMPLX_Matrix(HybrdBoundEdgeNum,TotExtMetalEdgeNum);
    Ccd = CMPLX_Matrix(TotExtMetalEdgeNum,HybrdBoundEdgeNum);
    Ccc = CMPLX_Matrix(TotExtMetalEdgeNum,TotExtMetalEdgeNum);

    Null_CMPLXMatrix(Cdd,HybrdBoundEdgeNum,HybrdBoundEdgeNum);
    Null_CMPLXMatrix(Cdc,HybrdBoundEdgeNum,TotExtMetalEdgeNum);
    Null_CMPLXMatrix(Ccd,TotExtMetalEdgeNum,HybrdBoundEdgeNum);
    Null_CMPLXMatrix(Ccc,TotExtMetalEdgeNum,TotExtMetalEdgeNum);

    Ddd = CMPLX_Matrix(HybrdBoundEdgeNum,HybrdBoundEdgeNum);
    Dcd = CMPLX_Matrix(TotExtMetalEdgeNum,HybrdBoundEdgeNum);
    Bdd = double_Matrix(HybrdBoundEdgeNum,HybrdBoundEdgeNum);

    Null_CMPLXMatrix(Ddd,HybrdBoundEdgeNum,HybrdBoundEdgeNum);
    Null_CMPLXMatrix(Dcd,TotExtMetalEdgeNum,HybrdBoundEdgeNum);
    Null_double_Matrix(Bdd,HybrdBoundEdgeNum,HybrdBoundEdgeNum);


    if( SourceType=='P' )
    {
        Einc =  CMPLX_Vector(TotBoundEdgeNum);
        Null_CMPLXVector(Einc,TotBoundEdgeNum);
    }

    ComputeParameters(); 
    
    mid_time=start_time=time(0);   /* start timer */
 
    for(ObservTrngle=0;ObservTrngle<TotTrngleNum;ObservTrngle++)
    {
        int RowNum, MultRow,RowRootIndxNum,RowAddlIndxNum, RowCurrIndxNum, 
            RowStart, RowEnd;
    
        double SgnM;

        double SrcTrngleEdgeLen[3], ObsTrngleEdgeLen[3], SrcTrngleNorm[3];
        int SourceTrngleNode[3], SourceTrngleEdge[3], ObservTrngleNode[3],
            ObservTrngleEdge[3];


        for(IEdge=0;IEdge<=2;IEdge++)
        {
            double Buff[3];

            ObservTrngleNode[IEdge] = TrngleNode[ObservTrngle][IEdge]-1;
            ObservTrngleEdge[IEdge] = TrngleEdge[ObservTrngle][IEdge];
            ObsTrngleEdgeLen[IEdge] = 
                       EdgeLength[abs(ObservTrngleEdge[IEdge])-1];
            VTXsub(TrngleCenTroid[ObservTrngle], 
                   NodeCord[ObservTrngleNode[IEdge]], Buff);
            RhoCtrd[ObservTrngle][IEdge][0] = Buff[0];
            RhoCtrd[ObservTrngle][IEdge][1] = Buff[1];
            RhoCtrd[ObservTrngle][IEdge][2] = Buff[2];
        }


        for(SourceTrngle=0;SourceTrngle<TotTrngleNum;SourceTrngle++)
        {
            double BE[3][3];

            for(JEdge=0;JEdge<=2;JEdge++)
            {
                SourceTrngleNode[JEdge]=TrngleNode[SourceTrngle][JEdge]-1;
                SourceTrngleEdge[JEdge]=TrngleEdge[SourceTrngle][JEdge];
                SrcTrngleEdgeLen[JEdge] 
                    =EdgeLength[abs(SourceTrngleEdge[JEdge])-1];
                SrcTrngleNorm[JEdge]=TrngleNormal[SourceTrngle][JEdge];
            }

 
            for(QuadPointM=0;QuadPointM<=TotQuadPoint-1;QuadPointM++)
            ComputeGaussQuadPoint( QuadPointM, ObservTrngleNode, 
                                   ObsPointArray[QuadPointM]);


            for(QuadPointN=0;QuadPointN<=TotQuadPoint-1;QuadPointN++)
                ComputeGaussQuadPoint( QuadPointN,SourceTrngleNode, 
                     SrcPointArray[QuadPointN]);

            ComputeCMatrix(ObservTrngle,     SourceTrngle, 
                           ObsTrngleEdgeLen, SrcTrngleEdgeLen, 
                           ObservTrngleNode, SourceTrngleNode,
                           ObservTrngleEdge, SourceTrngleEdge,
                           ObsPointArray);
 
            ComputeBddMatrix(ObservTrngle,     SourceTrngle,
                             ObsTrngleEdgeLen, SrcTrngleEdgeLen, 
                             ObservTrngleNode, SourceTrngleNode,
                             ObservTrngleEdge, SourceTrngleEdge,
                             SrcTrngleNorm,    BE);

            ComputeDMatrix(ObservTrngle,     SourceTrngle, 
                           ObsTrngleEdgeLen, SrcTrngleEdgeLen, 
                           ObservTrngleNode, SourceTrngleNode,
                           ObservTrngleEdge, SourceTrngleEdge,
                           ObsPointArray,    BE);

        } /* end of SourceTrngle --- the big loop */
 
        /****************Compute Eincident field Vector *************/
        if( SourceType=='P' )
        { 
            complex  QE[3],  EatCentroid[3];
            double MidPoint[3];

            MidPoint[0] = TrngleCenTroid[ObservTrngle][0];
            MidPoint[1] = TrngleCenTroid[ObservTrngle][1];
            MidPoint[2] = TrngleCenTroid[ObservTrngle][2];
            ComputeSourceVector(WaveNumber,MidPoint,EatCentroid);

            for(ObservEdge=0;ObservEdge<=2;ObservEdge++)
            {
                complex V1;

                SgnM   = Sign(ObservTrngleEdge[ObservEdge]);
                V1 = COMplex_Add2(
                         Real_Mul(RhoCtrd[ObservTrngle][ObservEdge][0], 
                                  EatCentroid[0]),
                         Real_Mul(RhoCtrd[ObservTrngle][ObservEdge][1], 
                                  EatCentroid[1]),
                         Real_Mul(RhoCtrd[ObservTrngle][ObservEdge][2], 
                                  EatCentroid[2]));
                         QE[ObservEdge]= Real_Mul( 
                                (ObsTrngleEdgeLen[ObservEdge]/2.0), V1);
                         QE[ObservEdge]=Real_Mul((SgnM),QE[ObservEdge]);
            }
 
            for(ObservEdge=0;ObservEdge<=2;ObservEdge++)
            {
  
                RowNum = abs(ObservTrngleEdge[ObservEdge]);
                SgnM = Sign(ObservTrngleEdge[ObservEdge]);
                MultRow= GlobalEdgeEnds[RowNum-1][2];
                RowRootIndxNum= GlobalEdgeEnds[RowNum-1][3];
                RowAddlIndxNum= GlobalEdgeEnds[RowNum-1][5];
                RowCurrIndxNum= RowRootIndxNum+RowAddlIndxNum;
                
                /* don't need to compute this edge. Thus jump to the next */
                if(MultRow<=0) continue;
                
                /*  the following codes are for MultRow >0 */
                if( SgnM >=0.0 )
                {
                    for(IEdge=0;IEdge<PlusTrngleIndex[RowNum-1];IEdge++)
                        if ((PlusTrngleDetect[RowNum-1][IEdge]-1)==ObservTrngle)
                        {
                            RowStart = RowCurrIndxNum+IEdge; 
                            /* revise on Jan 27 */
                            /*RowEnd = RowStart + MultRow;*/
                            RowEnd = RowStart +1;
                            break;
                        }
                        else continue;
                }
                else
                {
                    for(IEdge=0;IEdge<MinusTrngleIndex[RowNum-1];IEdge++)
                        if ((MinusTrngleDetect[RowNum-1][IEdge]-1)
                             ==ObservTrngle )
                        {
                            RowStart = RowCurrIndxNum; 
                            RowEnd = RowStart + 1; break;
                        }
                        else continue;
                }
 
                if( MultRow>0 )
                    for(IEdge=RowStart;IEdge<RowEnd;IEdge++)
                    {
                        Einc[IEdge] = COMplex_Add(Einc[IEdge],QE[ObservEdge]);
                        Einc[IEdge] = Real_Mul(PlaneWaveE_Mag, Einc[IEdge] );
                    }  /* for(IEdge */
                       
            }      /* for ObservEdge */ 
        }          /* if(!strcmp)  */
 
    }  /* for(ObservTrngle=0 */ 

    end_time=time(0);
    fprintf(OutF1, "Time Usage Report:\n\n");
    fprintf(OutF1, "\tcomputing Bdd, C, D matrices\n");
    fprintf(OutF1, "\t\ttime used: %d sec\n\n", (end_time-mid_time));
    mid_time=end_time;

    VcRowNum = INT_Vector(VSourceNum);
    VdVector =  CMPLX_Vector(HybrdBoundEdgeNum);
    GdVector =  CMPLX_Vector(HybrdBoundEdgeNum);
    UdVector =  CMPLX_Vector(HybrdBoundEdgeNum); 
    MdVector =  CMPLX_Vector(HybrdBoundEdgeNum);
    McVector =  CMPLX_Vector(TotExtMetalEdgeNum);
    VcVector =  CMPLX_Vector(TotExtMetalEdgeNum);
 
    Null_CMPLXVector(VdVector,HybrdBoundEdgeNum);
    Null_CMPLXVector(GdVector,HybrdBoundEdgeNum);
    Null_CMPLXVector(UdVector,HybrdBoundEdgeNum); 
    Null_CMPLXVector(MdVector,HybrdBoundEdgeNum);
    Null_CMPLXVector(McVector,TotExtMetalEdgeNum);
    Null_CMPLXVector(VcVector,TotExtMetalEdgeNum);

    if( SourceType=='P' )
    { 
        for(IEdge=0;IEdge<HybrdBoundEdgeNum;IEdge++) 
            VdVector[IEdge] = Einc[IEdge];
        for(IEdge=HybrdBoundEdgeNum;IEdge<TotBoundEdgeNum;IEdge++)
            VcVector[IEdge-HybrdBoundEdgeNum] = Einc[IEdge];
    } 

    if( SourceType=='V' )
    {     
        int flag=0;

        for(IEdge=0;IEdge<VSourceNum;IEdge++)
        {
        
            for(JEdge=0;JEdge<TotBoundEdgeNum;JEdge++)
 	    {
 		if( BoundEdgeStat[JEdge]==VSourceEdge[IEdge] )
                {

                    /* vsource can not be defined on the hybridboudary */
                   
                    if( JEdge<HybrdBoundEdgeNum )
                    {
                        printf("\nvsource can not be defined \
on the hybrid boudary\n Sorry, exit!!");
                        exit(1);

                    }   /* if (JEdge */

                    VcRowNum[IEdge]=JEdge-HybrdBoundEdgeNum;

                    flag=1;
                    break;

                }  /* if( BoundEdgeStat */
            } /*  for(JEdge) */
         
            /* Can't find the Vsource edge */
            if( flag==0 ) 
            {
                printf("Sorry: Check your Vsource edge\n");
                exit(1);

            } /* if(flag) */
    
  
        } /* for(IEdge) */
    
        for(IEdge=0;IEdge<VSourceNum;IEdge++)           
            VcVector[VcRowNum[IEdge]]= COMplex_Cmplx( 
                     VSourceMag[IEdge]*EdgeLength[VSourceEdge[IEdge]-1],0.0);
          
    }  /* if( !strcmp */   


    ComputeCorrectionTerm();

    end_time=time(0);
    fprintf(OutF1, "\tcomputing correction term\n");
    fprintf(OutF1, "\t\ttime used: %d sec\n", (end_time-mid_time));
    mid_time=end_time;

    /**************************************************************/
    /********* FINITE ELEMENT PART OF THE HYBRID CODE *************/
    /**************************************************************/

    GBLmatColIndex=INT_Vector(TotEdgeNum);
    GBLmat= CMPLX_Matrix(TotEdgeNum,MaxAmatElementNum);
    GBLmatColAddr=INT_Matrix(TotEdgeNum,MaxAmatElementNum);
    Null_INTVector(GBLmatColIndex,TotEdgeNum);
    Null_CMPLXMatrix(GBLmat,TotEdgeNum,MaxAmatElementNum);

    for(IEdge=0;IEdge<TotEdgeNum;IEdge++)  
        for(JEdge=0;JEdge<MaxAmatElementNum;JEdge++) 
            GBLmatColAddr[IEdge][JEdge] = -1;
 
    FemMatrixCompute(TotTetElement, MaxAmatElementNum, TetEdge, 
              GlobalEdgeEnds, NodeCord, TParam, Epsilon, GBLmatColAddr,
              GBLmat, GBLmatColIndex,TetVolume);

    /*   consider the resistor   */
    /*   the general form is L^2/ZL   */
    {  
        int i;
   
        for(i=0; i<ResistorEdgeNum; i++)
        {
            int j,flag, RowIndex, ColIndex;

            /* find the row index to the FEM matrix of the resistive edge */ 
            flag=0;
            for(j=0;j<TotInnerEdgeNum;j++) 
                if( InnerEdgeStat[j]==ResistorEdgeStat[i] )
                { 
		    flag=1; 
                    RowIndex=j;  /* row index to the FEM matrix */
                    break;
                }

            if( !flag ) 
            {   fprintf(stderr, "The resistive edge can be only \
defined within the FEM region\n");
                exit(1);
            }

            /* the column index */
            ColIndex=SearchNonZeroElement(RowIndex, RowIndex); 

            GBLmat[RowIndex][ColIndex]=COMplex_Add(GBLmat[RowIndex][ColIndex], 
              COMplex_Cmplx( 
                pow( EdgeLength[ResistorEdgeStat[i]-1],2.0) /ResistorEdgeValue[i], 0.0) ); 

        }  /* for(i=0  */

    } /* resistor */

 
    /**************CREATE RHS VECTOR ********************/

    RHSVector =  CMPLX_Vector(TotMatrixEqnNum);
    Null_CMPLXVector(RHSVector,TotMatrixEqnNum);
 
    end_time=time(0);
    fprintf(OutF1, "\n\tcomputing A matrix\n");
    fprintf(OutF1, "\t\ttime used: %d sec\n", (end_time-mid_time));
    mid_time=end_time;
 
    /*revise here */
    CreateRHSVector(TotInnerEdgeNum,HybrdBoundEdgeNum,TotISourceEdgeNum,
                    ISourceEdgeMag,RHSVector);  


    /**************SOLVING HYBRID MATRIX ********************/
 
    ConjugateSolver(TotMatrixEqnNum, GBLmatColIndex,GBLmatColAddr, GBLmat, 
                    RHSVector,TotMatrixEqnNum, TOLERANCE);

    end_time=time(0);
    fprintf(OutF1, "\n\tsolving the hybrid matrix equation\n");
    fprintf(OutF1, "\t\ttime used: %d sec\n", (end_time-mid_time));
    mid_time=end_time;


    /******* COMPUTATION OF SURFACE FIELDS **************/

    EdVector =  CMPLX_Vector(HybrdBoundEdgeNum);
    JdVector =  CMPLX_Vector(HybrdBoundEdgeNum);
    JcVector =  CMPLX_Vector(TotExtMetalEdgeNum);
    Null_CMPLXVector(EdVector,HybrdBoundEdgeNum);
    Null_CMPLXVector(JdVector,HybrdBoundEdgeNum);
    Null_CMPLXVector(JcVector,TotExtMetalEdgeNum);

    SurfaceFieldCompute();
    end_time=time(0);
    fprintf(OutF1,"\n\tcomputing surface currents:\n");
    fprintf(OutF1, "\t\ttime used: %d sec\n", (end_time-mid_time));
    mid_time=end_time;

    PrintOutput();

    end_time=time(0);
    fprintf(OutF1, "\n        total time used:   %d sec\n",   
                   (end_time-start_time) );

    fprintf(OutF1, "\n%s\n", log_mesg);
           
    fclose(OutF1);
    fclose(InF);

    return 0;


} /*  end of main()  */





/****************************************************************************   
Prototype:  int  SearchNonZeroElement(int RowNum, int ColNum) 
Description:  To search an element in the global FEM matrix. To conserve
              memory, only none-zero elements are stored in the FEM matrix 
              using row-indexed scheme.  If the element is found, return the  
              column index to the global FEM matrix. Otherwise, return -1. 
Input value: 
    int RowNum --- row number, which is assigned to the observing edge 
    int ColNum ---  column number, which is assigned to the source edge. 
Return value:  the column index to GBLmat. 
Global value used: int *GBLmatColIndex, int *GBLmatColAddr 
Global value modified: none 
Subroutines called: none 
*****************************************************************************/
int SearchNonZeroElement(int RowNum,int ColNum)
{
    int Count_i,Count_j;
 
    for(Count_i=0;Count_i<GBLmatColIndex[RowNum];Count_i++)
        if( ColNum==(Count_j=GBLmatColAddr[RowNum][Count_i]) ){ 
            return Count_i;
            break;
        }
    return -1;       
}




/****************************************************************************
Prototype:  void  ComputeGaussQuadPoint(int QuadPoint, int *TrngleNode,
                                            double *SrcPointCol) 
Description: To compute the coordinates of 7-point Gauss nodes of 
             a triangular patch. 
Input value: 
    int QuadPoint  --- node index, it can be from 0 to 6. 
    int *TrngleNode  ---  the three nodes of a tringular patch. 
    double *SrcPointCol --- where to store the results 
Return value: none 
Global value used: NodeCord, Qpnt 
Global value modified:    none 
Subroutines called:    none 
Note: Not very sure.            *****************************************************************************/
void ComputeGaussQuadPoint(int QuadPoint, int *TrngleNode, double *SrcPointCol)
{
    SrcPointCol[0] =   Qpnt[QuadPoint][0]*NodeCord[TrngleNode[0]][0]
                     + Qpnt[QuadPoint][1]*NodeCord[TrngleNode[1]][0]
                     + Qpnt[QuadPoint][2]*NodeCord[TrngleNode[2]][0];

    SrcPointCol[1] =   Qpnt[QuadPoint][0]*NodeCord[TrngleNode[0]][1]
                     + Qpnt[QuadPoint][1]*NodeCord[TrngleNode[1]][1]
                     + Qpnt[QuadPoint][2]*NodeCord[TrngleNode[2]][1];

    SrcPointCol[2] =   Qpnt[QuadPoint][0]*NodeCord[TrngleNode[0]][2]
                     + Qpnt[QuadPoint][1]*NodeCord[TrngleNode[1]][2]
                     + Qpnt[QuadPoint][2]*NodeCord[TrngleNode[2]][2];
}
  



/****************************************************************************
Prototype: void ComputeSingul(int *Node, double *MidPoint, double AREA)
Description:  To evaluate the singularity of the MOM integral analytically when
              the observing triangle coincides with the source triangle. The  
              results are stored in global variables RealIno, RealIxi and 
              RealIeta. 
Input value: 
    int *Node --- the nodes of the triangle 
    double *MidPoint --- middle point of the triangle 
    double AREA --- the area of the triangle 
Return value: none 
Global value used: double **NodeCord 
Global value modified: none
Subroutines called: none 
*****************************************************************************/
void ComputeSingul(int *Node, double *MidPoint, double AREA,                 
                   double *RealIno, double *RealIxi, double *RealIeta)
{
    double RIO,RIP,RIQ;
    double RO[3],RJ[3],RHO[3],HATL[3][3],SMLRO[3][3];
    double X,Y,Z,R,RX,RY,RZ,A,B,C,D,E,F,ABE,ACF,BDF,RJ1,RJ2,RJ3,RJ4,RD;
    int I1,I2,I3,L,L1,K,K1;

    I1 = Node[0];
    I2 = Node[1];
    I3 = Node[2];
    RIO = 0.0;
    RIP = 0.0;
    RIQ = 0.0;
    RX = MidPoint[0];
    RY = MidPoint[1];
    RZ = MidPoint[2];
 
    X = (NodeCord[I2][1]-NodeCord[I1][1])*(NodeCord[I3][2]-NodeCord[I1][2])
         - (NodeCord[I3][1]-NodeCord[I1][1])*(NodeCord[I2][2]-NodeCord[I1][2]);
    Y = (NodeCord[I3][0]-NodeCord[I1][0])*(NodeCord[I2][2]-NodeCord[I1][2])
         - (NodeCord[I2][0]-NodeCord[I1][0])*(NodeCord[I3][2]-NodeCord[I1][2]);
    Z = (NodeCord[I2][0]-NodeCord[I1][0])*(NodeCord[I3][1]-NodeCord[I1][1])
         - (NodeCord[I3][0]-NodeCord[I1][0])*(NodeCord[I2][1]-NodeCord[I1][1]);
 
    D = fabs( X*(RX-NodeCord[I1][0])+Y*(RY-NodeCord[I1][1])
          + Z*(RZ-NodeCord[I1][2]))/sqrt(X*X+Y*Y+Z*Z);
 
    for(L=0;L<=2;L++)
    {
        K = Node[L];
        if( L<2 )  K1 = Node[L+1];
        if( L==2 ) K1 = Node[0];
        R = sqrt( pow( (NodeCord[K1][0]-NodeCord[K][0]),  2.0)
                + pow( (NodeCord[K1][1]-NodeCord[K][1]),  2.0)
                + pow( (NodeCord[K1][2]-NodeCord[K][2]),  2.0) );

        HATL[L][0]=(NodeCord[K1][0]-NodeCord[K][0])/R;
        HATL[L][1]=(NodeCord[K1][1]-NodeCord[K][1])/R;
        HATL[L][2]=(NodeCord[K1][2]-NodeCord[K][2])/R;

        X=(RY-NodeCord[K][1])*(RZ-NodeCord[K1][2]) 
           -(RY-NodeCord[K1][1])*(RZ-NodeCord[K][2]);
        Y=(RX-NodeCord[K1][0])*(RZ-NodeCord[K][2]) 
           -(RX-NodeCord[K][0])*(RZ-NodeCord[K1][2]);
        Z=(RX-NodeCord[K][0])*(RY-NodeCord[K1][1]) 
           -(RX-NodeCord[K1][0])*(RY-NodeCord[K][1]);

        RO[L] = sqrt(X*X+Y*Y+Z*Z)/R;
        RHO[L]= sqrt(RO[L]*RO[L]-D*D);
        RJ[L] = sqrt( pow( (RX-NodeCord[K][0]), 2.0)
                    + pow( (RY-NodeCord[K][1]), 2.0)
                    + pow( (RZ-NodeCord[K][2]), 2.0) );
        R = (RX-NodeCord[K][0])*HATL[L][0]
          + (RY-NodeCord[K][1])*HATL[L][1]
          + (RZ-NodeCord[K][2])*HATL[L][2];
        SMLRO[L][0]=NodeCord[K][0]+R*HATL[L][0];
        SMLRO[L][1]=NodeCord[K][1]+R*HATL[L][1];
        SMLRO[L][2]=NodeCord[K][2]+R*HATL[L][2];
    }

    RIO = -2.0*M_PI*D;
    for(L=0;L<=2;L++)
    {
        L1=L+1;
        if(L == 2) L1=0;
        K=Node[L];

        if( L<2 ) K1=Node[L+1];
        if( L==2 ) K1=Node[0];
        X = HATL[L][0]*(NodeCord[K1][0]-SMLRO[L][0])
            + HATL[L][1]*(NodeCord[K1][1]-SMLRO[L][1])
            + HATL[L][2]*(NodeCord[K1][2]-SMLRO[L][2]);
        Y = HATL[L][0]*(NodeCord[K][0]-SMLRO[L][0])
            + HATL[L][1]*(NodeCord[K][1]-SMLRO[L][1])
            + HATL[L][2]*(NodeCord[K][2]-SMLRO[L][2]);
        RIO = RIO+(1.0/2.0)*RHO[L]*
              log((RJ[L1]+X)/(RJ[L1]-X)*(RJ[L]-Y)/(RJ[L]+Y) )
            + D*asin(D*X/RO[L]/sqrt(RJ[L1]*RJ[L1]-D*D))
            - D*asin(D*Y/RO[L]/sqrt(RJ[L] *RJ[L] -D*D));
    }

    RIO = RIO/(2.0*AREA);
 
    A =   pow( (NodeCord[I1][0]-NodeCord[I3][0]), 2.0 )
        + pow( (NodeCord[I1][1]-NodeCord[I3][1]), 2.0 )
        + pow( (NodeCord[I1][2]-NodeCord[I3][2]), 2.0 );
    B =   pow( (NodeCord[I2][0]-NodeCord[I3][0]), 2.0 )
        + pow( (NodeCord[I2][1]-NodeCord[I3][1]), 2.0 )
        + pow( (NodeCord[I2][2]-NodeCord[I3][2]), 2.0 );

    C=-2.0*( (RX-NodeCord[I3][0])*(NodeCord[I1][0]-NodeCord[I3][0])
         +(RY-NodeCord[I3][1])*(NodeCord[I1][1]-NodeCord[I3][1])
         +(RZ-NodeCord[I3][2])*(NodeCord[I1][2]-NodeCord[I3][2]));
    D=-2.0*( (RX-NodeCord[I3][0])*(NodeCord[I2][0]-NodeCord[I3][0])
         +(RY-NodeCord[I3][1])*(NodeCord[I2][1]-NodeCord[I3][1])
         +(RZ-NodeCord[I3][2])*(NodeCord[I2][2]-NodeCord[I3][2]));
    E= 2.0*( (NodeCord[I1][0]-NodeCord[I3][0])*(NodeCord[I2][0]-NodeCord[I3][0])
         +(NodeCord[I1][1]-NodeCord[I3][1])*(NodeCord[I2][1]-NodeCord[I3][1])
         +(NodeCord[I1][2]-NodeCord[I3][2])*(NodeCord[I2][2]-NodeCord[I3][2]));

    F=   pow((RX-NodeCord[I3][0]), 2.0)
       + pow((RY-NodeCord[I3][1]), 2.0)
       + pow((RZ-NodeCord[I3][2]), 2.0);

    ABE=sqrt(A+B-E);
    ACF=sqrt(A+C+F);
    BDF=sqrt(B+D+F);

    RJ1 = ( (2*B-C+D-E)*BDF+(2*A+C-D-E)*ACF )/4.0/(A+B-E)
         +(  4*(A+C)*(B+D+F)+4*F*(B-C-E)-(C+D+E)*(C+D+E) )/(8.0*(A+B-E)*ABE)
         *log(fabs( (2*ABE*BDF+(2*B-C+D-E))/(2*ABE*ACF-(2*A+C-D-E)) ));

    RJ2 = ( (2*B+D)*BDF-D*sqrt(F) )/(4.0*B)+(4*B*F-D*D)/(8*B*sqrt(B))
         *log(fabs( (2*sqrt(B)*BDF+2*B+D)/(2*sqrt(B*F)+D) ));
    RJ3= ( (2*A+C-D-E)*ACF+(2*B-C+D-E)*BDF )/(4.0*(A+B-E))
        +( 4*(A+C)*(B+D+F)+4*F*(B-C-E)-(C+D+E)*(C+D+E) )/(8.0*(A+B-E)*ABE)
        *log(fabs( (2*ABE*ACF+(2*A+C-D-E))/(2*ABE*BDF-(2*B-C+D-E)) ));
    RJ4 = ( (2*A+C)*ACF-C*sqrt(F) )/(4.0*A)+(4*A*F-C*C)
         /(8*A*sqrt(A))*log(fabs( (2*sqrt(A)*ACF+2*A+C)/(2*sqrt(A*F)+C) ));

    RD=4*A*B-E*E;
    RIP=( 4*B*(RJ1-RJ2)-2*E*(RJ3-RJ4)-(2*B*C-E*D)*RIO )/RD;
    RIQ=( 4*A*(RJ3-RJ4)-2*E*(RJ1-RJ2)-(2*A*D-E*C)*RIO )/RD;

    *RealIno = RIO;   
    *RealIxi = RIP;    
    *RealIeta = RIQ;

}    
    
    

/****************************************************************************
Prototype: void ComputeNonSingul(int *TrngleNode, double One,
                                 double *MidPoint, complex *Ino_PQ, 
                                 complex *Ixi_PQ,
                                 complex *Ieta_PQ, complex *Izeta_PQ) 
Description: To evaluate the MOM surface integral numerically when the
             observation triangle and the source triangle are different. 
Input value: 
    int *TrngleNode --- nodes of the source triangle 
    double One  --- a factor 
    double *MidPoint --- middle point of the observing trianlge 
Return value: 
    complex *Ino_PQ   --- to store the value of Ino
    complex *Ixi_PQ   --- to store the value of Ixi
    complex *Ieta_PQ  --- to store the value of Ieta
    complex *Izeta_PQ --- to store the value of Izeta.
    The above values are used in ComputeCMatrix().
Global value used: none 
Global value modified: none 
Subroutines called: COMplex_add(), COMplex_Cmplx(), COMplex_Expon(),
                    Real_Mul(), COMplex_Sub() 
****************************************************************************/
void ComputeNonSingul(int *TrngleNode, double One, double *MidPoint,
                      complex *Ino_PQ, complex *Ixi_PQ, complex *Ieta_PQ, 
                      complex *Izeta_PQ)
{
    double SrcPoint[3], Rcp;
    int QuadPointN;
    complex Eprj, V1;

    /* initialize the four complex variables to zero */

    for(QuadPointN=0;QuadPointN<TotQuadPoint;QuadPointN++)
    {
        ComputeGaussQuadPoint(QuadPointN,TrngleNode,SrcPoint);
 
        Rcp = VTXmag(MidPoint,SrcPoint);
        if( fabs(Rcp) <= 1.0e-07 )
        {

            Eprj = COMplex_Cmplx(0.0, -(Wght[QuadPointN]*WaveNumber) );
            *Ino_PQ = COMplex_Add(*Ino_PQ, Eprj);
            *Ixi_PQ = COMplex_Add(*Ixi_PQ, Real_Mul(Qpnt[QuadPointN][0],Eprj));
            *Ieta_PQ = COMplex_Add(*Ieta_PQ,
                                   Real_Mul(Qpnt[QuadPointN][1],Eprj) );

        }  /* if( fabs(Rcp */

        else
        {

            V1= COMplex_Cmplx(0.0,(Rcp*WaveNumber));
            V1= COMplex_Expon(-1.0,V1);
            V1 = COMplex_Cmplx( (Real(V1) - One),Aimag(V1));
            Eprj = Real_Mul( (Wght[QuadPointN]/Rcp),V1);
            *Ino_PQ = COMplex_Add(*Ino_PQ, Eprj);
            *Ixi_PQ = COMplex_Add(*Ixi_PQ, Real_Mul(Qpnt[QuadPointN][0],Eprj) );
            *Ieta_PQ = COMplex_Add(*Ieta_PQ, 
                                   Real_Mul(Qpnt[QuadPointN][1],Eprj) );

        }  /* else */

    }   /* for(QuadPointN */

    *Izeta_PQ = COMplex_Sub(*Ino_PQ,COMplex_Add(*Ixi_PQ , *Ieta_PQ) );

    
}  /* end of ComputeNonSingul( )  */
 


 
/**************************************************************************** 
Prototype:    void    ComputeNonSingul1(int *TrngleNode, double *MidPoint) 
Description:  To evaluate the MOM surface integral numerically when the
              observation triangle and the source triangle are 
              different. Results are stored in the global variables JxiPQ,  
              JetaPQ and JzetaPQ. 
Input value: 
    int *TrngleNode --- nodes of the source triangle 
    double *MidPoint --- middle point of the observing triangle 
Return value: none 
    Global value used: TotQuadPoint, Qpnt,  WaveNumber 
Global value modified: none
Subroutines called: ComputeGaussQuadPoint(), VTXmag(),
    COMplex_Cmplx(), Real_Mul(), COMplex_Mul(), COMplex_Add(),
    COMplex_Sub(), COMplex_Expon(), Aimag() 
*****************************************************************************/
void ComputeNonSingul1(int *TrngleNode, double *MidPoint, complex *JnoPQ,
                       complex *JxiPQ, complex *JetaPQ, complex *JzetaPQ )
{
    double SrcPoint[3], Rcp;
    int QuadPointN;

   
    for(QuadPointN=0;QuadPointN<TotQuadPoint;QuadPointN++)
    {
        ComputeGaussQuadPoint(QuadPointN,TrngleNode,SrcPoint);
        Rcp = VTXmag(MidPoint,SrcPoint);

        if( fabs(Rcp) <= 1.0e-07 )
        {
            fprintf(OutF1,"Check observation point \n"); 
            exit(1);
        }  /* if( fabs(Rcp) */
        else
        {
            complex Eprj, CFactor;

            CFactor = COMplex_Cmplx(0.0,(Rcp*WaveNumber));
            Eprj  = Real_Mul((Wght[QuadPointN]/(Rcp*Rcp*Rcp)),
                      COMplex_Cmplx((1.0+Real(CFactor)),Aimag(CFactor)));
            Eprj  = COMplex_Mul(Eprj,COMplex_Expon(-1.0,CFactor));
            *JnoPQ = COMplex_Add(*JnoPQ, Eprj);
            *JxiPQ = COMplex_Add(*JxiPQ,  Real_Mul(Qpnt[QuadPointN][0],Eprj) );
            *JetaPQ = COMplex_Add(*JetaPQ, Real_Mul(Qpnt[QuadPointN][1],Eprj) );

        }  /* else */

    }   /* for(QuadPoint */

    *JzetaPQ = COMplex_Sub(*JnoPQ,COMplex_Add(*JxiPQ, *JetaPQ) );

}
 


 
/****************************************************************************
Prototype:    void    ComputeParameters() 
Description:    To computes parameters(TetVolume, TrngleArea, TrngleCentroid
                and etc) which are necessary to build the FEM matrix 
                and MOM matrix, 
Input value:    none 
Return value:    none 
Global value used:     AbsPermeable, AbsPermitt, CFactor1, CFactor2,
    EdgeLength, FreeSpaceVel, ImpeDance, NodeCord, OperateFreq,
    GlobalEdgeEnds, TetVolume, TotEdgeNum, TotTetElement, TotTrnlgeNum,
    TParam, TrngleArea, TrngleCentroid, WaveLength, WaveNumber, 
Global value modified:     CFactor1, CFactor2, FreeSpaceVel,
    TetVolume, TrngleArea, TrngleCentroid 
Subroutines called:     COMplex_Cmplx(), ComputeTetraHeronVolume(),
                        ComputeTrngleCentroid(), Real_Mul(),  
                        TetfaceAreaNormal(), VTXmag() 
****************************************************************************/
void ComputeParameters()
{
    int II;

    OperateFreq   = OperateFreq * 1.0E+06;

    FreeSpaceVel  = 3.0E+08;
    AbsPermeable  = 1.25663706144E-06;
    AbsPermitt    = 8.8542E-12;

    WaveLength = FreeSpaceVel/OperateFreq;
    WaveNumber = 2.0*M_PI /WaveLength;
    ImpeDance  = WaveNumber/(2.0*M_PI*OperateFreq*AbsPermitt);

    CFactor1 = Real_Mul((-1.0/(8.0*M_PI)),
                        COMplex_Cmplx(0.0,(WaveNumber*ImpeDance)));
    CFactor2 = Real_Mul((-1.0/(2.0*M_PI)),
                        COMplex_Cmplx(0.0,-(ImpeDance/WaveNumber)));


    TParam[0]  = COMplex_Cmplx(0.0,-(1.0/(ImpeDance*WaveNumber)));
    TParam[1]  = COMplex_Cmplx(0.0,(WaveNumber/ImpeDance));
 
    TetVolume = double_Vector(TotTetElement);
    for(II=0;II<TotTetElement;II++)
        ComputeTetraHedronVolume(II,NodeCord,TetNode,TetVolume);
 
    EdgeLength = double_Vector(TotEdgeNum);
    for(II=0;II<TotEdgeNum;II++)  
        EdgeLength[II] = VTXmag(NodeCord[GlobalEdgeEnds[II][1]-1],
                         NodeCord[GlobalEdgeEnds[II][0]-1]);

    TrngleCenTroid = double_Matrix(TotTrngleNum,3);
    for(II=0;II<TotTrngleNum;II++)
        ComputeTrngleCentroid(II,NodeCord,TrngleNode,TrngleCenTroid);

    TrngleArea     = double_Vector(TotTrngleNum);
    TrngleNormal= double_Matrix(TotTrngleNum,3);
 
 
    for(II=0;II<TotTrngleNum;II++) 
        TetFaceAreaNormal(II,NodeCord,TrngleNode,TrngleNormal,TrngleArea);
 
}




/****************************************************************************
Prototype:  void  ComputeTetraHedronVolume(int TetHedNum, double**Cord, 
                                 int **TGNodeNum,double *TVolume ) 
Description:  To calculate the volume of a tetrahedron 
Input value: 
    int TetHedNum --- the number of  the tetrahedron 
    double **cord  --- the global node table 
    int **TGNodeNum --- the global node numbers of the tetrahedron 
    double *TVolume, -- where to store the results. 
Return value:     none 
Global value used:     none 
Global value modified:     none 
Subroutines called:     VTXsub1(), VTXcross() 
****************************************************************************/
void ComputeTetraHedronVolume(int TetHedNum, double **Cord, int **TGNodeNum,  
                              double *TVolume)
{
    int Count_i,Count_j;
    double Buff[4][3],Buff1[3],Buff2[3],Buff3[3],Buff4[3];

    for(Count_i=0;Count_i<=3;Count_i++)
        for(Count_j=0;Count_j<=2;Count_j++)
            Buff[Count_i][Count_j]=
                Cord[TGNodeNum[TetHedNum][Count_i]-1][Count_j];

    VTXsub1(1,0,Buff,Buff1);
    VTXsub1(2,0,Buff,Buff2);
    VTXsub1(3,0,Buff,Buff3);
    VTXcross(Buff1,Buff2,Buff4);

    TVolume[TetHedNum] =VTXdot(Buff3,Buff4);

}
 



/****************************************************************************
Prototype: void ComputeTrngleCentroid( int TrngleNum, double **Cord,
                                 int **BFaceNode, double **CenTroid) 
Description: To compute the centroid of a triangular patch. 
Input value: 
    int TrngleNum  --- the number of the triangle 
    double **Cord  --- the gloabl node table 
    int BFaceNode  --- the nodes of the triangle 
    double **Centroid -- where to store the results 
Return value:    none 
Global value used:    none 
Global value modified:    none 
Subroutines called:     VTXadd2() 
*****************************************************************************/
void ComputeTrngleCentroid(int TrngleNum, double **Cord, int **BFaceNode, 
                           double **CenTroid )
{
    int Count_k; double Buff[3],Buff1[3],Buff2[3],Buff3[3];
 
    for(Count_k=0;Count_k<=2;Count_k++)
    {
        Buff1[Count_k]=Cord[BFaceNode[TrngleNum][0]-1][Count_k];
        Buff2[Count_k]=Cord[BFaceNode[TrngleNum][1]-1][Count_k];
        Buff3[Count_k]=Cord[BFaceNode[TrngleNum][2]-1][Count_k];
    }

    VTXadd2(Buff1,Buff2,Buff3,Buff);
    CenTroid[TrngleNum][0] = (1.0/3.0) * Buff[0];
    CenTroid[TrngleNum][1] = (1.0/3.0) * Buff[1];
    CenTroid[TrngleNum][2] = (1.0/3.0) * Buff[2];

}
 






/****************************************************************************
Prototype: void TetFaceAreaNormal(int TrngleNum, double **Cord, int **BFaceNode,
                                  double **Normal, double *Area) 
Description: To compute the unit normal and the area of each face of the 
             tetrahedron 
Input value: 
    int TrngleNum --- index of the tetrahedron
    double **Cord --- the global node table
    int **BFaceNode  --- the nodes of each face
    double **Normal -- where to store the centroid of each face
    double *Area  --- where to store the area of each face
Return value:    none 
Global value used:     none 
Global value modified:     none 
Subroutines called:     VTXcross1(), VTXsub1(), VTXcross()
****************************************************************************/
void TetFaceAreaNormal(int TrngleNum, double **Cord, int **BFaceNode, 
                        double **Normal, double *Area)
{
    int j,k;
    double Buff[3][3],Buff1[3],Buff2[3],Buff3[3],MagAreaFace;
    double FaceNorm[3];

    for(j=0;j<=2;j++)
        for(k=0;k<=2;k++)
            Buff[j][k]= Cord[BFaceNode[TrngleNum][j]-1][k];
 
    VTXcross1(1,2,Buff,Buff1);
    VTXcross1(2,0,Buff,Buff2);
    VTXcross1(0,1,Buff,Buff3);

    FaceNorm[0] = (Buff1[0] + Buff2[0] + Buff3[0]);
    FaceNorm[1] = (Buff1[1] + Buff2[1] + Buff3[1]);
    FaceNorm[2] = (Buff1[2] + Buff2[2] + Buff3[2]);
    MagAreaFace = sqrt( FaceNorm[0]*FaceNorm[0] +
               FaceNorm[1]*FaceNorm[1] + FaceNorm[2]*FaceNorm[2]);
    Normal[TrngleNum][0] = FaceNorm[0] / MagAreaFace ;
    Normal[TrngleNum][1] = FaceNorm[1] / MagAreaFace ;
    Normal[TrngleNum][2] = FaceNorm[2] / MagAreaFace ;

    VTXsub1(0,1,Buff,Buff1);
    VTXsub1(2,1,Buff,Buff2);
    VTXcross(Buff1,Buff2,Buff3);
    Area[TrngleNum] = 0.5*sqrt(VTXdot(Buff3,Buff3));

}
  

  
/****************************************************************************
Prototype: void ComputeSourceVector(double WaveNum, 
                                    double*MidPoint, complex *EiC) 
Description:  To compute the E fields of an incident plane wave at the 
              boundary surface. 
Input value: 
    double WaveNum --- wavenumber of the plane wave. 
    double MidPoint --- midpoint of the observing triangle. 
    complex *EiC  --- where to store the results 
Return value:  none 
Global value used:  none 
Global value modified:  none 
Subroutines called:   none *****************************************************************************/
void ComputeSourceVector(double WaveNum, double *MidPoint, complex *EiC)
{
    double R,PolK[3],PolE[3];
    complex RE;

    /* the Progagation vector PolK[] in coordinate system */

    PolK[0] = sin(K_Theta*M_PI/180.0)*cos(K_Phi*M_PI/180.0);
    PolK[1] = sin(K_Theta*M_PI/180.0)*sin(K_Phi*M_PI/180.0);
    PolK[2] = cos(K_Theta*M_PI/180.0);

    /* the E polarization vector PolE[] in coordiante system */
    PolE[0] = sin(E_Theta*M_PI/180.0)*cos(E_Phi*M_PI/180.0);
    PolE[1] = sin(E_Theta*M_PI/180.0)*sin(E_Phi*M_PI/180.0);
    PolE[2] = cos(E_Theta*M_PI/180.0);

    /* RE=e(-jk*R)  vector productor, the phase term for every point */
    R = PolK[0]*MidPoint[0] + PolK[1]*MidPoint[1] + PolK[2]*MidPoint[2];
    RE= COMplex_Expon(-1.0,COMplex_Cmplx(0.0,(R*WaveNum)));
 
    /* E field including the phase term */
    EiC[0] = Real_Mul(PolE[0],RE);
    EiC[1] = Real_Mul(PolE[1],RE);
    EiC[2] = Real_Mul(PolE[2],RE);

}  /* end of ComputeSourceVector( ) */
 



/***************************************************************************
Prototype: void InvertMatrix(complex **Qmat, int MatrixSize) 
Description:    To invert a complex matrix. The results are stored in the 
                input matrix, thus no additional memory needed. 
Input value: 
    complex **Qmat --- a complex matrix 
    int matrixSize --- the size of the matrix 
Return value:     none 
Global value used:     none 
Global value modified:     none 
Subroutines called:    COMplex_Abs(), COMplex_Div(), COMplex_Cmplx(),
                       COMplex_Sub(), COMplex_Mul(), INT_Matrix(), 
                       free_INT_Matrix() 
****************************************************************************/
void InvertMatrix(complex **Qmat, int MatrixSize)
{
    int **Inter,i,j,k,l,KK,JJ,Row,Col;
    double Val1,Val2; complex tmp;

 
    Inter  = INT_Matrix(MatrixSize,2);

    for(i=0;i<MatrixSize;i++) 
        for(j=0;j<2;j++)  Inter[i][j]=0;
 
    for(k=0;k<MatrixSize;k++)  
    {
        JJ=k;
        if( k!=MatrixSize-1 )
        {
            KK= k+1; 
            Val1 = COMplex_Abs(Qmat[k][k]);
            for(i=KK;i<MatrixSize;i++)
            {
                Val2 = COMplex_Abs(Qmat[i][k]);
                if( (Val1 - Val2)<0.0 ) 
                {
                    Val1 = Val2; 
                    JJ=i;
                }     /* if( (Val1-Val2) */
            }         /* for(i=KK */
  
        }   /* if( k!=  */

        Inter[k][0] = k; Inter[k][1] = JJ;

        if( JJ!=k )
        {
            for(j=0;j<MatrixSize;j++)
            {
                tmp = Qmat[JJ][j]; 
                Qmat[JJ][j] = Qmat[k][j]; 
                Qmat[k][j] = tmp;

            }   /* for(j=0 */

            for(j=0;j<MatrixSize;j++)
                if(j!=k) Qmat[k][j]=COMplex_Div(Qmat[k][j],Qmat[k][k] );
            
        }    /* if( JJ!=K ) */
        else
        {
            for(j=0;j<MatrixSize;j++)
                if(j!=k)Qmat[k][j] = COMplex_Div( Qmat[k][j],Qmat[k][k] );
            
        }   /* else */


        Qmat[k][k] = COMplex_Div(COMplex_Cmplx(1.0,0.0),Qmat[k][k]);

        for(i=0;i<MatrixSize;i++)  
        {       
            if(i!=k)
            {
                for(j=0;j<MatrixSize;j++)
                    if( j!=k ) 
                    {      
                        Qmat[i][j]=COMplex_Sub(Qmat[i][j],
                        COMplex_Mul(Qmat[i][k],Qmat[k][j]));
                    }  /* if(j!=k) */              
            }  /* if( i!=k ) */
        }  /* for(i=0; */

        for(i=0;i<MatrixSize;i++)  
            if(i!=k) Qmat[i][k] = Real_Mul(-1.0, 
                          COMplex_Mul(Qmat[i][k],Qmat[k][k]));
           
   }   /* for(k=0;  */

   for(l=0;l<MatrixSize;l++)
   {
       k = MatrixSize - 1 - l;
       Row = Inter[k][0];
       Col = Inter[k][1];
       if( Row != Col ) 
           for(i=0;i<MatrixSize;i++)  
           {
               tmp = Qmat[i][Col]; 
               Qmat[i][Col]=Qmat[i][Row]; 
               Qmat[i][Row]=tmp;

           }   /* for(i=0 */
        
   }   /* for(l=0 */

   free_INT_Matrix(Inter,MatrixSize,2);

}    /*  end of InvertMatrix( )  */
 
 
 
/*****************************************************************************
Prototype: void ComputeCorrectionTerm() 
Description:  To couple the MoM matrix equation to the FEM matrix equation.
              The main idea is to obtain Jd  from C and D. 
              Then couple to the FEM matrix. 
Input value:     none 
Return value:     none 
Global value used:     Ccc, Ccd, Cdc, Cdd, Bdd, TotExtMetalEdgeNum,
                       HybrdBoundEdgeNum, Gd 
                       VcVertor, Vdvector
                       McVector,MdVetor
Global value modified:   none.  
Subroutines called:     InvertMatrix(), COMplex_Add() 
******************************************************************************/
void ComputeCorrectionTerm()
{
    int II,JJ,KK;
    complex CmpBuff, *TempBuff;


    TempBuff =  CMPLX_Vector(TotBoundEdgeNum);

    Null_CMPLXVector(UdVector,HybrdBoundEdgeNum);

    /*   Ccc <= Ccc^[-1]   */
    InvertMatrix(Ccc,TotExtMetalEdgeNum);
 
   
    for(II=0;II<TotExtMetalEdgeNum;II++) 
        McVector[II] = COMplex_Null();

    /*    Mc <=  Ccc * Vc ,  Mc =  Ccc^[-1] * Vc */
    for(JJ=0;JJ<TotExtMetalEdgeNum;JJ++)
    {
        CmpBuff = VcVector[JJ];
        for(II=0;II<TotExtMetalEdgeNum;II++)
            McVector[II] =COMplex_Add(McVector[II],  
                                COMplex_Mul(Ccc[II][JJ],CmpBuff));

    }   /* for(JJ=0 */

 
    for(KK=0;KK<HybrdBoundEdgeNum;KK++)
    {
        for(II=0;II<TotExtMetalEdgeNum;II++) 
            TempBuff[II] = COMplex_Null();

        for(JJ=0;JJ<TotExtMetalEdgeNum;JJ++)
        {
            CmpBuff = Ccd[JJ][KK];
            for(II=0;II<TotExtMetalEdgeNum;II++)
                TempBuff[II] =COMplex_Add(TempBuff[II],  
                                    COMplex_Mul(Ccc[II][JJ],CmpBuff));
        }  /* for(JJ=0 */

        for(II=0;II<TotExtMetalEdgeNum;II++) Ccd[II][KK] = TempBuff[II];

    }   /* for(KK=0 */

    for(KK=0;KK<HybrdBoundEdgeNum;KK++)
    {
        for(II=0;II<TotExtMetalEdgeNum;II++) 
            TempBuff[II] = COMplex_Null();

        for(JJ=0;JJ<TotExtMetalEdgeNum;JJ++)
        {
            CmpBuff = Dcd[JJ][KK];
            for(II=0;II<TotExtMetalEdgeNum;II++)
                TempBuff[II] = COMplex_Add(TempBuff[II], 
                                    COMplex_Mul(Ccc[II][JJ],CmpBuff));
        }  /* for(JJ=0 */

        for(II=0;II<TotExtMetalEdgeNum;II++) 
            Dcd[II][KK] = TempBuff[II];

    }   /* for(KK=0 */

    for(II=0;II<HybrdBoundEdgeNum;II++) 
        UdVector[II] = COMplex_Null();

    for(JJ=0;JJ<TotExtMetalEdgeNum;JJ++)
    {
        CmpBuff = McVector[JJ];

        for(II=0;II<HybrdBoundEdgeNum;II++)
            UdVector[II] = COMplex_Add(UdVector[II], 
                                COMplex_Mul(Cdc[II][JJ],CmpBuff));

    }   /* for(JJ=0 */

    /* Ud <=  Cdc * Mc - Vd */
    for(II=0;II<HybrdBoundEdgeNum;II++) 
        UdVector[II] = COMplex_Sub(UdVector[II],VdVector[II]);

    for(KK=0;KK<HybrdBoundEdgeNum;KK++)
    {
        for(II=0;II<HybrdBoundEdgeNum;II++) 
            TempBuff[II] = COMplex_Null();

        for(JJ=0;JJ<TotExtMetalEdgeNum;JJ++)
        {
            CmpBuff = Ccd[JJ][KK];
            for(II=0;II<HybrdBoundEdgeNum;II++)
                TempBuff[II] =COMplex_Add(TempBuff[II], 
                    COMplex_Mul(Cdc[II][JJ],CmpBuff));

        }  /* for(JJ=0 */

        for(II=0;II<HybrdBoundEdgeNum;II++)
            Cdd[II][KK] = COMplex_Sub(Cdd[II][KK],TempBuff[II]); /*Cdd=Cdd'*/

     }   /* for(KK=0 */

     for(KK=0;KK<HybrdBoundEdgeNum;KK++)
     {
         for(II=0;II<HybrdBoundEdgeNum;II++) 
             TempBuff[II] = COMplex_Null();

         for(JJ=0;JJ<TotExtMetalEdgeNum;JJ++)
         {
             CmpBuff = Dcd[JJ][KK];

             for(II=0;II<HybrdBoundEdgeNum;II++)
                 TempBuff[II] = COMplex_Add(TempBuff[II], 
                                    COMplex_Mul(Cdc[II][JJ],CmpBuff));

          }   /* for(JJ=0 */

          for(II=0;II<HybrdBoundEdgeNum;II++)
              Ddd[II][KK] = COMplex_Sub(Ddd[II][KK],TempBuff[II]); /*Ddd=Ddd' */

     }  /* for(KK=0  */


     /*  Cdd <= Cdd^[-1]   */
     InvertMatrix(Cdd,HybrdBoundEdgeNum);

     for(KK=0;KK<HybrdBoundEdgeNum;KK++)
     {
         for(II=0;II<HybrdBoundEdgeNum;II++) 
             TempBuff[II] = COMplex_Null();
 
         for(JJ=0;JJ<HybrdBoundEdgeNum;JJ++)
         {
             CmpBuff = Ddd[JJ][KK];
             for(II=0;II<HybrdBoundEdgeNum;II++)
                 TempBuff[II] = COMplex_Add( TempBuff[II], 
                       COMplex_Mul(Cdd[II][JJ],CmpBuff));
         }  /* for(JJ=0 */
         
         for(II=0;II<HybrdBoundEdgeNum;II++)
             Ddd[II][KK] = TempBuff[II]; /* Ddd=Cdd'^[-1]*Ddd' */

     }  /* for(KK=0 */

     for(II=0;II<HybrdBoundEdgeNum;II++) 
         MdVector[II] = COMplex_Null();

     /*  Md <= Cdd * Ud   */
     for(JJ=0;JJ<HybrdBoundEdgeNum;JJ++)
     {
         CmpBuff = UdVector[JJ];

         for(II=0;II<HybrdBoundEdgeNum;II++)
             MdVector[II] = COMplex_Add(MdVector[II], 
                                   COMplex_Mul(Cdd[II][JJ],CmpBuff));
     }   /* for(JJ=0 */

     for(KK=0;KK<HybrdBoundEdgeNum;KK++)
     {
         for(II=0;II<HybrdBoundEdgeNum;II++) 
             TempBuff[II] = COMplex_Null();

         for(JJ=0;JJ<HybrdBoundEdgeNum;JJ++)
         {
             CmpBuff = Ddd[JJ][KK];

             for(II=0;II<HybrdBoundEdgeNum;II++)
                 if( fabs(Bdd[II][JJ])>=1.0e-30 )
                     TempBuff[II] = COMplex_Add(TempBuff[II], 
                               Real_Mul(Bdd[II][JJ],CmpBuff));
                
         }   /* for(JJ=0  */

         for(II=0;II<HybrdBoundEdgeNum;II++)
             Cdd[II][KK] = Real_Mul(-1.0,TempBuff[II]);/* -Bdd*Cdd'^[-1]*Ddd'*/

     }  /* for(KK=0  */


     for(II=0;II<HybrdBoundEdgeNum;II++) 
         GdVector[II] = COMplex_Null();

     /*  Gd <= Bdd * Md  */
     for(JJ=0;JJ<HybrdBoundEdgeNum;JJ++)
     {
         CmpBuff = MdVector[JJ];
         for(II=0;II<HybrdBoundEdgeNum;II++)
             GdVector[II] = COMplex_Add(GdVector[II], 
                                 Real_Mul(Bdd[II][JJ],CmpBuff));
     }  /* for(JJ=0 */


     free_CMPLX_Vector(TempBuff, TotBoundEdgeNum);
      
 
}   /* end of ComputeCorrectionTerm */




/****************************************************************************
Prototype:  void  SurfaceFieldCompute() 
Description:  To compute the equivalent surface current density (Jc and Jd ) 
Input value:  none 
Return value:  none 
Global value used:  TotInnerEdgeNum EdVector HybrdBoundEdgeNum, 
                    JdVector, JcVector,TotExtMetalEdgeNum, Dcd, 
                    EdVector, 
                    MdVector,McVector
Global value modified:  none
Subroutines called:   COMplex_Add(),COMplex_Mul()
*****************************************************************************/

void SurfaceFieldCompute()
{
    int II,JJ;
    complex CmpBuff, *TempBuff;

    /*************************/
    for(II=TotInnerEdgeNum;II<TotMatrixEqnNum;II++)
        EdVector[II-TotInnerEdgeNum] = RHSVector[II];


    TempBuff =  CMPLX_Vector(TotBoundEdgeNum);
    for(II=0;II<HybrdBoundEdgeNum;II++) 
        TempBuff[II] = COMplex_Null();
    for(JJ=0;JJ<HybrdBoundEdgeNum;JJ++)
    {
        CmpBuff = EdVector[JJ];
        for(II=0;II<HybrdBoundEdgeNum;II++)
            TempBuff[II] = COMplex_Add( TempBuff[II], 
                     COMplex_Mul(Ddd[II][JJ],CmpBuff));
    }
    for(II=0;II<HybrdBoundEdgeNum;II++)
        JdVector[II] = COMplex_Add(TempBuff[II],MdVector[II]);
    /*  Jd = Ddd *Ed +Md */

    for(II=0;II<TotExtMetalEdgeNum;II++) 
        TempBuff[II] = COMplex_Null();
    for(JJ=0;JJ<HybrdBoundEdgeNum;JJ++)
    {
        CmpBuff = JdVector[JJ];
        for(II=0;II<TotExtMetalEdgeNum;II++)
        {
            TempBuff[II] = COMplex_Add(TempBuff[II], 
                 COMplex_Mul(Ccd[II][JJ],CmpBuff));
        }
    }
 

    for(II=0;II<TotExtMetalEdgeNum;II++)
        JcVector[II] = COMplex_Add(TempBuff[II],McVector[II]);
    for(II=0;II<TotExtMetalEdgeNum;II++)
        JcVector[II] = Real_Mul(-1.0,JcVector[II]);
   /*  Jc <= - Ccd * Jd - Mc, first assignment */

    for(II=0;II<TotExtMetalEdgeNum;II++) 
        TempBuff[II] = COMplex_Null();
    for(JJ=0;JJ<HybrdBoundEdgeNum;JJ++)
    {
        CmpBuff = EdVector[JJ];
        for(II=0;II<TotExtMetalEdgeNum;II++)
        {
            TempBuff[II] = COMplex_Add( TempBuff[II], 
                               COMplex_Mul(Dcd[II][JJ], CmpBuff));
        }
    }

    for(II=0;II<TotExtMetalEdgeNum;II++)
        JcVector[II] = COMplex_Add(TempBuff[II],JcVector[II]);
  /*  Jc <=  Dcd * Ed + Jc, second assignment  */

    free_CMPLX_Vector(TempBuff, TotBoundEdgeNum);


 }


/****************************************************************************
Prototype: void FemMatrixCompute(int NumOfTetElement, int MaxElementPerRow, 
                   int **TetEdge, int **GlobalEdgeEnds, 
                   double **NodeCord, complex *TParam, complex *Epsilon, 
                   int **GBLmatColAddr, complex **GBLmat, 
                   int *GBLmatColIndex, double *TetVolume) 
Description: To build the FEM matrix and store them using row-indexed scheme. 
Input value: 
    int NumOfTetElement --- total number of tetrahedron elements 
    int MaxElementPerRow --- the FEM matrix is a sparse matrix. Every row 
                             has a maximum number (MaxElementPerRow) of 
                             non-zero elemnet. Only non-zero elements are 
                             stored. 
    int **TetEdge --- the edge table of tetrahedron elements 
    int **GlobalEdgeEnds --- global edge table 
    double **NodeCord --- global node table 
    complex *TParam --- store constants used by FEM 
    complex *Epsilon --- the permittivity associated with each 
                          tetrahedron element. 
    int **GBLmatColAddr, complex **GBLmat, int *GBLmatColIndex --- the FEM   
                          matrix stored using row-indexed  scheme. 
    double *TetVolume ---the volume of tetrahedron elements 
Return value:     void 
Global value used:     none 
Global value modified: none 
Subroutines called:     VTXsub(), VTXcross(), Sign() 
*****************************************************************************/
void FemMatrixCompute(int NumOfTetElement, int MaxElementPerRow,
       int **TetEdge, int **GlobalEdgeEnds, double **NodeCord,     
       complex *TParam,complex *Epsilon, int **GBLmatColAddr, 
       complex **GBLmat,int *GBLmatColIndex, double *TetVolume )
{
    int    Count_i,Count_j,II,CountElmnt,RowNum,ColNum,IEdge,
           JEdge,INode,I5,QuadPointM,EdgeEnd1,EdgeEnd2,VVal;
    double  RRT[6][6],RIT[6][6],TetEdgeLen[6],TetNodeCord[4][3],Sgn,
            ObsPoint[3],TetWght[5],TetQuad[4][5],CurlTrm[6][3],wVector[6][3],
            gVector[6][3],fVector[6][3],VolumeDet,Buff[3],SgnM,SgnN,TmpVal;
    complex S[6][6],TFactor1,TFactor2;
 

    for(Count_i=0;Count_i<=3;Count_i++) 
        TetQuad[Count_i][0] = 1.0/4.0;
 
    TetWght[0] = -4.0/5.0;
 
    for(Count_j=0;Count_j<=3;Count_j++) 
    {
        for(Count_i=0;Count_i<=3;Count_i++)
            TetQuad[Count_i][Count_j+1] = 1.0/6.0;
     
        TetQuad[Count_j][Count_j+1] = 1.0/2.0;
        TetWght[Count_j+1] = 9.0/20.0;
     
    }  /* for(Count_j=0; */
 
    for(CountElmnt=0;CountElmnt<NumOfTetElement;CountElmnt++)
    {
        for(IEdge=0;IEdge<=5;IEdge++)
            for(JEdge=0;JEdge<=5;JEdge++)
            {
                RRT[IEdge][JEdge] = 0.0;
                RIT[IEdge][JEdge] = 0.0;
                S[IEdge][JEdge] = COMplex_Null();
            } /* for(JEdge=0; */
          
        for(IEdge=0;IEdge<=5;IEdge++) 
        {
            II = abs(TetEdge[CountElmnt][IEdge]) - 1;
            EdgeEnd1 = GlobalEdgeEnds[II][0] - 1;
            EdgeEnd2 = GlobalEdgeEnds[II][1] - 1;
            TetEdgeLen[IEdge]= VTXmag(NodeCord[EdgeEnd1],NodeCord[EdgeEnd2]);
        } /* for(IEdge=0; */
 
        for(INode=0;INode<=3;INode++) 
        {
            II = TetNode[CountElmnt][INode] - 1;
            TetNodeCord[INode][0] = NodeCord[II][0];
            TetNodeCord[INode][1] = NodeCord[II][1];
            TetNodeCord[INode][2] = NodeCord[II][2];
        } /* for(INode=0; */

        for(IEdge=0;IEdge<=5;IEdge++) 
        {
            I5 = abs(TetEdge[CountElmnt][5-IEdge]) - 1;
            Sgn = Sign(TetEdge[CountElmnt][5-IEdge]);
          
            if( Sgn>0.0) 
            {
                EdgeEnd1 = GlobalEdgeEnds[I5][0] - 1;
                EdgeEnd2 = GlobalEdgeEnds[I5][1] - 1;
            }  /* if(Sgn > 0.0) */
          
            else{
                EdgeEnd1 = GlobalEdgeEnds[I5][1] - 1;
                EdgeEnd2 = GlobalEdgeEnds[I5][0] - 1;
            }  /* else */
             
            VTXsub(NodeCord[EdgeEnd2],NodeCord[EdgeEnd1],gVector[IEdge]);
            VTXcross(NodeCord[EdgeEnd1],NodeCord[EdgeEnd2],fVector[IEdge]);
          
        }  /* for(IEdge=0; */
 
        VolumeDet = 6.0 * TetVolume[CountElmnt];
 
        for(QuadPointM=0;QuadPointM<=4;QuadPointM++) 
        {
      
            for(IEdge=0;IEdge<3;IEdge++)
                ObsPoint[IEdge] = 0.0;
          
            for(INode=0;INode<=3;INode++)  
            {
                ObsPoint[0] = ObsPoint[0] + 
                       TetQuad[INode][QuadPointM]*TetNodeCord[INode][0];
                ObsPoint[1] = ObsPoint[1] + 
                       TetQuad[INode][QuadPointM]*TetNodeCord[INode][1];
                ObsPoint[2] = ObsPoint[2] + 
                       TetQuad[INode][QuadPointM]*TetNodeCord[INode][2];
            }   /* for(INode=0; */
          
            for(IEdge=0;IEdge<=5;IEdge++)    
            {
                VTXcross(gVector[IEdge],ObsPoint,Buff);
                VTXadd(fVector[IEdge],Buff,wVector[IEdge]);
              
                CurlTrm[IEdge][0] = 2.0 * gVector[IEdge][0];
                CurlTrm[IEdge][1] = 2.0 * gVector[IEdge][1];
                CurlTrm[IEdge][2] = 2.0 * gVector[IEdge][2];
              
            }  /* for(IEdge=0; */
          
         for(IEdge=0;IEdge<=5;IEdge++) 
             for(JEdge=0;JEdge<=5;JEdge++)  
             {
                 RRT[IEdge][JEdge] = RRT[IEdge][JEdge] + TetWght[QuadPointM]
                       * VTXdot(CurlTrm[IEdge],CurlTrm[JEdge]) / VolumeDet ;
                 RIT[IEdge][JEdge] = RIT[IEdge][JEdge] + TetWght[QuadPointM]
                       * VTXdot(wVector[IEdge],wVector[JEdge])/VolumeDet;
                       
              }  /*  for(IEdge=0;  */
         }         /* for(QuadPointM=0; */ 
 
         TFactor1 = TParam[0];
         TFactor2 = COMplex_Mul(TParam[1],Epsilon[CountElmnt]);
       
         for(IEdge=0;IEdge<=5;IEdge++)
             for(JEdge=0;JEdge<=5;JEdge++)
             {
           
                 SgnM = Sign(TetEdge[CountElmnt][IEdge]);
                 SgnN = Sign(TetEdge[CountElmnt][JEdge]);
                 TmpVal = SgnM * TetEdgeLen[IEdge] * SgnN * TetEdgeLen[JEdge] ;
               
                 S[IEdge][JEdge] = COMplex_Add(
                        Real_Mul(RRT[IEdge][JEdge],TFactor1),
                        Real_Mul(RIT[IEdge][JEdge],TFactor2));

                 S[IEdge][JEdge] = Real_Mul(TmpVal,S[IEdge][JEdge]);
             }  /* for(JEdge=0; */

         for(IEdge=0;IEdge<=5;IEdge++) 
             for(JEdge=0;JEdge<=5;JEdge++)
             {

                 int i,flag;

                 /* find the Row number */
                 flag=0;

                 /* search the inner edge table */
                 for(i=0;i<TotInnerEdgeNum;i++) 		      
                     if( InnerEdgeStat[i] ==  
                         abs(TetEdge[CountElmnt][IEdge]))
                     { 
		         flag=1; 
                         RowNum=i;
                         break;
                     }

                 if( !flag )   /* not an inner edge number */
                     for(i=0;i<HybrdBoundEdgeNum;i++) /* search the boundary 
                                                            edge table */
		         if( BoundEdgeStat[i] ==  
                             abs(TetEdge[CountElmnt][IEdge]))
                         { 
		             flag=1; 
                             RowNum=i+TotInnerEdgeNum;
                             break;
                         }  
               
                 if(!flag) continue;

                 /* find the Column number */
                 flag=0;
                 for(i=0;i<TotInnerEdgeNum;i++) 
		   if( InnerEdgeStat[i] ==  
                       abs(TetEdge[CountElmnt][JEdge]) )
                   { 
		       flag=1; 
                       ColNum=i;
                       break;
                   }

                 if( !flag ) /* not an inner edge number */
                 for(i=0;i<HybrdBoundEdgeNum;i++) 
		       if( BoundEdgeStat[i] ==  
                           abs(TetEdge[CountElmnt][JEdge]))
                       { 
		           flag=1; 
                           ColNum=i+TotInnerEdgeNum;
                           break;
                       }
               
                 if(!flag) continue;
               
                 /* get one non-zero element */
                 if( COMplex_Abs(S[IEdge][JEdge]) >=1.0E-15 )
                 {
                     /* try to find whether it is already stored */
                     for(Count_i=0;Count_i<MaxElementPerRow;Count_i++)
                     {
                   
                         /* already stored in GBLmatColAddr */
                         if( GBLmatColAddr[RowNum][Count_i] == ColNum ) 
                         {
                             VVal=Count_i;
                             break;
                         }  /* if(GBLmatColAddr[RowNum][Count_i] == ColNum) */
                              
                         else VVal=-1; /* means "not available" in GBLmat */
                          
                     }   /* for(Count_i=0; */
                  
                  
                     /* a new member */        
                     if( VVal < 0)
                     {
                         GBLmatColAddr[RowNum][GBLmatColIndex[RowNum]]=ColNum;
                         GBLmat[RowNum][GBLmatColIndex[RowNum]]=S[IEdge][JEdge];
                         GBLmatColIndex[RowNum]++;
                     }  /*if(VVal < 0) */
                       
                     else GBLmat[RowNum][VVal] = 
                            COMplex_Add(GBLmat[RowNum][VVal],S[IEdge][JEdge]);
                              
                 }    /* if(COMplex_Abs(S[IEdge] */
                  
             }  /* for(JEdge=0; */
 
        } /*  for(CountElmnt=0; */
 
} /* end of the subroutine */



/***************************************************************************
Prototype: void  CreateRHSVector(int InnerEdgeNum, int BoundEdgeNum, 
                     int ForcdEdgeNum, complex *ISourceEdgeMag,
                     complex *RHSVector) 
Description: To use boundary conditions to create the right-hand side   
             vector. The final matrix equation is  [LHS][E]=[RHS]. After
             RHS is known,  [E] can be solved by using  ConjugateSolver() . 
Input value: 
    int InnerEdgeNum, -- number of total inner edges 
    int BoundEdgeNum, ---number of total surface edges 
    int ForcdEdgeNum---number of total forced edges
Return value: none 
Global value used:     SourceType ,GdVector, TotInnerEdgeNum 
Global value modified:     none 
Subroutines called:     COMplex_Mul(), COMplex_Add(), COMplex_Null(), 
                        Real_Mul() 
***************************************************************************/
void CreateRHSVector(int InnerEdgeNum, int BoundEdgeNum, int ForcdEdgeNum,
                     double *ISourceEdgeMag, complex *RHSVector)
{
 
    int IEdge, JEdge;

    if( SourceType=='I' )
    {
        for(IEdge=0; IEdge<TotISourceEdgeNum;IEdge++) 
        {

            /* find the index to RHSVector */
            for(JEdge=0; JEdge<TotInnerEdgeNum; JEdge++)
                if( InnerEdgeStat[JEdge]== ISourceEdgeStat[IEdge] ) break;

            if( JEdge==TotInnerEdgeNum )
            {
                fprintf(stderr, "Can't find the J source edge\n");
                exit(1);

            }   /* if( JEdge== ) */

            RHSVector[JEdge] = COMplex_Cmplx(                   
              EdgeLength[ISourceEdgeStat[IEdge]-1]*ISourceEdgeMag[IEdge],0.0 );

        }   /* for(IEdge */
    
    }       /*  if( !strcmp(SourceType  */


    /*************************************************/
    if( SourceType=='V' )
    {
        for(IEdge=0;IEdge<TotMatrixEqnNum;IEdge++)
        {
            RHSVector[IEdge] = COMplex_Null();
            if( IEdge>=TotInnerEdgeNum )
               RHSVector[IEdge] = GdVector[IEdge-TotInnerEdgeNum];
        }  /* for(IEdge */

    }   /* if(!strcmp(SourceType  */


    /************************************************/
    if( SourceType=='P')
    {  
        for(IEdge=0;IEdge<TotMatrixEqnNum;IEdge++)
        {
            RHSVector[IEdge] = COMplex_Null();
            if(IEdge>=TotInnerEdgeNum)
               RHSVector[IEdge] = GdVector[IEdge-TotInnerEdgeNum];

        }   /* for(IEdge */
    }       /* if( !strcmp */

}  /*  end of CreateRHSVector()  */



  
/****************************************************************************
Prototype:    double VectorNorm(complex *Vec, int Size)
Description:     To compute the the Euclidean norm of a complex vector
Input value: 
    complex *Vec  --- pointer to a vector of complex type
    int Size  ---  the length of the vector
Return value:    the Euclidean norm of the vector. If the vector is 
                 (x1, x2, ..., xn), the norm is defined as 
                 sqrt( |x1|^2 +|x2|^2+...+|xn|^2).
Global value used:     none 
Global value modified:     none 
Subroutines called:     COMplex_Abs().
*****************************************************************************/
double VectorNorm(complex *Vec, int Size)
{
    int II; 
    double sum=0,tmp=0;

  
    /* get the maxium component */
    for(II=0; II<Size;II++)
       if(COMplex_Abs(Vec[II])>tmp) tmp=COMplex_Abs(Vec[II]);

    /* return tmp; */

    for(II=0;II<Size;II++)
        sum += pow(COMplex_Abs(Real_Div(Vec[II],tmp)), 2.0);
    
    sum=tmp*sqrt(sum);

    return sum;
}

    

/******************************************************************************     
Prototype:     complex InnerProduct(complex *AVec,  complex *BVec, int Size)
Description:    To compute the inner product of  two vectors of complex type. 
Input value:     
    complex *AVector, *BVector --- two vectors of complex type, they should have the same length
    int Size  --- the length of the vectors.
Return value:   return the inner product of the two vectors.
                If B=(x1, x2, ... xn), A=(y1, y2, ..., yn),
                The inner product is defined as   x1y1* + x2y2* + ... +xnyn*, 
                where * denotes conjugate.
Global value used:  none
Global value modified:     none 
Subroutines called:   COMplex_Add().
*******************************************************************************/
complex InnerProduct(complex *AVec, complex *BVec, int Size)
{
    int II;
    complex I_P;

    I_P = COMplex_Null();
 
    for(II=0;II<Size;II++)
        I_P = COMplex_Add(I_P,COMplex_Mul(COMplex_Conjg(AVec[II]),BVec[II]));

    return I_P;
}


                  
/*******************************************************************************
Prototype:  void MatrixVectorProduct(char S, int Size, complex *XVec, 
              complex *AVec, int *RIIndex,int **RICol,complex **RIDat)
Description:    To multiply the final matrix equation { [FEM] +[Cdd] }
                with a vector. 
Input value: 
    complex **RIDat, int RIIndex, int **RICol --- the FEM matrix stored using the row-indexed scheme.
    complex *XVec --- pointer to a vector of complex type
    complex *AVec --- where the results are stored.
    int Size  --- the dimension of the vectors.  The vector and the matrix  
                  should have the same dimension.
    char S --- operation type. If s== '  ',  return [A]={[FEM] +[Cdd] }[X]. 
               If S='*', return [A]={[FEM] +[Cdd] }*[X], 
               where * denoteds conjugate.
Return value:   none
Global value used:  Cdd
Global value modified:     none 
Subroutines called:   COMplex_Add(), COMplex_Null(), COMplex_Mul(), 
                      COMplex_Conjg(), 
*******************************************************************************/
void MatrixVectorProduct(char S, int Size, complex *XVec, complex *AVec, 
                         int *RIIndex,int **RICol,complex **RIDat)
{   
    int II;
    int IEdge, JEdge;
  
    for(II=0;II<Size;II++)
        AVec[II]=COMplex_Null();
 
    if( S==' ' )
    {
        for(IEdge=0; IEdge<Size; IEdge++) 
            for(JEdge=0; JEdge<RIIndex[IEdge]; JEdge++)
                AVec[IEdge] = COMplex_Add(AVec[IEdge],
                      COMplex_Mul(RIDat[IEdge][JEdge],
                                  XVec[ RICol[IEdge][JEdge]]) );
     
        /* moment mehtod part is added to the right-bottom corner */
        for(IEdge=0;IEdge<HybrdBoundEdgeNum;IEdge++)  
            for(JEdge=0;JEdge<HybrdBoundEdgeNum;JEdge++)  
                 AVec[IEdge+TotInnerEdgeNum] =   COMplex_Add( 
                     AVec[IEdge+TotInnerEdgeNum],
                     COMplex_Mul(Cdd[IEdge][JEdge], XVec[JEdge+TotInnerEdgeNum]));
   }   /* if(S==' '  */
   
   
   else  /* S== '*'  */
   {
   
       for(IEdge=0; IEdge<Size; IEdge++) 
           for(JEdge=0; JEdge<RIIndex[IEdge]; JEdge++) 
               AVec[IEdge] = COMplex_Add(AVec[IEdge],
                   COMplex_Mul(COMplex_Conjg(RIDat[IEdge][JEdge]),
                               XVec[ RICol[IEdge][JEdge]]) );
      
       /* moment method part Cdd is added to the right-bottom corner */
       for(IEdge=0;IEdge<HybrdBoundEdgeNum;IEdge++)   
           for(JEdge=0;JEdge<HybrdBoundEdgeNum;JEdge++)               
               AVec[IEdge+TotInnerEdgeNum] = COMplex_Add( 
                    AVec[IEdge+TotInnerEdgeNum], 
                    COMplex_Mul(COMplex_Conjg(Cdd[JEdge][IEdge]), 
                                XVec[JEdge+TotInnerEdgeNum]));
      
   }  /*  else */
   
} /* end of MatrixVectorProduct() */
 


/*****************************************************************************
Prototype: void ConjugateSolver(int MatrixSize, int *RowIIndex, 
                    int **RowICol, complex **RowIDat, complex *RHSVec, 
                    int ITmax, double TOL) 
Description:  To solve a linear system using complex bi-conjugate 
                gradient method. 
Input value: 
    int MatrixSize--- the size of the matrix equation. 
    complex **RowIDat, int **RowICol, int *RowIIndex --- the FEM matrix
                    stored using the row-indexed scheme 
    complex *RHSVec----  the right-hand side vector, the results are 
                         stored in this vector. 
    int ITmax  --- the maximum iteration number. If the solver runs the maximum
                   iterations but still has not reach the tolerance, the 
                   solver will return and store the best results. 
    double TOL  --- the tolerance. If the solver reaches in the tolerance, 
                    it will stop and store the results. 
Return value:     none 
Global value used:     none 
Global value modified:     none 
Subroutines called:     none 
******************************************************************************/
void ConjugateSolver( int MatrixSize, int *RowIIndex,int **RowICol, 
             complex **RowIDat,complex *RHSVec,int ITmax, double TOL)
  
{   
    complex  *X, *P,*PP,*R,*RR,*A_P,*AH_PP,Alpha,Beta,temp_var, *tmp;
    int      II,IterIndex, Num, flag1, flag2, flag3, flag4, counter;
    double   Residue,Vnrm, Least;

    if(MatrixSize==0) return;   /* there is no coupling between MoM and FEM */


    X     = CMPLX_Vector(MatrixSize);  /* H */
    P     = CMPLX_Vector(MatrixSize);
    PP    = CMPLX_Vector(MatrixSize); /* P bar */
    R     = CMPLX_Vector(MatrixSize); /* R */
    RR    = CMPLX_Vector(MatrixSize); /* R bar */
    A_P   = CMPLX_Vector(MatrixSize); /* A*P */
    AH_PP = CMPLX_Vector(MatrixSize);  /* A^H *P bar */
    tmp   = CMPLX_Vector(MatrixSize); 
 
    for(II=0;II<MatrixSize;II++)
    {     
        X[II]     = COMplex_Cmplx(0.0,0.0);
        P[II]     = COMplex_Null();
        PP[II]    = COMplex_Null();
        R[II]     = COMplex_Null();
        RR[II]    = COMplex_Null();
        A_P[II]   = COMplex_Null();
        AH_PP[II] = COMplex_Null();
    }
         

    MatrixVectorProduct(' ',MatrixSize,X,R,RowIIndex,RowICol,RowIDat);
    
    for(II=0;II<MatrixSize;II++)
    {
        R[II]  = COMplex_Sub(RHSVec[II],R[II]);
        P[II]  = R[II];
        RR[II] = COMplex_Conjg(R[II]);  
        PP[II] = RR[II];
    }
   
 
    Vnrm = VectorNorm(RHSVec,MatrixSize);
 
    IterIndex = 0;

    Least = VectorNorm(R,MatrixSize)/Vnrm;
    for(II=0;II<MatrixSize;II++)  
        tmp[II]=X[II];


    counter=flag1=flag2=flag3=flag4=0;

 
    /* Actual Iteration */

    for(;;)
    {
   
        MatrixVectorProduct(' ',MatrixSize,P,A_P,RowIIndex,RowICol,RowIDat);
        MatrixVectorProduct('*',MatrixSize,PP,AH_PP,RowIIndex,RowICol,RowIDat);
  
        Alpha = InnerProduct(RR,R,MatrixSize);
   
        temp_var = InnerProduct(PP,A_P,MatrixSize);

        /* Alpha is the step length parameter */
        Alpha = COMplex_Div(Alpha,temp_var); 

        for(II=0;II<MatrixSize;II++)
        {
            /* New Solution Estimate */
            X[II]  = COMplex_Add(X[II],COMplex_Mul(Alpha,P[II]));
  
            /* Residual */
            R[II]  = COMplex_Sub(R[II],COMplex_Mul(Alpha,A_P[II])); 

            /* Bi-residual */
            RR[II] = COMplex_Sub(RR[II],
                    COMplex_Mul(COMplex_Conjg(Alpha),AH_PP[II]));       
        }
   
   
        Beta = InnerProduct(AH_PP,R,MatrixSize);
        Beta = COMplex_Div(Beta,temp_var);

        /* Beta is the Bi-Conjugacy coefficient  */
        Beta = Real_Mul(-1.0,Beta);      


        for(II=0;II<MatrixSize;II++)
        {
            /* Direction */
            P[II]  = COMplex_Add(R[II],COMplex_Mul(Beta,P[II]));   

            /*  Bi-Direction */
            PP[II] = COMplex_Add(RR[II], 
                         COMplex_Mul( COMplex_Conjg(Beta),  PP[II]));    
        }
    
        Residue = VectorNorm(R,MatrixSize)/Vnrm;
   
        if( Residue<Least )
        {
            Num= IterIndex;
            Least=Residue;
            for(II=0;II<MatrixSize;II++)  
                tmp[II]=X[II];
        }
   

        /* if it reaches convergence, break the loop */
        if( Residue<=TOL )  
        {
            /* Test termination condition */
            fflush(stdout);
            break;
        }

        /* if the error reaches within 0.2, restart the algorithm */
        if( (!flag1) &&(Residue<2.0e-1) )
        {
           
            /* re-initialize P,RR and PP */ 
            for(II=0;II<MatrixSize;II++) 
            {
                P[II]=R[II];
                RR[II]=COMplex_Conjg(R[II]);
                PP[II]=RR[II];
            }
            /* set flag1 to 1, thus it won't restart here later */
            flag1=1;
        } 

        /* if the error reaches within 0.05, restart the algorithm */
        if( (!flag2) &&(Residue<5.0e-2) )
        {
            /* re-initialize P,RR and PP */ 
            for(II=0;II<MatrixSize;II++) 
            {
                P[II]=R[II];
                RR[II]=COMplex_Conjg(R[II]);
                PP[II]=RR[II];
            }
            /* set flag2 to 1, thus it won't restart here later */
            flag2=1;      
        } 

        /* if the error reaches within 0.001, restart the algorithm */
        if( (!flag3) &&(Residue<1.0e-3) ){
           for(II=0;II<MatrixSize;II++) 
           {

               P[II]=R[II];
               RR[II]=COMplex_Conjg(R[II]);
               PP[II]=RR[II];
           }
       
           flag3=1;
       } 
   
   
       if( IterIndex==ITmax ) {    
           /* No. of iterations has exceeded the limit */
           printf("Iteration Exceeds\n");
           fflush(stdout);
           break;
       }
      
       IterIndex++;
   
   }/* End of for loop*/

   for(II=0;II<MatrixSize;II++) 
       RHSVec[II] = tmp[II];

   sprintf(log_mesg, "\nBi-Conjugate Gradient Solver Report:\n\n");
   sprintf(log_mesg+strlen(log_mesg), "\ttolerence:           %.2e\n", TOL);
   sprintf(log_mesg+strlen(log_mesg), "\toriginal RHS:        %.2e\n", Vnrm);
   sprintf(log_mesg+strlen(log_mesg), "\tthe least residue:   %.2e\n", Least);
   sprintf(log_mesg+strlen(log_mesg), "\tachieved when iterate %d times\n", 
                                      Num);
   sprintf(log_mesg+strlen(log_mesg), "\ttotal iteration number: %d\n",
           IterIndex );

}   /*  end of ConjugateSolver( ) */  
     





/****************************************************************************
Prototype:     void    ComputeCMatrix() 
Description:   Fill the Moment Method matrix C, partition it to 
               Ccc, Ccd, Cdd and  Cdc. 
Input value:   none 
Return value:  none 
Global value used:     TrngleNode,TrngleNormal,TotQuadPoint, TrngleCenTroid,
                       TrngleCenTroid, NodeCord, Ccc, Ccd,  Cdd, Cdc. 
Global value modified:     none 
Subroutines called:     ComputeGaussQuadPoint(),ComputeSingul(),  
                        COMplex_Cmplx(), ComputeNonSingul(), COMplex_Add2(),  
                        COMplex_Sub(), COMplex_Null() 
*****************************************************************************/
void ComputeCMatrix(int ObservTrngle, int SourceTrngle, 
                    double *ObsTrngleEdgeLen, double *SrcTrngleEdgeLen, 
                    int *ObservTrngleNode, int *SourceTrngleNode,
                    int *ObservTrngleEdge, int *SourceTrngleEdge,
                    double **ObsPointArray)
{
    int Count_i, Count_j, SourceEdge, ObservEdge, QuadPointM;
    int RowNum, MultRow,RowRootIndxNum,
        ColNum, MultCol,ColRootIndxNum;    
    int SrcNodeNum, SrcEdgeNum;
	double SgnM,  SgnN;
    double One, Factor1, MidPoint[3];
    complex V_x, V_y, V_z,  VPQ[3], CE[3][3],
            AVectorPQX[3], AVectorPQY[3], AVectorPQZ[3],
            InoPQ, IxiPQ, IetaPQ, IzetaPQ;
    double RealIno, RealIxi, RealIeta;


    One       = 0.0;
    RealIno   = 0.0;
    RealIxi   = 0.0;
    RealIeta  = 0.0;
 
    MidPoint[0] = TrngleCenTroid[ObservTrngle][0];
    MidPoint[1] = TrngleCenTroid[ObservTrngle][1];
    MidPoint[2] = TrngleCenTroid[ObservTrngle][2];
 
    if( ObservTrngle==SourceTrngle )
    {
        ComputeSingul(SourceTrngleNode,MidPoint,TrngleArea[SourceTrngle],
                      &RealIno, &RealIxi, &RealIeta);
        One = 1.0;
    }
 
    for(Count_i=0;Count_i<=2;Count_i++)
        for(Count_j=0;Count_j<=2;Count_j++)
            CE[Count_i][Count_j] = COMplex_Null();

    for(QuadPointM=0; QuadPointM<TotQuadPoint; QuadPointM++) 
    {
        InoPQ   = COMplex_Cmplx(RealIno, 0.0);
        IxiPQ   = COMplex_Cmplx(RealIxi, 0.0);
        IetaPQ  = COMplex_Cmplx(RealIeta,0.0);

        ComputeNonSingul(SourceTrngleNode,One,ObsPointArray[QuadPointM],
                          &InoPQ, &IxiPQ, &IetaPQ, &IzetaPQ);

        V_x = COMplex_Add2(Real_Mul( NodeCord[SourceTrngleNode[0]][0], IxiPQ),
                  Real_Mul( NodeCord[SourceTrngleNode[1]][0], IetaPQ),
                  Real_Mul( NodeCord[SourceTrngleNode[2]][0],IzetaPQ));
        V_y = COMplex_Add2(Real_Mul( NodeCord[SourceTrngleNode[0]][1], IxiPQ),
                  Real_Mul( NodeCord[SourceTrngleNode[1]][1], IetaPQ),
                  Real_Mul( NodeCord[SourceTrngleNode[2]][1],IzetaPQ));
        V_z = COMplex_Add2(Real_Mul( NodeCord[SourceTrngleNode[0]][2], IxiPQ),
                  Real_Mul( NodeCord[SourceTrngleNode[1]][2], IetaPQ),
                  Real_Mul( NodeCord[SourceTrngleNode[2]][2],IzetaPQ));
 
        for(SourceEdge=0;SourceEdge<=2;SourceEdge++)
        {
            SrcNodeNum = SourceTrngleNode[SourceEdge];
            SrcEdgeNum = abs(SourceTrngleEdge[SourceEdge]) - 1;

            AVectorPQX[SourceEdge] = 
                COMplex_Sub(V_x,Real_Mul(NodeCord[SrcNodeNum][0],InoPQ));
            AVectorPQY[SourceEdge] = 
                COMplex_Sub(V_y,Real_Mul(NodeCord[SrcNodeNum][1],InoPQ));
            AVectorPQZ[SourceEdge] = 
                COMplex_Sub(V_z,Real_Mul(NodeCord[SrcNodeNum][2],InoPQ));
            VPQ[SourceEdge] = InoPQ;

        }  /* for(SourceEdge=0 */


        for(ObservEdge=0;ObservEdge<=2;ObservEdge++)
        { 
            double  ObsQuadVector[3];

            /* ObservTrngleNode[] is indexed from 0 */              
            VTXsub(ObsPointArray[QuadPointM],  
                NodeCord[ObservTrngleNode[ObservEdge]],ObsQuadVector);

            for(SourceEdge=0;SourceEdge<=2;SourceEdge++)
            {
                complex V1, V2;

                V1 = COMplex_Add2( 
                          Real_Mul( ObsQuadVector[0],AVectorPQX[SourceEdge]), 
                          Real_Mul(ObsQuadVector[1], AVectorPQY[SourceEdge]), 
                          Real_Mul(ObsQuadVector[2],AVectorPQZ[SourceEdge]));
  
                V1 = COMplex_Mul(CFactor1,V1);
                V2 = COMplex_Mul(CFactor2,VPQ[SourceEdge]);
                SgnM = Sign(ObservTrngleEdge[ObservEdge]);
                SgnN    = Sign(SourceTrngleEdge[SourceEdge]);

                /* the weight is half of ordinary value, or Ali put 
                   a 2 factor into CFactor1 or CFactor2, here times 2 to 
                   get the correct answer */
                Factor1 = SgnN*SgnM * SrcTrngleEdgeLen[SourceEdge]*            
                    ObsTrngleEdgeLen[ObservEdge]*Wght[QuadPointM]*2;

                CE[ObservEdge][SourceEdge] = COMplex_Add( 
                     CE[ObservEdge][SourceEdge], 
                     Real_Mul(Factor1, COMplex_Add(V1,V2)));
            } /* for(SourceEdge=0; */
        }   /* for(ObservEdge=0; */

    }  /* for(QuadPointM=0 */

    for(ObservEdge=0;ObservEdge<=2;ObservEdge++)
    {
        /* RowNum actually refer to the edge number indexed to 
           GlobalEdgeEnds */
        /* it may be different from RowRootIndxNum */

        RowNum = abs(ObservTrngleEdge[ObservEdge]);
        SgnM   = Sign(ObservTrngleEdge[ObservEdge]);
        MultRow= GlobalEdgeEnds[RowNum-1][2];
        RowRootIndxNum= GlobalEdgeEnds[RowNum-1][3];

        if( MultRow<=0 ) continue;
        
        /* the following codes are for MulRow >0, that means this
           edge is inlcuded in C matrix, thus need to compute  */

        for(SourceEdge=0;SourceEdge<=2;SourceEdge++)
        {
            ColNum = abs(SourceTrngleEdge[SourceEdge]);
            SgnN = Sign(SourceTrngleEdge[SourceEdge]);
            MultCol= GlobalEdgeEnds[ColNum-1][2];
            ColRootIndxNum= GlobalEdgeEnds[ColNum-1][3];

            /* if the SourceEdge is not included in C matrix, jump to next
               SourceEdge */
            if( MultCol<=0 ) continue;
            
            /* the following codes are for SourceEdge which is included in 
               C matrix */
   
            /* another implementation of C[Row][Col] += CE[Ob][Sor] */
            if( ColRootIndxNum<HybrdBoundEdgeNum ) {
                if( RowRootIndxNum<HybrdBoundEdgeNum ) 
                    Cdd[RowRootIndxNum][ColRootIndxNum] = COMplex_Add(
                        Cdd[RowRootIndxNum][ColRootIndxNum],  
                        CE[ObservEdge][SourceEdge]);
                else Ccd[RowRootIndxNum-HybrdBoundEdgeNum][ColRootIndxNum] =  
                        COMplex_Add(
                          Ccd[RowRootIndxNum-HybrdBoundEdgeNum][ColRootIndxNum],
                          CE[ObservEdge][SourceEdge]);
        
            }   /* if(ColRootIndxNum<HybrdBoundEdgeNum) */

            else{
                /* ColRootIndxNum >= HybrdBoundEdgeNum */
           
                int col=ColRootIndxNum-HybrdBoundEdgeNum, 
                    row=RowRootIndxNum-HybrdBoundEdgeNum;
                
                if( RowRootIndxNum<HybrdBoundEdgeNum )
                    Cdc[RowRootIndxNum][col] = COMplex_Add(  
                        Cdc[RowRootIndxNum][col], 
                        CE[ObservEdge][SourceEdge]) ;
                else  Ccc[row][col] = COMplex_Add(Ccc[row][col], 
                                           CE[ObservEdge][SourceEdge]);
            }   /* else */
    
     
            /* whether a juction exist */
            for(Count_i=0; Count_i<JunctionNum; Count_i++)
            {
                int i;
        
                /* the index to C matrix of the external metal edge 
                   Junction{Count_i]-1 is the edge number. 
                   J vector is partitioned into Jd and Jc according to 
                   the edge status. Thus, the index to C  matrix is not 
                   necessary eqal to the edge number */

                i=GlobalEdgeEnds[ Junction[Count_i][1]-1 ][3];
        
                if( ObservTrngle==Junction[Count_i][2] && 
                    RowNum==Junction[Count_i][0] ) 
                {
        
                    /* C[i][ColRootIndxNum] += CE[ObservEdge][SourceEdge] */
           
                    if( ColRootIndxNum<HybrdBoundEdgeNum ) 
                    {
                        if( i<HybrdBoundEdgeNum ) 
                            Cdd[i][ColRootIndxNum] = COMplex_Add(
                                Cdd[i][ColRootIndxNum], 
                                CE[ObservEdge][SourceEdge]);
                              
                        else Ccd[i-HybrdBoundEdgeNum][ColRootIndxNum] =     
                                 COMplex_Add(  
                                     Ccd[i-HybrdBoundEdgeNum][ColRootIndxNum],
                                     CE[ObservEdge][SourceEdge]);
        
                    } /* if(ColRootIndxNum<HybrdBoundEdgeNum) */

                    else{
                        int col=ColRootIndxNum-HybrdBoundEdgeNum, 
                            row=i-HybrdBoundEdgeNum;
                
                        if( i<HybrdBoundEdgeNum )
                            Cdc[i][col]=COMplex_Add(Cdc[i][col], 
                                            CE[ObservEdge][SourceEdge]) ;
                        else  Ccc[row][col]=COMplex_Add(Ccc[row][col], 
                                                CE[ObservEdge][SourceEdge]);
                    }   /* else */
              
                }       /* if(ObservTrngle */

       
                if( SourceTrngle==Junction[Count_i][2] && 
                    ColNum==Junction[Count_i][0] ) 
                {
            
                    /* C[RowRootIndxNum][i] =+ CE[ObservEdge][SourceEdge] */
        
                    if( RowRootIndxNum<HybrdBoundEdgeNum ) {
                        if( i<HybrdBoundEdgeNum) 
                            Cdd[RowRootIndxNum][i] = COMplex_Add(
                              Cdd[RowRootIndxNum][i], 
                              CE[ObservEdge][SourceEdge]);
                              
                        else Cdc[RowRootIndxNum][i-HybrdBoundEdgeNum] 
                                =COMplex_Add(
                                     Cdc[RowRootIndxNum][i-HybrdBoundEdgeNum],
                                     CE[ObservEdge][SourceEdge]);
        
                    } /* if(RowRootIndxNum<HybrdBoundEdgeNum) */

                    else{  
                         /* RowRootIndxNum >= HybrdBoundEdgeNum */
                    
                         int col=i-HybrdBoundEdgeNum, 
                             row=RowRootIndxNum-HybrdBoundEdgeNum;
                
                         if( i<HybrdBoundEdgeNum )
                             Ccd[row][i] = COMplex_Add(
                                 Ccd[row][i], CE[ObservEdge][SourceEdge]) ;
                            else Ccc[row][col] = COMplex_Add(
                                 Ccc[row][col], CE[ObservEdge][SourceEdge]);
                    }   /* else */
              
                }        /* if(SourceTrngle */

       
                /* count junction on junction */
                for(Count_j=Count_i; Count_j<JunctionNum; Count_j++)
                {     
                    int flag1, flag2, j;
          
                    j=GlobalEdgeEnds[ Junction[Count_j][1]-1 ][3];
                    flag1=(SourceTrngle==Junction[Count_i][2] &&    
                           ColNum==Junction[Count_i][0] &&
                           ObservTrngle==Junction[Count_j][2] &&  
                           RowNum==Junction[Count_j][0]);
             
                    flag2=(SourceTrngle==Junction[Count_j][2] &&        
                           ColNum==Junction[Count_j][0] &&
                           ObservTrngle==Junction[Count_i][2] && 
                           RowNum==Junction[Count_i][0]);
               
           
                    /* should not count twice if flag1==flag2==1 */
                    if( flag1 ) 
                    {
           
                        /* C[i][j] =+ CE[ObservEdge][SourceEdge] */
              
                        if( j<HybrdBoundEdgeNum ) 
                        {
                            if( i<HybrdBoundEdgeNum ) 
                                Cdd[i][j] = COMplex_Add(Cdd[i][j], 
                                                CE[ObservEdge][SourceEdge]);
                            else Ccd[i-HybrdBoundEdgeNum][j] = COMplex_Add(
                                     Ccd[i-HybrdBoundEdgeNum][j], 
                                     CE[ObservEdge][SourceEdge]);
        
                        } /* if(ColRootIndxNum<HybrdBoundEdgeNum) */

                        else{
                            /* j >= HybrdBoundEdgeNum */
                            int row=i-HybrdBoundEdgeNum, 
                                col=j-HybrdBoundEdgeNum;
                
                            if( i<HybrdBoundEdgeNum )
                                Cdc[i][col] = COMplex_Add(Cdc[i][col],  
                                       CE[ObservEdge][SourceEdge]) ;
                            else  Ccc[row][col] = COMplex_Add(Ccc[row][col],  
                                       CE[ObservEdge][SourceEdge]);
                      
                        }   /* else  */
                    }        /* flag1 */

          
                    else  if(flag2)  
                    {     
               
                        /* C[j][i]=+ CE[ObservEdge][SourceEdge]) */
                     
                        if( i<HybrdBoundEdgeNum ) 
                        {
                            if( j<HybrdBoundEdgeNum ) 
                                Cdd[j][i] = COMplex_Add(Cdd[j][i], 
                                                CE[ObservEdge][SourceEdge]);
                            else Ccd[j-HybrdBoundEdgeNum][i] = COMplex_Add(
                                     Ccd[j-HybrdBoundEdgeNum][i], 
                                     CE[ObservEdge][SourceEdge]);
        
                        } /* if(ColRootIndxNum<HybrdBoundEdgeNum) */

                        else{  
                            /* i >= HybrdBoundEdgeNum */
                            int row=j-HybrdBoundEdgeNum,  
                                col=i-HybrdBoundEdgeNum;
                
                            if( j<HybrdBoundEdgeNum )
                                Cdc[i][col] = COMplex_Add(Cdc[i][col], 
                                           CE[ObservEdge][SourceEdge]) ;
                                                                              
                            else  Ccc[row][col] = COMplex_Add(Ccc[row][col],  
                                           CE[ObservEdge][SourceEdge]);
                        }   /* else  */
                    }       /* flag2 */

                    else continue;
                         
                } /* end for Count_j*/
       
            }    /* end for Count_i*/
        }   /* for(SourceEdge */
    }  /* for(ObservEdge  */
    
}   /* end of ComputeCMatrix()  */




/****************************************************************************
Prototype:   void  ComputeBddMatrix() 
Description:    Compute the B matrix generated by FEM. Only Bdd has non-zero   
                elements, thus only store Bdd. 
Input value:     none 
Return value:     none 
Global value used:  NodeCord, Junction, Bdd, 
                    GlobalEdgeEnds, TrngleCenTroid, 
Global value modified:      none  
Subroutines called:    sign(),  TrngleArea(), VTXdot()
*****************************************************************************/
void ComputeBddMatrix(int ObservTrngle, int SourceTrngle,
                      double *ObsTrngleEdgeLen, double *SrcTrngleEdgeLen, 
                      int *ObservTrngleNode, int *SourceTrngleNode,
                      int *ObservTrngleEdge, int *SourceTrngleEdge,
                      double *SrcTrngleNorm, double BE[][3])
{
    int Count_i, Count_j, ObservEdge, SourceEdge, NodeNum,
        RowNum, MultRow, RowRootIndxNum,
        ColNum, MultCol, ColRootIndxNum;

    double Value1, RFactor, SgnM, SgnN;

    for(Count_i=0;Count_i<=2;Count_i++)
        for(Count_j=0;Count_j<=2;Count_j++)
            BE[Count_i][Count_j] = 0.0;

    if( ObservTrngle==SourceTrngle ) 
    {
        double Buff[3], Buff0[3], Buff1[3], Buff2[3], Buff3[3], 
               Buff4[3], Buff5[3];

        RFactor = 4.0 * TrngleArea[ObservTrngle];
        
        Buff0[0] =TrngleCenTroid[ObservTrngle][0];
        Buff0[1] =TrngleCenTroid[ObservTrngle][1];
        Buff0[2] =TrngleCenTroid[ObservTrngle][2];
        
        for(ObservEdge=0;ObservEdge<=2;ObservEdge++)
        {
            SgnM = Sign(ObservTrngleEdge[ObservEdge]);
            NodeNum = ObservTrngleNode[ObservEdge];
           
            Buff1[0] = NodeCord[NodeNum][0];
            Buff1[1] = NodeCord[NodeNum][1];
            Buff1[2] = NodeCord[NodeNum][2];
           
            for(SourceEdge=0;SourceEdge<=2;SourceEdge++)
            {
                SgnN  = Sign(SourceTrngleEdge[SourceEdge]);
                NodeNum = SourceTrngleNode[SourceEdge];
               
                Buff2[0] = NodeCord[NodeNum][0];
                Buff2[1] = NodeCord[NodeNum][1];
                Buff2[2] = NodeCord[NodeNum][2];
 
                VTXcross(Buff1,Buff2,Buff3);
                VTXcross(Buff0,Buff1, Buff4);
                VTXcross(Buff2, Buff0,Buff5);
                VTXadd2(Buff3,Buff4,Buff5,Buff);
                Value1 = VTXdot(SrcTrngleNorm,Buff);
               
                BE[ObservEdge][SourceEdge] = ObsTrngleEdgeLen[ObservEdge]* 
                     SrcTrngleEdgeLen[SourceEdge]*(SgnM*SgnN/RFactor)*Value1;
                        
            }  /* for(SourceEdge=0; */
        }     /* for(ObservEdge=0  */
    }     /* if(ObservTrngle==  */

    for(ObservEdge=0; ObservEdge<=2;ObservEdge++){

        /* if the source is a metal edge, flag=0, other wise equal to 1 */
        int flag,i,j;
   
        RowNum = abs(ObservTrngleEdge[ObservEdge]);
        SgnM   = Sign(ObservTrngleEdge[ObservEdge]);
        MultRow= GlobalEdgeEnds[RowNum-1][2];
        RowRootIndxNum= GlobalEdgeEnds[RowNum-1][3];
   
        /* revised here on Feb 2 */
        /* Wn(r) should vanish in metal patch */
   
        /* to find whether the source triangle is a metal or diectrical patch */
        for(i=j=0;j<3;j++) 
            if( GlobalEdgeEnds[ abs(ObservTrngleEdge[j])-1][4]==2 ) i++;
      
        /* the source is a metal patch if metal edge number > 2*/
        /* if i==2, the patch has one edge locating on the hybrid contour */
       
        if( i>=2 && GlobalEdgeEnds[RowNum-1][4]==1 ) flag=0;
            else flag=1;

        /* don't need to compute this edge, thus jump to the next edge  */
        if( !(MultRow>0 && flag==1) ) continue;
       
        
        /* if MultRow <) && flag==1, continue the following codes */
        for(SourceEdge=0;SourceEdge<=2;SourceEdge++) {  
            ColNum = abs(SourceTrngleEdge[SourceEdge]);
            SgnN = Sign(SourceTrngleEdge[SourceEdge]);
            MultCol= GlobalEdgeEnds[ColNum-1][2];
            ColRootIndxNum= GlobalEdgeEnds[ColNum-1][3];
            
            /* don't need to compute this edge, thus jump to the next edge */
            if( MultCol<=0 )continue;

   
            /* if(MultCol>0), execute the following codes */

            if( RowRootIndxNum<HybrdBoundEdgeNum && 
                ColRootIndxNum <HybrdBoundEdgeNum )
            {
                Bdd[RowRootIndxNum][ColRootIndxNum] = 
                       Bdd[RowRootIndxNum][ColRootIndxNum] 
                     + BE[ObservEdge][SourceEdge];
            }

             /* whether a juction exist */
            for(Count_i=0; Count_i<JunctionNum; Count_i++){
                int i;
                
                /* the index to C matrix of the external metal edge 
                   Junction{Count_i]-1 is the edge number 
                   J vector is partitioned into Jd and Jc according to
                   the edge status. Thus, the index to C  matrix is not  
                   necessary eqal to the edge number */
                    
                i=GlobalEdgeEnds[ Junction[Count_i][1]-1 ][3];
        
                if( ObservTrngle==Junction[Count_i][2] && 
                    RowNum==Junction[Count_i][0] &&
                    i<HybrdBoundEdgeNum && 
                    ColRootIndxNum<HybrdBoundEdgeNum) 
                {
 
                    Bdd[i][ColRootIndxNum] = +BE[ObservEdge][SourceEdge];
                }
       
                if( SourceTrngle==Junction[Count_i][2] &&
                    ColNum==Junction[Count_i][0] && 
                    i<HybrdBoundEdgeNum  && 
                    RowRootIndxNum<HybrdBoundEdgeNum)
                {
                
                    Bdd[RowRootIndxNum][i] =+ BE[ObservEdge][SourceEdge];
                }
        
              
                /* count junction on junction */
                for(Count_j=Count_i; Count_j<JunctionNum; Count_j++){     
                    int flag1, flag2, j;
          
                    j=GlobalEdgeEnds[ Junction[Count_j][1]-1 ][3];
                    flag1=(SourceTrngle==Junction[Count_i][2] && 
                           ColNum==Junction[Count_i][0] &&
                           ObservTrngle==Junction[Count_j][2] && 
                           RowNum==Junction[Count_j][0]);
             
                    flag2=(SourceTrngle==Junction[Count_j][2] && 
                           ColNum==Junction[Count_j][0] &&
                           ObservTrngle==Junction[Count_i][2] && 
                           RowNum==Junction[Count_i][0]);
           
                    /* should not count twice if flag1==flag2==1 */
                    if( flag1 && i<HybrdBoundEdgeNum && j<HybrdBoundEdgeNum ) 
                        Bdd[i][j]=+BE[ObservEdge][SourceEdge];
                    else if( flag2 && i<HybrdBoundEdgeNum && 
                             j<HybrdBoundEdgeNum ) 
                             Bdd[j][i]=+BE[ObservEdge][SourceEdge];
                    else continue;
                                  
                } /*   for(Count_j=0     */
       
            }    /*   for(Count_i=0;    */

        }  /*   for(SourceEdge=0; */          
    }      /*   for(ObservEdge=0  */

}   /* end of ComputeBddMatrix() */





/*****************************************************************************
Prototype:  void  ComputeDMatrix() 
Description:    Fill D matrix generated by MOM. Partition it to Dcd and Ddd.
Input value:     none 
Return value:   none 
Global value used:   Junction,  NordCord, Ddd, Dcd, 
Global value modified:     none 
Subroutines called:   sign(),  Real_Mul(), computeNonSingul1(), 
                      COMplex_Add2(), COMplex_Null(),  
                      ComputeGaussQuadPoint(). 
*****************************************************************************/
void ComputeDMatrix(int  ObservTrngle, int SourceTrngle, 
                    double *ObsTrngleEdgeLen, double *SrcTrngleEdgeLen, 
                    int *ObservTrngleNode, int *SourceTrngleNode,
                    int *ObservTrngleEdge, int *SourceTrngleEdge,
                    double **ObsPointArray, double BE[][3]) 
{

    int Count_i, Count_j, ObservEdge, SourceEdge, flag, QuadPointM,
        RowNum, MultRow, RowRootIndxNum,
        ColNum, MultCol, ColRootIndxNum,
        SrcNodeNum, SrcEdgeNum;
    double Sgn, SgnM, SgnN, Factor1;
    complex V_x, V_y, V_z, DE[3][3],
            PPQX[3][13], PPQY[3][13], PPQZ[3][13],     
            JnoPQ, JxiPQ, JetaPQ, JzetaPQ;


    /* only need to compute Ddd and Dcd thus test source triangle 
       and observe triangle first */

    flag=0;  /* set to zero */
 
    for(SourceEdge=0;SourceEdge<=2;SourceEdge++)
    {
        ColNum = abs(SourceTrngleEdge[SourceEdge]);
        ColRootIndxNum= GlobalEdgeEnds[ColNum-1][3];
        if(ColRootIndxNum < HybrdBoundEdgeNum ) 
        {
            /* at least one edge needing to compute */
            flag=1;  
            break;
        }  /* if(ColRootIndxNum */
        
        /* one edge is a junction edge, it may need to compute */   
        for(Count_i=0; Count_i<JunctionNum; Count_i++)
            if( ColNum==Junction[Count_i][0]) 
            {
                flag=1;
                break;
            }  /* if( ColNum==Junction[Count_i][0]  */
    }   /* for(SourceEdge */
           
   
    /* only compute Ddd and Dcd */

    /* don't need to compute since the enforcement of boundary condition */
    if (!flag) return;
   
    
    if(ObservTrngle==SourceTrngle)
    {
        for(SourceEdge=0;SourceEdge<=2;SourceEdge++)
            for(QuadPointM=0;QuadPointM<TotQuadPoint;QuadPointM++)
            {
                PPQX[SourceEdge][QuadPointM] = COMplex_Null();
                PPQY[SourceEdge][QuadPointM] = COMplex_Null();
                PPQZ[SourceEdge][QuadPointM] = COMplex_Null();
            } /* for((QuadPointM=0 */
    }  /* if(ObservTrngle  */
   
    else 
    {
    
        for(QuadPointM=0;QuadPointM<=TotQuadPoint-1;QuadPointM++)
        {
            double ObsPoint[3];

            JnoPQ   = COMplex_Null();
            JxiPQ   = COMplex_Null();
            JetaPQ  = COMplex_Null();
            JzetaPQ = COMplex_Null();
 
            ObsPoint[0] = ObsPointArray[QuadPointM][0];
            ObsPoint[1] = ObsPointArray[QuadPointM][1];
            ObsPoint[2] = ObsPointArray[QuadPointM][2];
 
            ComputeNonSingul1(SourceTrngleNode,ObsPoint, &JnoPQ, &JxiPQ,  
                              &JetaPQ, &JzetaPQ);

 
            V_x = COMplex_Add2(
                      Real_Mul( NodeCord[SourceTrngleNode[0]][0],  JxiPQ),
                      Real_Mul( NodeCord[SourceTrngleNode[1]][0], JetaPQ),
                      Real_Mul( NodeCord[SourceTrngleNode[2]][0],JzetaPQ));

            V_y = COMplex_Add2(
                      Real_Mul( NodeCord[SourceTrngleNode[0]][1],  JxiPQ),
                      Real_Mul( NodeCord[SourceTrngleNode[1]][1], JetaPQ),
                      Real_Mul( NodeCord[SourceTrngleNode[2]][1],JzetaPQ));

            V_z = COMplex_Add2(Real_Mul( 
                      NodeCord[SourceTrngleNode[0]][2],  JxiPQ),
                      Real_Mul( NodeCord[SourceTrngleNode[1]][2], JetaPQ),
                      Real_Mul( NodeCord[SourceTrngleNode[2]][2],JzetaPQ));

            for(SourceEdge=0;SourceEdge<=2;SourceEdge++) 
            {
                double Buff[3];

                SrcNodeNum = SourceTrngleNode[SourceEdge];
                Sgn        = Sign(SourceTrngleEdge[SourceEdge]);
                SrcEdgeNum = abs(SourceTrngleEdge[SourceEdge]) - 1;
                Factor1    = Sgn * SrcTrngleEdgeLen[SourceEdge]/(4.0*M_PI);
 
         	Buff[0]    = NodeCord[SrcNodeNum][0];
         	Buff[1]    = NodeCord[SrcNodeNum][1];
         	Buff[2]    = NodeCord[SrcNodeNum][2];

         	PPQX[SourceEdge][QuadPointM]  = Real_Mul( Factor1, 
                    COMplex_Add(
                       Real_Mul((ObsPoint[1]*Buff[2] - ObsPoint[2]*Buff[1]), 
                                JnoPQ),
                    COMplex_Sub(   Real_Mul((Buff[1] - ObsPoint[1]) ,V_z),
                       Real_Mul((Buff[2] - ObsPoint[2]) ,V_y))));

         	PPQY[SourceEdge][QuadPointM ] = Real_Mul( Factor1, 
                    COMplex_Add(
                       Real_Mul((ObsPoint[2]*Buff[0] - ObsPoint[0]*Buff[2]), 
                                JnoPQ),
                    COMplex_Sub(   Real_Mul((Buff[2] - ObsPoint[2]) ,V_x),
                       Real_Mul((Buff[0] - ObsPoint[0]) ,V_z))));

         	PPQZ[SourceEdge][QuadPointM]  = Real_Mul( Factor1, 
                    COMplex_Add(
                        Real_Mul((ObsPoint[0]*Buff[1] - ObsPoint[1]*Buff[0]), 
                                 JnoPQ),
                    COMplex_Sub(   Real_Mul((Buff[0] - ObsPoint[0]) ,V_y),
                         Real_Mul((Buff[1] - ObsPoint[1]) ,V_x))));

    	    } /* for(SourceEdge=0 */
        }   /* for(QuadPointM=0; */
    } /* else  */
 
    for(Count_i=0;Count_i<=2;Count_i++)
        for(Count_j=0;Count_j<=2;Count_j++)
            DE[Count_i][Count_j] = COMplex_Null();
 
    for(ObservEdge=0;ObservEdge<=2;ObservEdge++)
    {
   	Sgn = Sign(ObservTrngleEdge[ObservEdge]);
   	for(SourceEdge=0;SourceEdge<=2;SourceEdge++)
        {
            complex V1;

   	    V1 = COMplex_Null();
   	    for(QuadPointM=0;QuadPointM<=TotQuadPoint-1;QuadPointM++)
                ComputeGaussQuadPoint(QuadPointM, ObservTrngleNode, 
                                      ObsPointArray[QuadPointM]);
 
   	    for(QuadPointM=0;QuadPointM<=TotQuadPoint-1;QuadPointM++)
            {
                double Buff[3];

                Buff[0] = (ObsPointArray[QuadPointM][0]
                            - NodeCord[ObservTrngleNode[ObservEdge]][0]);
                Buff[1] = (ObsPointArray[QuadPointM][1]
                            - NodeCord[ObservTrngleNode[ObservEdge]][1]);
                Buff[2] = (ObsPointArray[QuadPointM][2]
                            - NodeCord[ObservTrngleNode[ObservEdge]][2]);
                V1 = COMplex_Add(V1, Real_Mul( Wght[QuadPointM] ,
                         COMplex_Add2(
                             Real_Mul(Buff[0],PPQX[SourceEdge][QuadPointM]),
                             Real_Mul(Buff[1],PPQY[SourceEdge][QuadPointM]),
                             Real_Mul(Buff[2],PPQZ[SourceEdge][QuadPointM]))));

            } /* for(QuadPointM=0  */

            V1  = Real_Mul( (Sgn*ObsTrngleEdgeLen[ObservEdge]), V1);
            DE[ObservEdge][SourceEdge] = COMplex_Sub(V1,
                COMplex_Cmplx((BE[ObservEdge][SourceEdge]/2.0),0.0));

        }  /* for(SourceEdge=0  */

    }/* for(ObservEdge=0;  */
   
    for(ObservEdge=0;ObservEdge<=2;ObservEdge++)
    {
        RowNum = abs(ObservTrngleEdge[ObservEdge]);
        SgnM   = Sign(ObservTrngleEdge[ObservEdge]);
        MultRow= GlobalEdgeEnds[RowNum-1][2];
        RowRootIndxNum= GlobalEdgeEnds[RowNum-1][3];
        
        /* don't need to compute this observing edge */
        if(MultRow <=0 ) continue;
       
        /* the following codes are used since this observing edge is
           included in D matrix */
          
        for(SourceEdge=0;SourceEdge<=2;SourceEdge++) 
        {
   
            /* if the source is a metal edge, flag=0, other wise equal to 1 */
            int flag,i,j;
    
            ColNum = abs(SourceTrngleEdge[SourceEdge]);
            SgnN = Sign(SourceTrngleEdge[SourceEdge]);
            MultCol= GlobalEdgeEnds[ColNum-1][2];
            ColRootIndxNum= GlobalEdgeEnds[ColNum-1][3];
   
            /* revised here on Feb 2 */
            /* Wn(r) should vanish in metal patch */
   
            /* to find whether the source triangle is a metal 
               or diectrical patch */
            for(i=j=0;j<3;j++) 
                if( GlobalEdgeEnds[ abs(SourceTrngleEdge[j])-1][4]==2)i++;
      
            /* the source is a metal patch, Wn(r) should vanish here */
            if( i>=2 && GlobalEdgeEnds[ColNum-1][4]==1 )
                flag=0;
            else flag=1;
       
            if( MultCol>0 && flag==1 )
            {
                   
                if( ColRootIndxNum<HybrdBoundEdgeNum  ) 
                {     
                    if( RowRootIndxNum<HybrdBoundEdgeNum ) 
                        Ddd[RowRootIndxNum][ColRootIndxNum] = COMplex_Add(
                            Ddd[RowRootIndxNum][ColRootIndxNum],
                            DE[ObservEdge][SourceEdge]);
                                 
                    else Dcd[RowRootIndxNum-HybrdBoundEdgeNum][ColRootIndxNum]= COMplex_Add(
                        Dcd[RowRootIndxNum-HybrdBoundEdgeNum][ColRootIndxNum],
                                          DE[ObservEdge][SourceEdge]);
                                          
                }  /* if( ColRootIndxNum< */

                         
                /* whether a juction exist */
                for(Count_i=0; Count_i<JunctionNum; Count_i++)
                {
                   
                    int i,j;
        
                    /* the index to C matrix of the external metal edge */
                    /* Junction{Count_i]-1 is the edge number */
                    /* J vector is partitioned into Jd and Jc according
                       to the edge status */
                    /* thus, the index to C  matrix is not necessary eqal 
                       to the edge number */
                       
                    i=GlobalEdgeEnds[ Junction[Count_i][1]-1 ][3];
        
                    if( ObservTrngle==Junction[Count_i][2] && 
                        RowNum==Junction[Count_i][0] && 
                        ColRootIndxNum<HybrdBoundEdgeNum ) 
                    {
                        
                        if( i<HybrdBoundEdgeNum )
                            Ddd[i][ColRootIndxNum] = COMplex_Add(
                                 Ddd[i][ColRootIndxNum],
                                 DE[ObservEdge][SourceEdge]);
                        else
                        {  
                            Dcd[i-HybrdBoundEdgeNum][ColRootIndxNum]= COMplex_Add(             
                                     Dcd[i-HybrdBoundEdgeNum][ColRootIndxNum],  
                                     DE[ObservEdge][SourceEdge]);  
                        }  /* else */
                    }      /* if(ObservTrngle */
       
                    if( SourceTrngle==Junction[Count_i][2] &&  
                        ColNum==Junction[Count_i][0] && 
                        i< HybrdBoundEdgeNum) 
                    {
                        if( RowRootIndxNum<HybrdBoundEdgeNum ) 
                            Ddd[RowRootIndxNum][i] = COMplex_Add( 
                                Ddd[RowRootIndxNum][i],  
                                DE[ObservEdge][SourceEdge]);
                        else  
                            Dcd[RowRootIndxNum-HybrdBoundEdgeNum][i]= COMplex_Add(
                                    Dcd[RowRootIndxNum-HybrdBoundEdgeNum][i],
                                    DE[ObservEdge][SourceEdge]);
                    } /* if(SourceTrngle  */
        
                        
                    /* count junction on junction */
                    for(Count_j=Count_i; Count_j<JunctionNum; Count_j++)
                    { 
                        int flag1, flag2;
          
                        j=GlobalEdgeEnds[ Junction[Count_j][1]-1 ][3];
                        flag1=(SourceTrngle==Junction[Count_i][2] && 
                               ColNum==Junction[Count_i][0] &&
                               ObservTrngle==Junction[Count_j][2] && 
                               RowNum==Junction[Count_j][0]);
              
                        flag2=(SourceTrngle==Junction[Count_j][2] && 
                               ColNum==Junction[Count_j][0] &&
                               ObservTrngle==Junction[Count_i][2] && 
                               RowNum==Junction[Count_i][0]);
                
                        /* should not count twice if flag1==flag2==1 */
                        if( flag1 && j<HybrdBoundEdgeNum ) 
                        {   
                            if( i<HybrdBoundEdgeNum )
                                Ddd[i][j] = COMplex_Add( Ddd[i][j],  
                                DE[ObservEdge][SourceEdge]);
                            else 
                                Dcd[i-HybrdBoundEdgeNum][j]=
                                   COMplex_Add( Dcd[i-HybrdBoundEdgeNum][j],  
                                                DE[ObservEdge][SourceEdge]);
                        }  /* if( flag1 )  */
                                     
                        else  if( flag2 && i<HybrdBoundEdgeNum ) 
                              {
                                   if( j<HybrdBoundEdgeNum )
                                       Ddd[j][i] = COMplex_Add( Ddd[j][i],  
                                           DE[ObservEdge][SourceEdge]);
                                   else 
                                       Dcd[j-HybrdBoundEdgeNum][i]=COMplex_Add(
                                           Dcd[j-HybrdBoundEdgeNum][i],  
                                           DE[ObservEdge][SourceEdge]);
                               }  /*  if(flag2)     */
                        else continue;
                    } /* for(Count_j=Count_i */
		       
                }    /* for(Count_i=0;*/

            } /* if(MultCol>0 && flag==1) */
        }      /* for(SourceEdge=0;  */
    }            /* for(ObservEdge=0;   */
}   /* end of ComputeDMatrix() */





/****************************************************************************
Prototype:     double    Sign( int Value) 
Description:     To return the sign bit of an integer 
Input value: 
    int Value --- the input integer 
Return value: if Value>=0, return 1, otherwise return -1. 
*****************************************************************************/
double Sign(int Value)
{
 return (Value>0 ? 1:-1);
}


 

/****************************************************************************
Prototype:    void    ReadInputFile() 
Description:    To read data from an input  file, initialize global variables. Input value:     none 
Return value:     none 
Global value used: TotEdgeNum, TotTetElement, TotNodeNum, 
                   TotTrngleNum,TotBoundEdgeNum, TotInnerEdgeNum,      
                   DielBoundEdgeNum, HybrdBoundEdgeNum, TrngleNode, 
                   TrngleEdge,NodeCord, TetNode, TetNode, 
                   Epsilon,TetEdge , GlobalEdgeEnds, PlusTrngleDetect, 
                   PlusTrngleIndex,MinusTrngleDetect,MinusTrngleIndex, 
                   InnerEdgeStat, BoundEdgeStat,SourceType,VSourceMag, 
                   VSourceNum, OperateFreq,
                   VSourceEdge 
Global value modified:  the save as "Global value used".
Subroutines called:     INT_Matrix(), INT_Vector() 
*****************************************************************************/
void ReadInputFile()
{
    char buff[100]; /* a buffer used to store comment information */
    int i,j, IEdge;
    char c;
    /* read the global information table */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your global information table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
 
    fscanf(InF,"%d",&TotEdgeNum);
    fscanf(InF,"%d",&TotTetElement);
    fscanf(InF,"%d",&TotNodeNum);
    fscanf(InF,"%d",&TotTrngleNum);
 
    fscanf(InF,"%d",&TotBoundEdgeNum);
    fscanf(InF,"%d",&TotInnerEdgeNum);
    fscanf(InF,"%d",&DielBoundEdgeNum);
    fscanf(InF,"%d",&HybrdBoundEdgeNum);
    fscanf(InF,"%d",&TotISourceEdgeNum);

    /* read the global node table */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your global node table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    NodeCord = double_Matrix(TotNodeNum,3);
    for(i=0;i<TotNodeNum;i++)      
        for(j=0;j<=2;j++)
            fscanf(InF,"%lf",&NodeCord[i][j]);

    /* read the global edge table */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your global edge table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    GlobalEdgeEnds= INT_Matrix(TotEdgeNum,6);
    for(i=0;i<TotEdgeNum;i++)
        for(j=0;j<6;j++)
        {
            fscanf(InF," %d",&GlobalEdgeEnds[i][j]);
        }
            

    /* read the triangle table */
    fscanf(InF, " %c", &c);
    if ( c!='#' ) 
    {
       fprintf(stderr, "Check your triangle table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    TrngleNode = INT_Matrix(TotTrngleNum,3);
    TrngleEdge = INT_Matrix(TotTrngleNum,3);
    for(i=0;i<TotTrngleNum;i++) 
    {
        for(j=0;j<=2;j++) fscanf(InF,"%d",&TrngleNode[i][j]);
        for(j=0;j<=2;j++) fscanf(InF,"%d",&TrngleEdge[i][j]);
    }
  
  

    /* read tetrahedron table */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your tetrahedron table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    TetNode = INT_Matrix(TotTetElement,4);
    Epsilon = CMPLX_Vector(TotTetElement);
    for(i=0;i<TotTetElement;i++)
    {
        for(j=0;j<=3;j++)
            fscanf(InF,"%d",&TetNode[i][j]);
        fscanf(InF,"%lf",&Epsilon[i].x);
        fscanf(InF,"%lf",&Epsilon[i].y);
    }

    TetEdge = INT_Matrix(TotTetElement,6);
    for(i=0;i<TotTetElement;i++) 
        for(j=0;j<=5;j++)   
            fscanf(InF,"%d",&TetEdge[i][j]);


    /* read plus/minus triangle table */    
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your plus/minus triangle table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    PlusTrngleDetect = INT_Matrix(TotEdgeNum,3);
    PlusTrngleIndex = INT_Vector(TotEdgeNum);
    MinusTrngleDetect = INT_Matrix(TotEdgeNum,3);
    MinusTrngleIndex = INT_Vector(TotEdgeNum);
    for(i=0;i<TotEdgeNum;i++)
    {
       fscanf(InF,"%d",&PlusTrngleIndex[i]);
       for(j=0;j<PlusTrngleIndex[i];j++)
           fscanf(InF,"%d",&PlusTrngleDetect[i][j]);

       fscanf(InF,"%d",&MinusTrngleIndex[i]);
       for(j=0;j<MinusTrngleIndex[i];j++)
           fscanf(InF,"%d",&MinusTrngleDetect[i][j]);
    }

    /* read inner edge table */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your inner edge table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    InnerEdgeStat = INT_Vector(TotInnerEdgeNum);
    for(i=0;i<TotInnerEdgeNum;i++) 
        fscanf(InF,"%d",&InnerEdgeStat[i]);

    /* read boundary edge table */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your boundary edge table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    BoundEdgeStat = INT_Vector(TotBoundEdgeNum);
    for(i=0;i<TotBoundEdgeNum;i++)
        fscanf(InF,"%d",&BoundEdgeStat[i]);


    /* read source information */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your source information table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    fscanf(InF," %c", &SourceType);
    if ( SourceType!='I' && SourceType!='V' && SourceType!='P' )
    {
       fprintf(stderr, "Check your source type\n");
       exit(1);
    }


    /*Reading Voltage Source information*/
    if( SourceType=='V' )
    {
        fscanf(InF,"%lf",&OperateFreq);        
        fscanf(InF,"%d",&VSourceNum);

        for(IEdge=0;IEdge<VSourceNum;IEdge++)
        {
            fscanf(InF,"%d %lf",&VSourceEdge[IEdge],&VSourceMag[IEdge]);
        }   /* for(IEdge */
    }      

    /* Reading plane wave information */
    if( SourceType=='P' )
    {
        fscanf(InF,"%lf%lf%lf%lf%lf%lf", &OperateFreq,&E_Theta, 
                    &E_Phi, &K_Theta, &K_Phi,&PlaneWaveE_Mag);
    } 

    /* read the current source information */
    if( SourceType=='I' )
    {  
        fscanf(InF, "%lf %d", &OperateFreq, &TotISourceEdgeNum);
 
        ISourceEdgeStat=INT_Vector(TotISourceEdgeNum);
        ISourceEdgeMag=double_Vector(TotISourceEdgeNum);

        for(i=0; i<TotISourceEdgeNum; i++)
            fscanf(InF,"%d %lf", &ISourceEdgeStat[i], &ISourceEdgeMag[i]); 
    }    

    /* read junction table information */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your junctionn table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);
    fscanf(InF,"%d\n", &JunctionNum);
    for(i=0; i<JunctionNum; i++)
        fscanf(InF,"%d  %d  %d\n", &Junction[i][0],&Junction[i][1], 
               &Junction[i][2]);
 

    /* read the resistor table */ 
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your resistor table\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);  /* skip the comment line */
    fscanf(InF, "%s %d", buff, &ResistorEdgeNum);
    if( !strcmp(buff, "resistor") )
    {
        int i;
        for(i=0; i<ResistorEdgeNum; i++)
        {
            fscanf(InF, "%d %lf", ResistorEdgeStat+i, ResistorEdgeValue+i);
        }  /* for(i=0 */
    }      /* if( !strcmp */

    else 
    {
        fprintf(stderr, "Check the resistor table\n");
        exit(1);
    }   /* else */
        

}    /* end of ReadInputFile() */



/***************************************************************************
Prototype:     void   PrintOutput() 
Description:   To print out the field within the area specified by the 
               input file . 
Input value:     none 
Return value:   none 
Global value used:  InF, OutF, HybrdBoundEdgeNum, EdVector, JdVector, 
                    JcVector, TetGlobalEdgeEnds, HybrdBoundEdgeNum, 
                    TotInnerEdgeNum, NodeCord, RHSVector 
Global value modified:     none 
Subroutines called:  none 
****************************************************************************/
void PrintOutput()
{
    char buff[100];
    char c;
    int flag, II, JJ, IEdge, JEdge;

    /* skip the comment line before the output table  */
    fscanf(InF, " %c", &c);
    if (c!='#') 
    {
       fprintf(stderr, "Check your output information\n");
       exit(1);
    }
    fgets(buff, sizeof(buff)-1, InF);

    fscanf(InF, "%s", buff);  /* read "default_out", just skip it  */
    fscanf(InF, "%d", &flag);
 
    /* need to print out default output */
    if (flag==1) 
    { 
 
        fscanf(InF, " %s", buff);  /* get the real name */
        if( (OutF=fopen(buff, "w"))==NULL ) 
        {
            fprintf(stderr, "Can't open output file\n");
            exit(1);
        }

        /* equivalent surface magnetic currents */
        for(II=0; II<HybrdBoundEdgeNum;II++)    
            fprintf(OutF,"%e %e\n", EdVector[II].x, EdVector[II].y);

            /* E field vanish in the metal patch, don't need to print out */

            fprintf(OutF,"\n\n");

            /* equivalnet surface electric currents */
            for(II=0; II<HybrdBoundEdgeNum;II++)   
                fprintf(OutF,"%e %e\n",  JdVector[II].x, JdVector[II].y) ;     
     
            for(II=0; II<TotExtMetalEdgeNum; II++)
                fprintf(OutF,"%e  %e\n", JcVector[II].x, JcVector[II].y);

            fprintf(OutF,"\n\n");
            fclose(OutF);

    } /* if(flag */

    fscanf(InF, "%d", &JJ);
    for(II=0; II<JJ; II++) 
    {
        double start[3], end[3];
        double Min[3], Max[3];
        char  axis;
        
        char  name[20];
        int KK;
        int AxisFlag[3],Out_flip=0;
          
        fscanf(InF, " %s", name);
        
        if( (OutF=fopen(name, "a"))==NULL ) 
        {
            fprintf(OutF1, "Can't open file %s\n", name);
            exit(1);
        }
        
        /* get coordiantes that we are interested in */
        fscanf(InF, "%lf %lf %lf %lf %lf %lf %c", Min, Min+1, Min+2, 
               Max, Max+1, Max+2, &axis );
        fprintf(OutF, "\n\n J within (%lf, %lf, %lf) and (%lf, %lf, %lf)\n\n", Min[0], Min[1], Min[2], Max[0], Max[1], Max[2]);

        
        AxisFlag[0]=AxisFlag[1]=AxisFlag[2]=1;
        if( fabs(Min[0]-Max[0])<1e-6 ) AxisFlag[0]=0;  /* disable x axis */
        else if(fabs(Min[1]-Max[1])<1e-6) AxisFlag[1]=0;
        else if(fabs(Min[2]-Max[2])<1e-6) AxisFlag[2]=0;
        else 
        {
            fprintf(stderr, "Please specify a rectagnle area\n");
            exit(1);
        }

        /* make it bigger thus it is easy to test whether edges 
           are within the rectangle */
        for(IEdge=0; IEdge<3; IEdge++ ) 
        {
            Min[IEdge] -= 1e-6;
            Max[IEdge] += 1e-6;
        }
        for(IEdge=0; IEdge<TotBoundEdgeNum; IEdge++)
        {
            complex *Result;
            int OutAxis, Index, flag;
        

            for(JEdge=0; JEdge<3; JEdge++)
            { 
                int EdgeNum, StartNodeNum, EndNodeNum;

                EdgeNum = BoundEdgeStat[IEdge];
                StartNodeNum = GlobalEdgeEnds[ EdgeNum-1][0];
                EndNodeNum = GlobalEdgeEnds[ EdgeNum-1][1];

                /* get coordiantes of IEdge */
                start[JEdge]=NodeCord[StartNodeNum-1][JEdge];
                end[JEdge] = NodeCord[EndNodeNum-1][JEdge]; 
            }
            
            flag=0; /* test whether it is out of bound */
            
            for(KK=0; KK<3; KK++)
                if( start[KK]<Min[KK] || end[KK]>Max[KK]) 
                { 
                   /* out of bound */
                   flag=1;
                   break;
                }
                  
            /* don't need to check this edge, thus jump to the next edge */
            if( flag ) continue;


            /* those edges, which are not parallel to the axis we specified, 
               will be set to 0, thus won't be print out */

            flag=0;
            switch(axis)
            {
                case 'x': if( fabs(start[1]-end[1])<1e-6 &&  
                              fabs(start[2]-end[2])<1e-6 ) 
                          {
                              if( AxisFlag[1]==0 ) 
                                  OutAxis=2;  /* print out Z cooridnate */
                              else   OutAxis=1; /* print out Y coordinate */                               
                              flag=1;                       
                          }                     
                          break;

                case 'y': if( fabs(start[0]-end[0])<1e-6 && 
                              fabs(start[2]-end[2])<1e-6 ) 
                          {
                              if( AxisFlag[2]==0 )
                                  OutAxis=0;    /* print out Y coord */
                              else OutAxis=2;   /* print out Z coord */

                              flag=1;
	                  }      

                     
                          break;

                case 'z': if( fabs(start[0]-end[0])<1e-6 && 
                              fabs(start[1]-end[1])<1e-6)  
                          {
                              if( AxisFlag[0]==0 ) 
                                  OutAxis=1;   /* print out Y coordiante */
                              else OutAxis=0;  /* print out X coordinate */

                              flag=1;
                          }
	                  break;

                default: fprintf(stderr, "the output definition has problems\n");
                         exit(1);
            }    

            if (!Out_flip && flag==1) 
            {
               char axis[3]={'X', 'Y', 'Z'};
      
               Out_flip = !Out_flip;
              
               fprintf(OutF, "\n%c coordinate\tJ: Real\t\tJ:Imag\t(Unit A/m)\n\n", axis[OutAxis]);
            }

            if (flag==1 ) 
            {
                if( IEdge<HybrdBoundEdgeNum ) 
                {
                    /* it is a dielectric edge */
                    Result=JdVector;    
		    Index =IEdge;
                  
                }  /* if(IEdge */

                else  
                {
                    /* it is a mettalic edge */

                    Result=JcVector;
		    Index=IEdge-HybrdBoundEdgeNum;

                }  /* else */

                fprintf(OutF, "%lf\t%lf\t%lf\n", start[OutAxis],  
                              Result[Index].x, Result[Index].y);
            }  /* if(flag  */              
                       
        } /* for(IEdge=0; */


        /* plot out the field within the FEM region */
        for(IEdge=0; IEdge<TotInnerEdgeNum; IEdge++)
        {
            for(JEdge=0; JEdge<3; JEdge++)
            { 

                int EdgeNum, StartNodeNum, EndNodeNum;
                     
                EdgeNum = InnerEdgeStat[IEdge];
                StartNodeNum = GlobalEdgeEnds[ EdgeNum -1 ][0];
                EndNodeNum = GlobalEdgeEnds[ EdgeNum -1 ][1];

                /* get coordiantes of IEdge */
                start[JEdge]=NodeCord[ StartNodeNum -1 ] [JEdge];
                end[JEdge] = NodeCord[ EndNodeNum -1 ] [JEdge]; 

            }  /* for(JEdge */

            flag=0;  /* test whether it is out of bound */
            
            for(KK=0; KK<3; KK++)
                if( start[KK]<Min[KK] || end[KK]>Max[KK] ) 
                {
                    /* out of bound */
                    flag=1;
                    break;
                }
  
            if (flag) continue;

            /* those edges, which are not parallel to the axis we specified,
               will be set to 0, thus won't be print out */

            if (!Out_flip)
            {
 		fprintf(OutF, "\n%c coordinate(mid)\tE: Real\t\t E:Imag\t(Unit:V/m)\n\n", axis);
                Out_flip=!Out_flip;
            }
      

            switch(axis)
            {
                
                case 'x': if( fabs(start[1]-end[1])<1e-6 && 
                              fabs(start[2]-end[2])<1e-6 ) 
                              fprintf(OutF, "%lf\t%lf\t%lf\n",(start[0]+end[0])/2.0,  
                                               RHSVector[IEdge].x, RHSVector[IEdge].y );                         
                          break;

                case 'y': if( fabs(start[0]-end[0])<1e-6 && 
                              fabs(start[2]-end[2])<1e-6 )                
                          fprintf(OutF, "%lf\t%lf\t%lf \n", (start[1]+end[1])/2.0,
                                   RHSVector[IEdge].x, RHSVector[IEdge].y );
	                                          
                          break;

                case 'z': if( fabs(start[0]-end[0])<1e-6 && 
                              fabs(start[1]-end[1])<1e-6 )             
                              fprintf(OutF, "%lf\t%lf\t%lf\n", 
                                      (start[2]+end[2])/2.0, 
                                      RHSVector[IEdge].x, RHSVector[IEdge].y );
	                  break;
                default:  fprintf(stderr,"Check output specifications\n");
                          exit(1);
            }  /* switch */
        }  /* for(IEdge=0 */  

        fprintf(OutF,"\n\n");  
        fclose(OutF);

        Out_flip=0;    /* reset the output flip */
        
    }  /* for(II=0; */

}      /* end of PrintOutput() */



  
/*************************** END OF CODE ***********************************/

