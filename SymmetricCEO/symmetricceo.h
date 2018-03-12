#include <math.h>
#include <limits.h>

#define DEFAULT_M 10
#define DEFAULT_STEPS 32
#define DEFAULT_LRHO 0.05
#define DEFAULT_HRHO 0.5
#define DEFAULT_OUTPUT_DIR "outdata"

int get_ber_symmetric_binary_ceo_majority_criterion(double rho,int M,double* ber);

void help(int argc,char *argv[])
{
    
    printf("%s\n\n",argv[0]);

    printf("\tGrafica teorica y experimental de la probabilidad de error BER\n");
    printf("\tde estimar la entrada u0 de un bloque de M canales BSC, conociendo\n");
    printf("\tunicamente la salida. Es usado en la parte experimental el criterio \n");
    printf("\tde la mayoria visto en: http://dx.doi.org/10.1109/SPAWC.2008.4641566\n\n");

    printf("\t%s --sources 10 --steps-rho 32 --low-rho 0.05 --high-rho  0.5\n\n",argv[0]);

    printf("\t\t--sources       M\n");
    printf("\t\t                Número de fuentes en el test.\n");
    printf("\t\t                Por defecto: %d\n\n",DEFAULT_M);

    printf("\t\t--steps-rho     STEPS\n");
    printf("\t\t                Número valores de rho analizados. Estos serán\n");
    printf("\t\t                linearmente espaciados entre LRHO e HRHO.\n");
    printf("\t\t                Por defecto: %d\n\n",DEFAULT_STEPS);

    printf("\t\t--low-rho       LRHO\n");
    printf("\t\t                El menor valor de los rho analizados.\n");
    printf("\t\t                Por defecto: %f\n\n",DEFAULT_LRHO);

    printf("\t\t--high-rho      HRHO\n");
    printf("\t\t                El mayor valor de los rho analizados.\n");
    printf("\t\t                Por defecto: %f\n\n",DEFAULT_HRHO);

    printf("\t\t--output-dir    DIR\n");
    printf("\t\t                El directorio de salida.\n");
    printf("\t\t                Por defecto: %s\n\n",DEFAULT_OUTPUT_DIR);
}

////////////////////////////////////////////////////////////////////////////////
// BER Theorico
int get_theoretical_ber_symmetric_binary_ceo(const PdsVector* RHO,PdsRaNatural M,PdsVector* BER)
{
    return pds_vector_symetric_ber_bsc_model(RHO,M,BER);
}

////////////////////////////////////////////////////////////////////////////////
// BER Experimental
int get_experimental_ber_symmetric_binary_ceo(const PdsVector* RHO,int M,PdsVector* BER)
{
    int dat;
    PdsRaNatural i;
    double val=0;

    if(RHO==NULL)   return FALSE;
    if(BER==NULL)   return FALSE;
    if(M<1)         return FALSE;

    if(BER->Nel!=RHO->Nel)   return FALSE;

    for(i=0;i<RHO->Nel;i++)
    {
        dat=get_ber_symmetric_binary_ceo_majority_criterion(RHO->V[i], M,&val);
        if(dat<0)  return FALSE;
        BER->V[i]=val;
    }
    return TRUE;
}

double pds_log2(double x)
{
    return log(x)/log( 2 );
}

int get_ber_symmetric_binary_ceo_majority_criterion(double rho,int M,double* ber)
{
    unsigned long int err_count,counts;
    int d;
    
    PdsBscBlock *BLOCK=NULL;
    PdsCoin *U0=NULL;
    PdsBVector *U=NULL;

    PdsRvByte u0;
    PdsRvByte hat_u0;

    PdsBaNatural n;


    BLOCK=pds_bscblock_new_symmetric(M,rho);

    U0=pds_coin_new(0.5);

    U=pds_bvector_new(M);
    if( (BLOCK==NULL)||(U0==NULL)||(U==NULL) )
    {
        pds_bscblock_free(BLOCK);
        pds_coin_free(U0);
        pds_bvector_free(U);
        return -1;
    }

    err_count=0;
    counts=0;
    while(err_count<512)
    {
        pds_coin_get_bit (U0,&u0);
        
        pds_bscblock_evaluate_val(BLOCK,(PdsBaBit)u0,U);
        pds_bvector_weight_bvector (U, &n);
        if(n<(M/2))         hat_u0=0;
        else                hat_u0=1;
        /*else if(n>(M/2))    hat_u0=1;
        else
        {   
            d=rand();
            if(d>(RAND_MAX/2))
            hat_u0=1;
        }*/

        if(hat_u0!=u0)  
        {
            err_count=err_count+1;
            //printf("err_count:%4ld\n",err_count);
        }

        counts++;
        if( counts>=(ULONG_MAX-4) )
        {
            fprintf(stderr,"ERROR FATAL: La variable que cuenta intentos experimentales fue desbordada.\n");
            exit(0);
        }
        
    }
        
    *ber=(err_count*1.0)/counts;

    printf("Experimental\tcounts:%10ld -> rho:%f\tber:%e\n",counts,rho,*ber);
    return counts;
}


int pds_vector_set_hb(PdsVector* H,const PdsVector* P)
{
    int i;

    for(i=0;i<P->Nel;i++)   H->V[i]=pds_hb(P->V[i]);
    
}
