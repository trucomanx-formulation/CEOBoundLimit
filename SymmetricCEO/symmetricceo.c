/*
 * Asymmetricceo.c
 * 
 * Copyright 2011 Fernando Pujaico Rivera <fernando.pujaico.rivera@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>


#include <pds/pdsrv.h>
#include <pds/pdsra.h>
#include <pds/pdsit.h>
#include <pds/pdscm.h>
#include <pds/pdsba.h>
#include <pds/pdsdatafunc.h>
#include <pds/pdsoctplot.h>
#include <pds/pdsmath.h>
#include "symmetricceo.h"

int main(int argc, char** argv)
{

    PdsVector *T_BER=NULL;
    PdsVector *E_BER=NULL;
    PdsVector *RHO=NULL;

    PdsVector *T_HBER=NULL;
    PdsVector *T_HU0OMEGA=NULL;

    int M,STEPS,dat;
    double lrho,hrho;
    char *output_dir=NULL;

    M=DEFAULT_M;
    STEPS=DEFAULT_STEPS;
    lrho=DEFAULT_LRHO;
    hrho=DEFAULT_HRHO;

    if( (argc==1)||pds_exist_param(argc,argv,"--help")||pds_exist_param(argc,argv,"-h") )
    {
        help(argc,argv);
        return EXIT_SUCCESS;
    }

    pds_get_int_param (argc,argv,"--sources", &M);
    pds_get_int_param (argc,argv,"--steps-rho", &STEPS);
    pds_get_double_param (argc,argv,"--low-rho", &lrho);
    pds_get_double_param (argc,argv,"--high-rho", &hrho);
    pds_get_chars_param (argc,argv,"--output-dir", &output_dir);

    if(output_dir==NULL)
    {
        output_dir=(char*)calloc(16,sizeof(char));
        sprintf(output_dir,DEFAULT_OUTPUT_DIR);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Calculando BER


    //Creando vectores;
    RHO=pds_vector_new_linspace(lrho,hrho,STEPS);   

    T_BER=pds_vector_new(STEPS);
    E_BER=pds_vector_new(STEPS);

    // Calculando BER
    dat=get_theoretical_ber_symmetric_binary_ceo(RHO,M,T_BER);
    dat=get_experimental_ber_symmetric_binary_ceo(RHO,M,E_BER);

    // plotando datos
    mkdir(output_dir,S_IRWXU);
    pds_octave_semilogy_compare_vectors (RHO,T_BER,E_BER,
                                    "rho",pds_sprintf("BER, M=%d",M), 
                                    "Theoretical","Experimental",
                                    pds_sprintf("%s%cM%03d_ceo.m"  ,output_dir,pds_filesep(),M),
                                    pds_sprintf("%s%cM%03d_ceo.eps",output_dir,pds_filesep(),M));

    ////////////////////////////////////////////////////////////////////////////
    // Calculando entropia binaria.
    T_HBER=pds_vector_new(STEPS);
    T_HU0OMEGA=pds_vector_new(STEPS);


    pds_vector_set_hb(T_HBER,T_BER);
    pds_vector_symetric_entropy_u0_omega_bsc_model(RHO,M,T_HU0OMEGA);


    pds_octave_semilogy_compare_vectors (RHO,T_HBER,T_HU0OMEGA,
                                    "rho",pds_sprintf("H, M=%d",M), 
                                    "hb(BER)","H(U0|Omega)",
                                    pds_sprintf("%s%cM%03d_hbceo.m"  ,output_dir,pds_filesep(),M),
                                    pds_sprintf("%s%cM%03d_hbceo.eps",output_dir,pds_filesep(),M));
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////

