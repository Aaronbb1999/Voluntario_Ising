#include <stdio.h>
#include "math.h"
#include "gsl_rng.h"

gsl_rng *tau;

//Parametros que podemos cambiar
#define N 16
#define T 2.2861
#define pmc 1000000

int main(void)
{
    double p,et;
    int i,j,k,l,aux;
    int nd,na,md,ma,E,n,m;
    int s[N][N];
    extern gsl_rng *tau;

    double magn=0, magne=0, energ=0, energe=0, energsq=0, cn=0, aux1=0, aux2=0, aux3=0;
    double f[N];
    

    for(l=0; l<N; l++)
    {
        f[l]=0;
    }


    FILE *fout;
    fout=fopen("Resultados.dat","w");

    int semilla=198756;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);


    //Tomamos un estado inicial ordenado
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            s[i][j]=1;
        }
    }


    //Algoritmo de Metropolis
    for(i=0;i<pmc;i++)
    {
    
        for(j=0;j<(N*N);j++)
        {
            //Empezamos el algoritmo tomando una posición aleatoria (n,m)
            n=gsl_rng_uniform_int(tau,N);
            m=gsl_rng_uniform_int(tau,N);
    
    
            //Aplicamos las condiciones periódicas
            if(n==0)
            {
                nd=1;
                na=N-1;
            }
            else if(n==N-1)
            {
                nd=0;
                na=N-2;
            }
            else
            {
                nd=n+1;
                na=n-1;
            }
            if(m==0)
            {
                md=1;
                ma=N-1;
            }
            else if(m==N-1)
            {
                md=0;
                ma=N-2;
            }
            else
            {
                md=m+1;
                ma=m-1;
            }
    
            //Calculamos E y comprobamos si hay que cambiar o no la orientación del spin.
            E=2*s[n][m]*(s[nd][m]+s[na][m]+s[n][md]+s[n][ma]);
    
            p=1;
            if(exp(-E/T)<p)
            {
                p=exp(-E/T);
            }
    
            et=gsl_rng_uniform(tau);
            if(et<p)
            {
                s[n][m]=-s[n][m];
            }
        }

        if(i%100==0)
        {
            //Energía y magnetización
            aux1=0;
            aux2=0;
            for(j=0; j<N; j++)
            {
                for(k=0; k<N; k++)
                {
                    aux1=aux1+s[j][k];//Auxiliar para hallar magnetizacion

                    if(j==0)
                    {
                        nd=1;
                        na=N-1;
                    }
                    else if(j==N-1)
                    {
                        nd=0;
                        na=N-2;
                    }
                    else
                    {
                        nd=j+1;
                        na=j-1;
                    }
                    if(k==0)
                    {
                        md=1;
                        ma=N-1;
                    }
                    else if(k==N-1)
                    {
                        md=0;
                        ma=N-2;
                    }
                    else
                    {
                        md=k+1;
                        ma=k-1;
                    }
                    aux2=aux2-0.5*s[j][k]*(s[nd][k]+s[na][k]+s[j][md]+s[j][ma]);//Auxiliar para hallar energia
                }
            }
            magn=magn+abs(aux1)/(N*N);//Vamos sumando a la magnetizacion y energia los resultados auxiliares para luego promediar
            magne=magne+abs(aux1)/(N*N)*abs(aux1)/(N*N);
            energ=energ+aux2;
            energsq=energsq+(aux2*aux2);


            //Funcion de correlación
            for(l=0;l<N;l++)
            {
                aux3=0;
                for(j=0;j<N;j++)
                {
                    for(k=0;k<N;k++)
                    {
                        if(j+l>=N)
                        {
                            aux3=aux3+s[j][k]*s[j+l-N][k];
                        }
                        else{
                            aux3=aux3+s[j][k]*s[j+l][k];
                        }
                    }
                }
                f[l]=f[l]+aux3;
            }
        }
    }


    //Calculo de las magnitudes
    double medidas=(pmc/100);
    magn=magn/medidas;
    magne=(sqrt(magne/medidas-magn*magn)/sqrt(medidas))*2;
    aux1=(energ/medidas);
    energ=aux1/(2*N*N);
    energe=2*(sqrt((energsq/medidas-aux1*aux1)/(2*N*N)))/sqrt(medidas);
    cn=abs(energsq/medidas-aux1*aux1)/(N*N*T);

    for(l=0;l<N;l++)
    {
        f[l]=f[l]/medidas;
    }
    
    //Pasamos a ficheros las magnitudes
    fprintf(fout, "Magnetizacion = %lf    +- %lf\n", magn, magne);
    fprintf(fout, "Energia Media = %lf    +- %lf\n", energ, energe);
    fprintf(fout, "Calor especifico = %lf\n\n",cn);

    
    fprintf(fout,"Funcion de correlacion:\n");
    for(l=0; l<N; l++)
    {
        fprintf(fout, "f(%i) = %lf\n", l, f[l]/(N*N));
    }
    
    fclose(fout);

    return 0;
}