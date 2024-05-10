#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float p, V, T, R = 8.314;

/*Input functions*/

void getInputPressureVolumeTemperature()
{
    printf("\nEnter known Variable V(m3)=");
    scanf("%f", &V);
    printf(" T=");
    scanf("%f", &T);
}

void getInputPressureTemperature()
{
    printf("\nEnter known Variable p(Pa)=");
    scanf("%f", &p);
    printf(" T=");
    scanf("%f", &T);
}

void getInputTemperatureVolume()
{
    printf("\nEnter known Variable p=");
    scanf("%f", &p);
    printf(" V=");
    scanf("%f", &V);
}

void calculateLX(char sql)
{
    switch (sql)
    {
    case 'P':
    case 'p':
        getInputPressureVolumeTemperature();
        p = R * T / V;
        printf("\nCalculation formula p=R*T/V\nRequired Variable p=%f\n", p);
        break;
    case 'V':
    case 'v':
        getInputPressureTemperature();
        V = T * R / p;
        printf("\nCalculation formula V=T*R/p\nRequired Variable V=%f\n", V);
        break;
    case 'T':
    case 't':
        getInputTemperatureVolume();
        T = p * V / R;
        printf("\nCalculation formula T=p*V/R\nRequired Variable T=%f\n", T);
        break;
    default:
        printf("\n\nIncorrect input for the Variable requested\n");
    }
}

/*Virial equation*/

void calculateVL(char sql)
{
    float B, C, Z;
    float f[20], g[20];
    int i;
    f[0] = 1;
    printf("\nEnter Virial coefficient B=");
    scanf("%f", &B);
    printf("\nEnter Virial coefficient C (if none, please enter 0)=");
    scanf("%f", &C);
    if (C == 0)
    {
        switch (sql)
        {
        case 'P':
        case 'p':
            getInputPressureVolumeTemperature();
            Z = 1 + B / V;
            p = Z * R * T / V;
            printf("\nCalculation formula Z=1+B/V;\np=Z*R*T/V;\nCoefficient Z=%f\nRequired Variable p=%f\n", Z, p);
            break;
        case 'V':
        case 'v':
            getInputPressureTemperature();
            Z = 1 + B * p / R / T;
            V = Z * T * R / p;
            printf("\nCalculation formula Z=1+B*p/R/T;\nV=Z*T*R/p; \nCoefficient Z=%f\nRequired Variable V=%f\n", Z, V);
            break;
        default:
            printf("\n\nIncorrect input for the Variable requested\n ");
        }
    }
    if (C != 0)
    {
        switch (sql)
        {
        case 'P':
        case 'p':
            getInputPressureVolumeTemperature();
            Z = 1 + B / V + C / V / V;
            p = Z * R * T / V;
            printf("\nCalculation formula Z=1+B/V+C/V/V;\np=Z*R*T/V;\nCoefficient Z=%f\nRequired Variable p=%f\n", p);
            break;
        case 'V':
        case 'v':
            getInputPressureTemperature();
            for (i = 0; i < 20; i++)
            {
                g[i] = R * T / p * f[i];
                f[i + 1] = 1 + B / g[i] + C / g[i] / g[i];
                if (f[i + 1] - f[i] > -1e-6 && f[i + 1] - f[i] < 1e-6)
                    break;
            }
            Z = f[i + 1];
            V = Z * T * R / p;
            printf("\nAfter the %dth calculation,\nZ%d=%f;\nV%d=%f;", i, i, f[i + 1], i, g[i]);
            printf("\nCalculation formula\n V[i]=R*T/p*Z[i];\nV[i]=R*T/p*Z[i];Z[i+1]=1+B/g[i]+C/g[i]/g[i]; \nCoefficient Z=%f\nRequired Variable V=%f", Z, V);
            break;
        default:
            printf("\n\nIncorrect input for the Variable requested\n");
        }
    }
}

/*Van der Waals equation*/

void calculateVDW(char sql)
{
    float a, b;
    printf("\nEnter relevant parameter a=");
    scanf("%f", &a);
    printf("\nEnter relevant parameter b=");
    scanf("%f", &b);
    switch (sql)
    {
    case 'P':
    case 'p':
        getInputPressureVolumeTemperature();
        p = R * T / (V - b) - a / V / V;
        printf("\nCalculation formula\n p=R*T/(V-b)-a/V/V;\nRequired Variable p=%f\n", p);
        break;
    case 'V':
    case 'v':
        printf("\n\nNo relevant formula");
        break;
    default:
        printf("\n\nIncorrect input for the Variable requested");
    }
}

/*RK equation*/

void calculateRK(char sql)
{
    float Tc, Pc, Z, a, b, A, B;
    float f[20], g[20];
    int i;
    f[0] = 1;
    printf("\nEnter critical temperature Tc(K)=");
    scanf("%f", &Tc);
    printf("\nEnter critical pressure Pc(Pa)=");
    scanf("%f", &Pc);
    a = 0.42748 * R * R * Tc * Tc * sqrt(Tc) / Pc;
    b = 0.08664 * R * Tc / Pc;
    switch (sql)
    {
    case 'P':
    case 'p':
        getInputPressureVolumeTemperature();
        p = R * T / (V - b) - a / sqrt(T) / V / (V + b);
        printf("\nCalculation formula\n p=R*T/(V-b)-a/sqrt(T)/V/(V+b);\nRequired Variable p=%f\n", p);
        break;
    case 'V':
    case 'v':
        getInputPressureTemperature();
        A = a * p / R / R / T / T / sqrt(T);
        B = b * p / R / T;
        for (i = 0; i < 20; i++)
        {
            g[i] = B / f[i];
            f[i + 1] = 1 / (1 - g[i]) - A / B * g[i] / (1 + g[i]);
            printf("\nAfter the %dth calculation,\nZ%d=%f;\nh%d=%f;", i, i, f[i + 1], i, g[i]);
            if ((f[i + 1] - f[i]) / f[i + 1] > -1e-3 && (f[i + 

1] - f[i]) / f[i + 1] < 1e-3)
                break;
        }
        Z = f[i + 1];
        V = Z * T * R / p;
        printf("\nCalculation formula\n h[i]=B/Z[i];\nZ[i+1]=1/(1-h[i])-A/B*h[i]/(1+h[i]);\nCoefficient A=%f\nCoefficient B=%f\nCoefficient Z=%f\nRequired Variable V=%f\n", A, B, Z, V);
        break;
    default:
        printf("\n\nIncorrect input for the Variable requested\n");
    }
}

/*SRK(KRS) equation*/

void calculateSRK(char sql)
{
    float Tc, Pc, Tr, a, b, A, B, w, m, K;
    float f[20], g[20];
    int i;
    f[0] = 1;
    printf("\nEnter critical temperature Tc(K)=");
    scanf("%f", &Tc);
    printf("\nEnter critical pressure Pc(Pa)=");
    scanf("%f", &Pc);
    printf("\nEnter acentric factor w=");
    scanf("%f", &w);
    a = 0.42748 * R * R * Tc * Tc * sqrt(Tc) / Pc;
    b = 0.08664 * R * Tc / Pc;
    switch (sql)
    {
    case 'P':
    case 'p':
        getInputPressureVolumeTemperature();
        Tr = T / Tc;
        m = 0.480 + 1.574 * w - 0.176 * w * w;
        K = a * (1 + m * (1 - sqrt(Tc))) * (1 + m * (1 - sqrt(Tc)));
        p = R * T / (V - b) - K / V / (V + b);
        printf("\nCalculation formula\n p=R*T/(V-b)-a/ V/(V+b);\nCoefficient a=%f\nCoefficient b=%f\nCoefficient Tr=%f\nCoefficient K=%f\nRequired Variable p=%f\n", a, b, Tr, K, p);
        break;
    case 'V':
    case 'v':
        getInputPressureTemperature();
        A = a * p / R / R / T / T;
        B = b * p / R / T;
        printf("\nCalculation formula\n h[i]=B/Z[i];Z[i+1]=1/(1-h[i])-A/B*h[i]/(1+h[i]);\nRequired Variable V=%f\n", V);
        break;
    default:
        printf("\n\nIncorrect input for the Variable requested\n");
    }
}

/*Main function*/

int main()
{
    char Num, jg, sql;
    int i, n;
    printf("\n Liu Quan et al.， 数值软件提升化工热力学教学成效的量化研究. 化学教育（中英文），2024");
    printf("\n\nEquations of state:\n1 Ideal gas equation;\n2 Virial equation;\n3 Van der Waals equation;\n4 RK equation;\n5 SRK equation.");
    printf("\n\nEnter the required equation number Num=");
    scanf("%s", &Num);
    printf("\nSelect symbol for the Variable requested: p (pressure, Pa), V (volume, m3), T (temperature, K)\n\nEnter the requested Variable sql=");
    scanf("%s", &sql);
    switch (Num)
    {
    case '1':
        printf("\nYou chose to calculate using the Ideal Gas Equation\n");
        calculateLX(sql);
        break;
    case '2':
        printf("\nYou chose to calculate using the Virial Equation\n");
        calculateVL(sql);
        break;
    case '3':
        printf("\nYou chose to calculate using the Van der Waals Equation\n");
        calculateVDW(sql);
        break;
    case '4':
        printf("\nYou chose to calculate using the RK Equation\n");
        calculateRK(sql);
        break;
    case '5':
        printf("\nYou chose to calculate using the SRK Equation\n");
        calculateSRK(sql);
        break;
    default:
        printf("\nThe equation you selected does not exist\n");
    }
}
