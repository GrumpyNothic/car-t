#include <stdio.h>
#include <math.h>


int main(void) {
    /* Código para simulação do cenário de reprodução
    do Paciente B (B11-CLL) do poster do SIMMA (2024)
    */

    double beta = 2.4577, eta = 2.0e-8;
    double rmin = 1.0, p1 = 4.7136, p2, p3 = 8.5e1, A = 9.5, xi = 1.2584, theta = 6.0e-6, alpha = 5.5e-7;
    double epsilon = 4.1e-4, mu = 2.0e-2;
    double lambda = 1.0e-1;
    double r = 0.176, b=0.5e-12, gamma = 2.25, d= 3.05e-1, a = 1.0e3;
 
    double dt = 1.0e-6;
    double tk = 0;

    double C_D1, C_D2, C_D3, C_D4, C_T1, C_T2, C_T3, C_T4, C_M1, C_M2, C_M3, C_M4;
    double kt, Ce, prol, cit, dist, citt, func, func1, func2, cmest, cest, incl1, incl2, incl3, incl4, incl5, incl6;
    double k10, k20, k30, k40, k11, k12, k13[2001], k21, k22, k23[2001], k31, k32, k33[2001], k41, k42, k43[2001];
    double C_D_inicial = 1.86e8, C_Dose2 = 0, C_Dose3 = 0, C_D, C_T = 0.0, C_M = 0.0;
    
    int n = 1999, i, j, limiar;
    double xbar = 0.42, hx, aux[2001], T[2001], T1[2001], T2[2001], T3[2001], T4[2001]; // k_4 alterados para [2001] ali em cima 
    double carga = 1.0e7, xn = 0.2992, xp = 0.6113, sd = 1e-2, xises[2001], gx[2001];
    double rbar, rsc, rsc1, rsc2, rsc3, rsc4, sc=0, rc=0, sc1, sc2, sc3, sc4, rc1, rc2, rc3, rc4;
    int cont = 0, intervalo = 100;
    double g0 = 0, p = 0, q = 1e4;
    double var1 = 0, var2 = 0;

    FILE *arq, *d1, *d2, *d3, *d4, *d5, *d6, *d7, *d8;

    arq = fopen("pb2_dados.txt", "w");
    d1 = fopen("pb2_dados_tumor_t0.txt", "w");
    d2 = fopen("pb2_dados_tumor_t30.txt", "w");
    d3 = fopen("pb2_dados_tumor_t60.txt", "w");
    d4 = fopen("pb2_dados_tumor_t90.txt", "w");
    d5 = fopen("pb2_dados_tumor_t180.txt", "w");
    d6 = fopen("pb2_dados_heat_taxa.txt", "w");
    d7 = fopen("pb2_dados_heat_efeito.txt", "w");
    d8 = fopen("pb2_dados_barplot.txt", "w");
    
    p2=1.0/(pow(7, p3));
    C_D = C_D_inicial;
    kt = rmin + p1/(1+p2*pow(tk,p3));
    
    hx = (double)1.0/(n+1);
    for(i=0; i<n+2; i++){
        xises[i] = hx*i;
        // Carga Tumoral Inicial
        if (xises[i] <= xbar){ // resistentes
            aux[i] = 1e3/(sd*sqrt(2*M_PI)) * exp( -0.5 * pow((xises[i]-xn)/sd, 2));
            T[i] = 0;
        }
        else{ // sensiveis
            aux[i] = 0;
            T[i] = carga/(sd*sqrt(2*M_PI)) * exp( -0.5 * pow((xises[i]-xp)/sd, 2));
        }
        T[i] = T[i] + aux[i];
        if (T[i] < 1e-16){
            T[i] = 0;
        }
        gx[i] = (pow(xises[i], 100)) / (pow(xbar, 100) + pow(xises[i], 100));
    }
    
    for (i=0; i<n+2; i++){ // salvando a distribuição inicial do tumor
        fprintf(d1, "%e\t%e\n", xises[i], T[i]);
    }
    fclose(d1);
    
    limiar = (xbar/hx);
    // Sensiveis
    sc = T[limiar];
    for (i=limiar+1; i<n+1; i++) {
        sc = sc + (2*T[i]);
    }
    sc = sc + T[n+1];
    sc = (sc * hx)/2;
    // Resistentes
    rc = T[0];
    for (i=1; i<limiar; i++) {
        rc = rc + (2*T[i]);
    }
    rc = rc + T[limiar];
    rc = (rc * hx)/2;
    // Totais
    rsc = rc + sc;
    
    rbar = r*(1-b*rsc)*(1-((p+q)/(rsc+q)));
    
    var1 = alpha*rsc*C_T;
    var2 = (C_D+C_T)/(a+C_D+C_T+d*sc);
    fprintf(arq, "%e %e %e %e %e %e %e %e %e %e %e %e\n", tk, C_D, C_T, C_M, C_D+C_T+C_M, sc, rc, rsc, var1, var2, kt, rbar);

    while(tk < 180.0 + dt){
        kt = rmin + p1/(1+p2*pow(tk,p3));

        C_D1 = C_D;
        C_T1 = C_T; 
        C_M1 = C_M;
        for (i=0; i<n+2; i++){
            T1[i] = T[i];
        } 
        
        // Sensiveis
	sc1 = T1[limiar];
	for (i=limiar+1; i<n+1; i++) {
	    sc1 = sc1 + (2*T1[i]);
	}
	sc1 = sc1 + T1[n+1];
	sc1 = (sc1 * hx)/2;
	// Resistentes
	rc1 = T1[0];
	for (i=1; i<limiar; i++) {
	    rc1 = rc1 + (2*T1[i]);
	}
	rc1 = rc1 + T1[limiar];
	rc1 = (rc1 * hx)/2;
	// Totais
	rsc1 = rc1 + sc1;

        k10 = -C_D1*(beta+eta);
    	k11 = eta*C_D1 + kt*((sc1*C_T1)/(A+sc1)) - (xi+epsilon+lambda)*C_T1 + theta*sc1*C_M1 - alpha*rsc1*C_T1; //efec
    	k12 = epsilon*C_T1 - theta*sc1*C_M1 - mu*C_M1; //mem
        for (i=0; i<n+2; i++){
            k13[i] = r*T1[i]*(1-b*rsc1)*(1-((p+q)/(rsc1+q))) - ( (1-gx[i])*g0 + gamma*gx[i] ) *( (C_D1+C_T1) / (a+C_D1+C_T1+d*sc1) ) *T1[i]; //tum
            if (fabs(k13[i]) < 1e-16){ k13[i] = 0; }
        }
        
        kt = rmin + p1/(1+p2*pow((tk+0.5*dt),p3));
    			
    	C_D2 = C_D + 0.5*dt*k10;
    	if (C_D2 < 1e-16){ C_D2 = 0; }
    	C_T2 = C_T + 0.5*dt*k11;
    	if (C_T2 < 1e-16){ C_T2 = 0; }
    	C_M2 = C_M + 0.5*dt*k12;
    	if (C_M2 < 1e-16){ C_M2 = 0; }
    	for (i=0; i<n+2; i++){
            T2[i] = T[i]+ 0.5*dt*k13[i];
            if (T2[i] < 1e-16){ T2[i] = 0; }
        } 
        
        // Sensiveis
	sc2 = T2[limiar];
	for (i=limiar+1; i<n+1; i++) {
	    sc2 = sc2 + (2*T2[i]);
	}
	sc2 = sc2 + T2[n+1];
	sc2 = (sc2 * hx)/2;
	// Resistentes
	rc2 = T2[0];
	for (i=1; i<limiar; i++) {
	    rc2 = rc2 + (2*T2[i]);
	}
	rc2 = rc2 + T2[limiar];
	rc2 = (rc2 * hx)/2;
	// Totais
	rsc2 = rc2 + sc2;

        k20 = -C_D2*(beta+eta);
    	k21 = eta*C_D2 + kt*((sc2*C_T2)/(A+sc2)) - (xi+epsilon+lambda)*C_T2 + theta*sc2*C_M2 - alpha*rsc2*C_T2; //efec
    	k22 = epsilon*C_T2 - theta*sc2*C_M2 - mu*C_M2; //mem
        for (i=0; i<n+2; i++){
            k23[i] = r*T2[i]*(1-b*rsc2)*(1-((p+q)/(rsc2+q))) - ( (1-gx[i])*g0 + gamma*gx[i] ) *( (C_D2+C_T2) / (a+C_D2+C_T2+d*sc2) ) *T2[i]; //tum
            if (fabs(k23[i]) < 1e-16){ k23[i] = 0; }
        }
        
        C_D3 = C_D + 0.5*dt*k20;
        if (C_D3 < 1e-16){ C_D3 = 0; }
        C_T3 = C_T + 0.5*dt*k21; 
        if (C_T3 < 1e-16){ C_T3 = 0; }
    	C_M3 = C_M + 0.5*dt*k22;
    	if (C_M3 < 1e-16){ C_M3 = 0; }
    	for (i=0; i<n+2; i++){
            T3[i] = T[i]+ 0.5*dt*k23[i];
            if (T3[i] < 1e-16){ T3[i] = 0; }
        }
        
        // Sensiveis
	sc3 = T3[limiar];
	for (i=limiar+1; i<n+1; i++) {
	    sc3 = sc3 + (2*T3[i]);
	}
	sc3 = sc3 + T3[n+1];
	sc3 = (sc3 * hx)/2;
	// Resistentes
	rc3 = T3[0];
	for (i=1; i<limiar; i++) {
	    rc3 = rc3 + (2*T3[i]);
	}
	rc3 = rc3 + T3[limiar];
	rc3 = (rc3 * hx)/2;
	// Totais
	rsc3 = rc3 + sc3;

        k30 = -C_D3*(beta+eta);
    	k31 = eta*C_D3 + kt*((sc3*C_T3)/(A+sc3)) - (xi+epsilon+lambda)*C_T3 + theta*sc3*C_M3 - alpha*rsc3*C_T3; //efec
    	k32 = epsilon*C_T3 - theta*sc3*C_M3 - mu*C_M3; //mem
        for (i=0; i<n+2; i++){
            k33[i] = r*T3[i]*(1-b*rsc3)*(1-((p+q)/(rsc3+q))) - ( (1-gx[i])*g0 + gamma*gx[i] ) *( (C_D3+C_T3) / (a+C_D3+C_T3+d*sc3) ) *T3[i]; //tum
            if (fabs(k33[i]) < 1e-16){ k33[i] = 0; }
        }
        
        kt = rmin + p1/(1+p2*pow((tk+dt),p3));

        C_D4 = C_D + dt*k30;
        if (C_D4 < 1e-16){ C_D4 = 0; }
    	C_T4 = C_T + dt*k31;
    	if (C_T4 < 1e-16){ C_T4 = 0; }
    	C_M4 = C_M + dt*k32;
    	if (C_M4 < 1e-16){ C_M4 = 0; }
    	for (i=0; i<n+2; i++){
            T4[i] = T[i]+ dt*k33[i];
            if (T4[i] < 1e-16){ T4[i] = 0; }
        }
        
        // Sensiveis
	sc4 = T4[limiar];
	for (i=limiar+1; i<n+1; i++) {
	    sc4 = sc4 + (2*T4[i]);
	}
	sc4 = sc4 + T4[n+1];
	sc4 = (sc4 * hx)/2;
	// Resistentes
	rc4 = T4[0];
	for (i=1; i<limiar; i++) {
	    rc4 = rc4 + (2*T4[i]);
	}
	rc4 = rc4 + T4[limiar];
	rc4 = (rc4 * hx)/2;
	// Totais
	rsc4 = rc4 + sc4;

        k40 = -C_D4*(beta+eta);
    	k41 = eta*C_D4 + kt*((sc4*C_T4)/(A+sc4)) - (xi+epsilon+lambda)*C_T4 + theta*sc4*C_M4 - alpha*rsc4*C_T4; //efec
    	k42 = epsilon*C_T4 - theta*sc4*C_M4 - mu*C_M4; //mem
        for (i=0; i<n+2; i++){
            k43[i] = r*T4[i]*(1-b*rsc4)*(1-((p+q)/(rsc4+q))) - ( (1-gx[i])*g0 + gamma*gx[i] ) *( (C_D4+C_T4) / (a+C_D4+C_T4+d*sc4) ) *T4[i]; // sensiveis
            if (fabs(k43[i]) < 1e-16){ k43[i] = 0; }
        }
        
        // Solução do passo de tempo
        C_D = C_D + (dt/6.0)*(k10 + 2.0*k20 + 2.0*k30 + k40);
        if (C_D < 1e-16){ C_D = 0; }
        C_T = C_T + (dt/6.0)*(k11 + 2.0*k21 + 2.0*k31 + k41);
        if (C_T < 1e-16){ C_T = 0; }
    	C_M = C_M + (dt/6.0)*(k12 + 2.0*k22 + 2.0*k32 + k42);
    	if (C_M < 1e-16){ C_M = 0; }
	for (i=0; i<n+2; i++){
	    T[i] = T[i] + (dt/6.0) * (k13[i] + 2.0*k23[i] + 2.0*k33[i] + k43[i]);
	    if (T[i] < 1e-16){ T[i] = 0; }
	}
	
	
	// Sensiveis
	sc = T[limiar];
	for (i=limiar+1; i<n+1; i++) {
	    sc = sc + (2*T[i]);
	}
	sc = sc + T[n+1];
	sc = (sc * hx)/2;
	// Resistentes
	rc = T[0];
	for (i=1; i<limiar; i++) {
	    rc = rc + (2*T[i]);
	}
	rc = rc + T[limiar];
	rc = (rc * hx)/2;
	// Totais
	rsc = rc + sc;
        
        tk = tk + dt;
        kt = rmin + p1/(1+p2*pow(tk,p3));
        if (fabs(tk - 1.0) < 1e-7){ // infusão de doses posteriores
            C_D = C_D + C_Dose2;
        }
        if (fabs(tk - 2.0) < 1e-7){
            C_D = C_D + C_Dose3;
        }
        // Salvamento dos dados dia 30
        if (fabs(tk - 30.0) <= 1e-7){ 
            fprintf(d8, "%.0f\t%e\t%e\t%e\n", tk, sc, rc, rsc); // Barplot t30
            for (i=0; i<n+2; i++){ // salvando a distribuição do tumor
                fprintf(d2, "%e\t%e\n", xises[i], T[i]);
            }
        }
        // Salvamento dos dados dia 60
        if (fabs(tk - 60.0) <= 1e-7){
            fprintf(d8, "%.0f\t%e\t%e\t%e\n", tk, sc, rc, rsc); // Barplot t60
            for (i=0; i<n+2; i++){ // salvando a distribuição do tumor
                fprintf(d3, "%e\t%e\n", xises[i], T[i]);
            }
        }
        // Salvamento dos dados dia 60
        if (fabs(tk - 90.0) < 1e-6){
            fprintf(d8, "%.0f\t%e\t%e\t%e\n", tk, sc, rc, rsc); // Barplot t90
            for (i=0; i<n+2; i++){ // salvando a distribuição do tumor
                fprintf(d4, "%e\t%e\n", xises[i], T[i]);
            }
        }
        
        if (rsc > 1e-16){
            rbar = r*(1-b*rsc)*(1-((p+q)/(rsc+q)));
        }
        else{
            rbar = 0;
        }

        cont++;
        
        // Salvando ao longo do tempo
        if (cont % 100 == 0) {
            var1 = alpha*rsc;
            var2 = (C_D+C_T) / (a+C_D+C_T+d*sc);
            fprintf(arq, "%e %e %e %e %e %e %e %e %e %e %e %e\n", tk, C_D, C_T, C_M, C_D+C_T+C_M, sc, rc, rsc, var1, var2, kt, rbar);
        }
        
        // Salvando para o heatmap
        if (cont % 10000 == 0){
            for (i=0; i<n+2; i++){
                aux[i] = ( (1-gx[i])*g0 + gamma*gx[i] ) * var2;
                fprintf(d6, "%e %e %e\n", tk, xises[i], aux[i]); // heatmap da taxa
                aux[i] = aux[i] * T[i];
                fprintf(d7, "%e %e %e\n", tk, xises[i], aux[i]); // heatmap do efeito
            }
            fprintf(d6, "\n");
            fprintf(d7, "\n");
        }
        
        // Print no terminal
        printf ("%e Ct: %e sc: %e rc: %e Tt: %e rbar: %f\n", tk, C_T, sc, rc, rsc, rbar);
    }
    
    // Salvando a distribuição final do tumor
    for (i=0; i<n+2; i++){ 
        fprintf(d5, "%e\t%e\n", xises[i], T[i]);
    }
    fclose(d2);
    fclose(d3);
    fclose(d4);
    fclose(d5);
    
    // Barplot último instante
    fprintf(d8, "%.0f\t%e\t%e\t%e\n", tk, sc, rc, rsc);
    fclose(d8);
    
    // Fechando arquivo dados ao longo do tempo
    fclose(arq);
    // Fechando arquivo heatmap
    fclose(d7);
    return 0;
}
