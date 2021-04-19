///-----------------------------------------------------------------
///
/// @file      Source_Code.cpp
/// @author    Jeferson Diehl de Oliveira, Elaine Maria Cardoso and Jacqueline Biancon Copetti
/// Created:   5/18/2020 1:03:43 AM
///
///------------------------------------------------------------------

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main()
{
	
	printf("\n\n - - - - - - - - - - - - - - VFC - VOID FRACTION CALCULATOR - - - - - - - - - - - - - -\n\n");
	printf("\n");
	printf("VFC - Void Fraction Calculator is a registered software (Process Number: BR512020002495-5) \n\n");
	printf("\n");
	printf("Developed by (1)Jeferson Diehl de Oliveira, (2)Elaine Maria Cardoso and (3)Jacqueline Biancon Copetti\n\n");
	printf("\n");
	
	printf("(1)FSG - Centro Universitario, Campus Sede, Av. Rua Os Dezoito do Forte, 2366, 95020-472, Caxias do Sul, RS, Brasil.\n");
	printf("(2)UNESP - Sao Paulo State University, Post-Graduation Program in Mechanical Engineering, Av. Brasil, 56, 15385-000, Ilha Solteira, SP, Brazil.\n");
	printf("(3)University of Vale do Rio dos Sinos, Mechanical Engineering Graduate Program, Unisinos Avenue, 950, Sao Leopoldo, 93022-750, RS, Brazil.\n");
	printf("\nCreated: 05/18/2020\n");
	printf("Published: 05/20/2020\n");
	printf("\n");
	printf("How to cite: Oliveira, J.D., Cardoso, M.E., Copetti, J.M.,2020. VFC - VOID FRACTION CALCULATOR [Computer Software]. Retrieved from https://github.com/");
	printf("\n\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
	printf("\n");	
	
	
	 int i;
	 int N = 1000;
   
    double G; //Mass flow [kg/s]
    printf("\n\n Mass flux [kg/(m^2s)]: ");
    scanf("%f", &G);
    
    double D; //Hydraulic Diameter [m]
	printf("\n\n Diamter [m]: ");
    scanf("%f", &D);
	
	double g = 9.81 ; // Local gravity [m/s²]
	
	
	double P;  // Pressure or saturation pressure [kPa]
    printf("\n\n Pressure or Saturation pressure [kPa]: ");
    scanf("%f", &P);
	
	double rho_l; // Liquid-phase density  [kg/m³]
    printf("\n\n Liquid-phase density  [kg/m³]: ");
    scanf("%f", &rho_l);
	
	double rho_v; // Vapor-phase density  [kg/m³]
	printf("\n\n Vapor-phase density  [kg/m³]: ");
    scanf("%f", &rho_v);
	
	double mu_l;  // Liquid-phase viscosity  [Pa.s]
    printf("\n\n Liquid-phase Viscosity  [Pa.s]: ");
    scanf("%f", &mu_l);
    
	double mu_v;
	printf("\n\n Vapor-phase Viscosity  [Pa.s]: ");
   scanf("%f", &mu_v);
	
	double sigma; // Surface tension [N/m]
	printf("\n\n Surface tension  [N/m]: ");
    scanf("%f", &sigma);
	
	FILE * write;
        
     
	    
	double Re_lo = G*D/mu_l;
	double We_lo = pow(G,2)*D/(sigma*rho_l);
	double Fr_lo = pow(G,2)/(g*D*pow(rho_l,2));
	
	double x[N];
	
	double  b[N]; ///From Madsen (1975)
	double F_X_tt[N]; ///From Tandor et al. (1985)
    double Fr_tp[N]; ///From Guzhov et al. (1967)
    double Ft[N];     ///From Graham et al. (1997)
	double Re_l[N];
    double X_tt[N];
    double rho_tp[N];
    
    char quality[]= "Quality[]";
    char Alpha_Armand[] = "alpha_Armand[]";
    char Alpha_Bankoff[] = "Alpha_Bankoff[]";
    char Alpha_Baroczy[] = "Alpha_Baroczy[]";
    char Alpha_Chen[] = "Alpha_Chen[]";
    char Alpha_Chen86[] = "Alpha_Chen86[]";
    char Alpha_Chisholm[] = "Alpha_Chisholm[]";  
    char Alpha_Chisholm_II[] = "Alpha_Chisholm_II[]";
    char Alpha_Czop[] = "Alpha_Czop[]";
    char Alpha_Domanski[] = "Alpha_Domanski[]";
    char Alpha_Fauske[] = "Alpha_Fauske[]";
    char Alpha_Graham[] = "Alpha_Graham[]";
    char Alpha_Guzhov[] = "Alpha_Guzhov[]";  
    char Alpha_Hamersma[] = "Alpha_Hamersma[]";
    char Alpha_Harms[] = "Alpha_Harms[]";
    char Alpha_Homo[] = "Alpha_Homo[]";
    char Alpha_Huq[] = "Alpha_Huq[]";
    char Alpha_Kawahara[] = "Alpha_Kawahara[]";
    char Alpha_Kopke[] = "Alpha_Kopke[]";
    char Alpha_Laird[] = "Alpha_Laird[]";
    char Alpha_Lockhart[] = "Alpha_Lockhart[]";
    char Alpha_Madsen[] = "Alpha_Madsen[]";
    char Alpha_Massena[] = "Alpha_Massena[]";
    char Alpha_Nishino[] = "Alpha_Nishino[]";
    char Alpha_Premoli[] = "Alpha_Premoli[]";
    char Alpha_Rouhani[] = "Alpha_Rouhani[]";
    char Alpha_Smith[] = "Alpha_Smith[]";
    char Alpha_Spedding[] = "Alpha_Spedding[]";
    char Alpha_Tandon[] = "Alpha_Tandon[]";
    char Alpha_Thom[] = "Alpha_Thom[]";
    char Alpha_Turner[] = "Alpha_Turner[]";
    char Alpha_Wallis[] = "Alpha_Wallis[]"; 
    char Alpha_Xu[] = "Alpha_Xu[]";
    char Alpha_Yashar[] = "Alpha_Yashar[]";
    char Alpha_Zivi[] = "alpha_Zivi[]";
    
    

    double alpha_Armand[N]; //A.A. Armand, The resistance during the movement of a two-phase system in horizontal pipes, Izv. Vses. Tepl. Inst. 1 (1946) 16e23.
    double alpha_Bankoff[N]; // Bankoff, S.G., 1960. A variable density single fluid model for two phase flow with particular reference to steam water flow. Trans. ASME, J. Heat Transfer 82, 265–272.
	double alpha_Baroczy[N]; //C.J. Baroczy, Correlation of liquid fraction in two-phase flow with applications to liquid metals, Chem. Eng. Prog. Symp. Ser. 61 (1965) 179e191.	
	double alpha_Chen[N];  //J.J.J. Chen, P.L. Spedding, An extension of the LockharteMartinelli theory of two-phase pressure drop and holdup, Int. J. Multiphase Flow 7 (1981) 659e675.
    double alpha_Chen86[N]; // Chen, J.J.J., 1986. A further examination of void-fraction in annular two-phase flow. Int. J. Heat Mass Transfer 29, 1760–1763.
	double alpha_Chisholm[N];  //D. Chisholm, Pressure gradients due to friction during the flow of evaporating two-phase mixtures in smooth tubes and channels, Int. J. Heat Mass Transfer 16 (1973) 347e358.
	double alpha_Chisholm_II[N]; //D. Chisholm, Two Phase Flow in Pipelines and Heat Exchangers, Longman, New York, 1983.
	double alpha_Czop[N]; // Czop, V., Barbier, D., Dong, S., 1994. Pressure drop, void fraction and shear stress measurements in adiabatic two-phase flow in coiled tube. Nucl. Eng. Design 149, 323–333.
	double alpha_Domanski[N]; //P. Domanski, D. Didion, Computer Modeling of the Vapor Compression Cycle with Constant Flow Area Expansion Device, NBS Building Science Series 155, U.S. Department of Commerce & National Bureau of Standards, USA, 1983.
	double alpha_Fauske[N];  //H. Fauske, Critical two-phase, steamewater flows, in: Proceedings of 1961 Heat Transfer Fluid Mechanical Institute, Stanford University Press, California, 1961, pp. 79e89.
	double alpha_Graham[N]; //D.M. Graham, T.A. Newell, J.C. Chato, Experimental Investigation of Void Fraction during Refrigerant Condensation, ACRC TR-135, University of Illinois at Urbana-Champaign, USA, 1997.
	double alpha_Guzhov[N];  //A.L. Guzhov, V.A. Mamayev, G.E. Odishariya, A study of transportation in gas liquid systems, in: 10th International Gas Union Conference, 1967. Hamburg, Germany.
	double alpha_Hamersma[N]; // Hamersma, P.J., Hart, J., 1987. A pressure drop correlation for gas/liquid pipe flow with a small liquid holdup. Chem. Eng. Sci. 42, 1187–1196.
	double alpha_Harms[N]; // T.M. Harms, D. Li, E.A. Groll, J.E. Braun, A void fraction model for annular flow in horizontal tubes, Int. J. Heat Mass Transfer 46 (2003) 4051e4057.
	double alpha_Homo[N];  // Homogeneous model
	double alpha_Huq[N];  // R.H. Huq, J.L. Loth, Analytical two-phase flow void fraction prediction method, J. Thermophys 6 (1992) 139e144.
	double alpha_Kawahara[N]; // A. Kawahara, M. Sadatomi, K. Okayama, M. Kawaji, P.M.-Y. Chung, Effects of channel diameter and liquid properties on void fraction in adiabatic twophase flow through microchannels, Heat Transfer Eng. 26 (2005) 13e19.
	double alpha_Kopke[N]; // H.R. Kopke, T.A. Newell, J.C. Chato, Experimental Investigation of Void Fraction during Refrigerant Condensation in Horizontal Tubes, ACRC TR-142, University of Illinois at Urbana-Champaign, USA, 1998.
    double alpha_Laird[N]; // Chisholm, D., Laird, A.D.K., 1958. Two phase flow in rough tubes. Trans. ASME 80, 276–286.
	double alpha_Lockhart[N]; // R.W. Lockhart, R.C. Martinelli, Proposed correlation of data for isothermal two-phase, two-component flow in pipes, Chem. Eng. Prog. 45 (1949) 39e48.
	double alpha_Madsen[N]; // Madsen, N., 1975. A void fraction correlation for vertical and horizontal bulk-boiling of water. AIChE J. 21, 607–608. 
	double alpha_Massena[N]; // W.A. Massena, SteameWater Pressure Drop and Critical Discharge Flow e a Digital Computer Program, 1960. HW-65706.
	double alpha_Nishino[N]; // H. Nishino, Y. Yamazaki, A new method of evaluating steam volume fractions in boiling systems, J. Soc. Atom. Energy Japan 5 (1963) 39e59.
	double alpha_Premoli[N]; // A. Premoli, D. Francesco, A. Prina, A dimensionless correlation for determining the density of two-phase mixtures, La Termotecnica 25 (1971) 17e26.
	double alpha_Rouhani[N]; //  Rouhani, S.Z., Axelsson, E., 1970. Calculation of void volume fraction in the sub cooled and quality boiling regions. Int. J. Heat Mass Transfer 13, 383–393.
	double alpha_Smith[N]; // S.L. Smith, Void fractions in two phase flow: a correlation based upon an equal velocity head model, Proc. Inst. Mech. Eng. 36 (1969) 647e664.
	double alpha_Spedding[N]; // Spedding, P.L., Chen, J.J.J., 1984. Holdup in two phase flow. Int. J. Multiphase Flow 10, 307–339.
	double alpha_Tandon[N]; // Tandon, T.N., Varma, H.K., Gupta, C.P., 1985. A void fraction model for annular two-phase flow. Int. J. Heat Mass Transfer 28, 191–198.
	double alpha_Thom[N]; // J.R.S. Thom, Prediction of pressure drop during forced circulation boiling of water, Int. J. Heat Mass Transfer 7 (1964) 709e724.
	double alpha_Turner[N]; // J.M. Turner, G.B. Wallis, The Separate-cylinders Model of Two-phase Flow, NYO-3114-6, Thayer’s School Eng., Dartmouth College, Hanover, New Hampshire, USA, 1965.
	double alpha_Wallis[N]; // G.B. Wallis, One Dimensional Two-phase Flow, McGraw-Hill Inc., New York, 1969.
	double alpha_Xu[N]; // Y. Xu, X. Fang, Correlations of void fraction for two-phase refrigerant flow in pipes, Applied Thermal Eng. 64 (2014) 242e251. 
	double alpha_Yashar[N]; //D.A. Yashar, D.M. Graham, M.J. Wilson, J.C. Chato, H.R. Kopke, T.A. Newell, Investigation of refrigerant void fraction in horizontal, microfin tubes, HVAC&R Res. 7 (2001) 67e82.
	double alpha_Zivi[N]; //S.M. Zivi, Estimation of steady state steam void fraction by means of the principle of minimum entropy production, J. Heat Transfer 86 (1964) 247e252
	
		
		
	///***************************************************
	//// Parameters of Premoli et al. (1971)
	double E1 = 1.578*pow(Re_lo,-0.19)*pow(rho_l/rho_v,0.22);
	double E2 = 0.0273*We_lo*pow(Re_lo,-0.51)*pow(rho_l/rho_v,-0.08);
    double y[N];
	double S_Premoli[N];
		////****************************************************


  
   double S[N];

   x[0]=0.001;
   X_tt[0] = 0.001;
   alpha_Homo[0]=0.001;
   S[0]=1.0;

	for(i=1;i<N;i++)
	{
				
	x[i] = float(i)/1000+0.0001;	    
	  
	alpha_Homo[i] = 1/((1 + S[i-1]*(1-x[i])/(x[i])*(rho_v/rho_l)));
	  
	S[i] = (rho_l*x[i]*(1 - alpha_Homo[i]))/(rho_v*(1-x[i])*alpha_Homo[i]); 
	
	X_tt[i] = pow((1-x[i])/(x[i]+0.00001),0.9)*sqrt(rho_v/rho_l)*pow(mu_l/mu_v,0.1);
				
	}
	
	
    for(i=0;i<N;i++)
	{	
	    Ft[i] = sqrt( pow(G,2)*pow(x[i],3)/((1-x[i])*pow(rho_v,2)*g*D)); ///From Graham et al. (1997)
	    //	F_X_tt[i] =  0.15*(pow(X_tt[i],-1) + 2.85*pow(X_tt[i], -0.476)); ///From Tandor et al. (1985)
	    Re_l[i] = G*(1-x[i])*D/mu_l;
		S_Premoli[i] =  1 + E1*sqrt( (y[i]/(1+y[i]*E2)) - y[i]*E2 );

		y[i]= alpha_Homo[i]/(1.0-alpha_Homo[i]);	
		rho_tp[i] = pow((1-x[i])/rho_l + x[i]/rho_v, -1); ///From Guzhov et al. (1967)	
		Fr_tp[i] = pow(G,2)/(g*D*pow(rho_tp[i],2)); ///From Guzhov et al. (1967)
					
	}
	   
	  		  
	  for(i=0;i<N;i++)
	{
			  
	  alpha_Thom[i] = pow(( 1 + (1-x[i])/(x[i])*pow(rho_v/rho_l,0.89)*pow(mu_l/mu_v,0.18)),-1); 
	
	  alpha_Zivi[i] =  pow(( 1 + (1-x[i])/(x[i])*pow(rho_v/rho_l,0.6667)),-1); 
	  
	  alpha_Smith[i] =  pow(1+(1-x[i])/x[i]*(rho_v/rho_l)*(0.4 + (1 - 0.4)*sqrt( (rho_l/rho_v) + 0.4*((1-x[i])/x[i])/( 1+0.4*(1-x[i])/x[i]) ) ) , -1);     
	  
	  alpha_Premoli[i] =  pow((1 + (1-x[i])/(x[i])*(rho_v/rho_l)*S_Premoli[i]),-1);
	  	  
	  alpha_Fauske[i] = pow(( 1 + (1-x[i])/(x[i])*pow(rho_v/rho_l,0.5)),-1); 
	  
	  alpha_Chisholm[i] = pow( 1+ (1-x[i])/x[i]*(rho_v/rho_l)*sqrt(1-x[i]*(1-rho_l/rho_v)),-1); 
	  
	  alpha_Turner[i] = pow( 1 + pow( (1-x[i])/x[i] ,0.72)*pow( rho_v/rho_l,0.4)*pow(mu_l/mu_v,0.08),-1);
	
	  alpha_Lockhart[i] = pow( 1+0.28*pow(X_tt[i],0.71), -1);
	  
	  alpha_Baroczy[i] = pow( 1+pow((1-x[i])/(x[i]+0.0001), 0.74)*pow(rho_v/rho_l, 0.65)*pow(mu_l/mu_v,0.13),-1);
	   
	  alpha_Harms[i] = pow( 1.0 - 10.06*pow(Re_l[i],-0.875)*pow(1.74+0.104*sqrt(Re_l[i]), 2)*pow(1.376+7.242/pow(X_tt[i],1.655),-0.5),2);
	
	/////// ++++++++++++++  Domanski and Didion (1983)
	   if(X_tt[i]<=10){
	   	alpha_Domanski[i] = pow(1+pow(X_tt[i],0.8),-0.38);
	   }
	    else{
	    alpha_Domanski[i] =	0.823 - 0.157*log(X_tt[i]);
	    }
	/////// ++++++++++++++
	 
	   alpha_Graham[i]= 1 - exp( -1 - 0.3*log(Ft[i]) - 0.0328*pow(log(Ft[i]),2)); 
	
	   alpha_Yashar[i] = pow(1 + 1/Ft[i] + X_tt[i],-0.321);
	   
	   alpha_Wallis[i] = pow(1 + pow(X_tt[i],0.8),-0.38);
	   
	   alpha_Chen[i]= 3.5/(3.5 +pow(X_tt[i], 0.6667));
	   
	   alpha_Huq[i] = 1 - 2*pow(1-x[i],2)/( 1 - 2*x[i] + sqrt(1 + 4*x[i]*(1 - x[i])*(rho_l/rho_v - 1))); 
	
	   if(Ft[i]<0.044){
	   	alpha_Kopke[i] =  alpha_Homo[i];
	   }
	   else{
	   	alpha_Kopke[i] = 1.045 - exp(- 1 - 0.342*log(Ft[i]) - 0.0268*pow(log(Ft[i]),2) + 0.00597*pow(log(Ft[i]),3));
	   }
	   
	   alpha_Chisholm_II[i] = alpha_Homo[i]/(alpha_Homo[i] + sqrt(1 - alpha_Homo[i]));
	   
	   alpha_Armand[i] = 0.833*alpha_Homo[i];
	   
	   if(alpha_Homo[i]<0.9){
	    alpha_Massena[i] = 0.833*alpha_Homo[i];
	   }
	   else{
	   	alpha_Massena[i] = (0.833 + (1 - 0.833)*x[i])*alpha_Homo[i];
	   }
	   
	   alpha_Nishino[i] = 1 - sqrt((1-x[i])/x[i]*(rho_v/rho_l))*sqrt(alpha_Homo[i]);
	   
	   alpha_Guzhov[i] = 0.81*(1-exp(-2.2*sqrt(Fr_tp[i])))*alpha_Homo[i];	
	   
	  if(D == 0.0001){
	   alpha_Kawahara[i] =  0.03*sqrt(alpha_Homo[i])/(1-0.97*sqrt(alpha_Homo[i]));
	  }
	   if(D == 0.00005){
	   alpha_Kawahara[i] =  0.02*sqrt(alpha_Homo[i])/(1-0.98*sqrt(alpha_Homo[i]));
	  }
	   if(D > 0.00025){
        alpha_Kawahara[i] = 0.833*alpha_Homo[i];
	   }
	   
	alpha_Xu[i] = pow(1 + 2*pow(Fr_lo,-0.2)*pow(alpha_Homo[i],3.5)*(1-x[i])/x[i]*rho_v/rho_l,-1);
	
	alpha_Rouhani[i] = x[i]/rho_v*pow( (1 + 0.12*(1 - x[i]))*(x[i]/rho_v + (1 - x[i] )/rho_l) + 1.18*(1-x[i])*pow(g*sigma*(rho_l - rho_v)/(G*sqrt(rho_l)),0.25),-1);
	
     ////////+++++++ Madsen (1975) 
	b[i] = 1 + log(rho_l/rho_v)*pow(log((1-x[i])/x[i]),-1);
	alpha_Madsen[i] = pow(1 + pow((1-x[i])/x[i],b[i])*sqrt(rho_l/rho_v), -1);
	////////++++++++++++++++++++++++++++++
	
	alpha_Spedding[i] = pow(1 + 2.22*pow((1-x[i])/x[i],0.65)*pow(rho_v/rho_l,0.65), -1);
	
	alpha_Chen86[i] = pow( 1 + 0.18*pow( (1-x[i])/x[i] ,0.6)*pow(rho_v/rho_l,0.33)*pow(mu_l/mu_v,0.07),-1);
	
	alpha_Hamersma[i] = pow(1 + 0.26*pow((1-x[i])/x[i],0.67)*pow(rho_v/rho_l,0.33), -1);
	
    alpha_Bankoff[i] = (0.71 + 0.0145*P/1000)*alpha_Homo[i];
	
	alpha_Czop[i] = -0.285 + 1.097*alpha_Homo[i];
	
	alpha_Laird[i] =  1 + pow(0.8/(1 + 21/X_tt[i] + 1/pow(X_tt[i],2)),1.75);
	
	    if(Re_lo > 50 && Re_lo < 1125)
	    {
	    	alpha_Tandon[i] =  1 - 1.928*pow(Re_lo,-0.315)*pow(0.15*(1/X_tt[i] + 2.85/pow(X_tt[i],0.476)),-1) + 0.9293*pow(Re_lo,-0.63)*pow(0.15*(1/X_tt[i] + 2.85/pow(X_tt[i],0.476)),-2);
	    }
	     else
	     {
	     	alpha_Tandon[i] =  1 - 38*pow(Re_lo,-0.088)*pow(0.15*(1/X_tt[i] + 2.85/pow(X_tt[i],0.476)),-1) + 0.0361*pow(Re_lo,-0.176)*pow(0.15*(1/X_tt[i] + 2.85/pow(X_tt[i],0.476)),-2);
	     }
	
	}

	  
	  write = fopen("results.csv", "w");
	 
	 fprintf(write,"%s   %s   %s   %s   ", quality,Alpha_Homo, Alpha_Thom, Alpha_Smith );
     fprintf(write,"%s   %s   %s   %s   ", Alpha_Premoli, Alpha_Fauske, Alpha_Chisholm,Alpha_Turner );
     fprintf(write,"%s   %s   %s   %s   ", Alpha_Lockhart, Alpha_Baroczy, Alpha_Massena, Alpha_Bankoff );
     fprintf(write,"%s   %s   %s   %s   ", Alpha_Harms, Alpha_Domanski, Alpha_Graham, Alpha_Yashar);
     fprintf(write,"%s   %s   %s   %s   %s   %s   ", Alpha_Wallis, Alpha_Chen, Alpha_Huq, Alpha_Kopke, Alpha_Chisholm_II, Alpha_Armand );
     fprintf(write,"%s   %s   %s   %s   ", Alpha_Massena, Alpha_Nishino, Alpha_Guzhov, Alpha_Xu);
     fprintf(write,"%s   %s   %s   %s   %s   %s  ", Alpha_Rouhani, Alpha_Madsen, Alpha_Spedding, Alpha_Chen86, Alpha_Czop, Alpha_Tandon);
	 
	 
	  for(i=0;i<N;i++)
	{
	 fprintf(write,"\n%f   %f   %f   %f   ", x[i],alpha_Homo[i], alpha_Thom[i], alpha_Smith[i] );
	 fprintf(write,"%f   %f   %f   %f   ", alpha_Premoli[i], alpha_Fauske[i], alpha_Chisholm[i],alpha_Turner[i] );	
     fprintf(write,"%f   %f   %f   %f   ", alpha_Lockhart[i], alpha_Baroczy[i], alpha_Massena[i], alpha_Bankoff[i] );
     fprintf(write,"%f   %f   %f   %f   ", alpha_Harms[i], alpha_Domanski[i], alpha_Graham[i], alpha_Yashar[i]);
	 fprintf(write,"%f   %f   %f   %f   %f   %f   ", alpha_Wallis[i],  alpha_Chen[i], alpha_Huq[i], alpha_Kopke[i], alpha_Chisholm_II[i], alpha_Armand[i] );
	 fprintf(write,"%f   %f   %f   %f   ", alpha_Massena[i], alpha_Nishino[i], alpha_Guzhov[i], alpha_Xu[i]); 
	 fprintf(write,"%f   %f   %f   %f   %f   %f   ", alpha_Rouhani[i], alpha_Madsen[i], alpha_Spedding[i], alpha_Chen86[i], alpha_Czop[i], alpha_Tandon[i]);
	}
	 
	fclose(write);

	
	printf("\n\n ----------------All results can be found in results.csv.-------------\n");
	printf("-------------To format the file, follow the instructions in README.pdf.--------------\n\n");
	
	
return(0);	
}
