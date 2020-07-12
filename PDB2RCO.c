/* Progrma to calculate RCO for a PDB structure */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>


float get_cord(char In[], int a, int b);

int get_num(char In[], int a, int b);

int main (int argc, char *argv[])
{
	FILE *fpinp;
	FILE *fpout;
	FILE *fpout2;
	char inp[1000], tmp[15] ,In_pdb_file[200], Out_csv_file[200], Out_txt_file[200], at_name[50000][4]; 
	char Rsname[10000][10];
	char MRsname[100][10];
	char chain_ID[3];
	float X[10000],Y[10000],Z[10000],  t_int[4000], intcoe[500];
	int i, j, MResid[100], Resid[10000], multi_conf, Mcnt, cnt, cnt2, cnt3, CA_cnt[1000], file_CHK, dum;
	cnt=0;
	multi_conf = 0;
	Mcnt=0;
	cnt2=0;
	double dist=0;

	if(argc<5)
	{
		printf(" In correct usage\n");
		printf(" Correct usage: ./PDB2RCO file_ID chain_ID start_residue end_residue\n");
		printf(" Example: ./PDB2RCO 1an0 A 5 175 \n");
		exit(0);
	} 

	sprintf(In_pdb_file, "%s.pdb", argv[1]);
	sprintf(Out_csv_file, "%s.csv", argv[1]);
	sprintf(Out_txt_file, "%s.txt", argv[1]);

	file_CHK = access (In_pdb_file, F_OK);

	if ( file_CHK == 0 ) 
	{
                //printf("%s exists!!\n",In_pdb_file);
		dum=1;
        } else {
                printf("ERROR: %s doesn't exist!\n",In_pdb_file);
		exit(0);
        }



	for(i=0;i<14000;i++)
	{
		t_int[i]=0;
		intcoe[i]=0;
	}


	fpinp = fopen(In_pdb_file, "r");
	fpout = fopen(Out_csv_file, "w");
	fpout2 = fopen(Out_txt_file, "w");

//Get the coordinates and residue numbers
	while(fgets(inp,1000,fpinp) != NULL)
	{
		chain_ID[0]=inp[21];
                chain_ID[1]='\0';

		if(((inp[0]=='A' && inp[1]=='T' && inp[2]=='O') || (inp[0]=='H' && inp[1]=='E' &&inp[2]=='T') ) && strcmp(chain_ID, argv[2])==0)
		{
			for(i=17;i<20;i++)
				tmp[i-17]=inp[i];
			tmp[i-17]='\0';



			if( (strcmp(tmp, "ALA")==0 || strcmp(tmp, "ARG")==0 || strcmp(tmp, "ASN")==0 || strcmp(tmp, "ASP")==0 || strcmp(tmp, "CYS")==0 || strcmp(tmp, "GLU")==0 || strcmp(tmp, "GLN")==0 || strcmp(tmp, "GLY")==0 || strcmp(tmp, "HIS")==0 || strcmp(tmp, "ILE")==0 || strcmp(tmp, "LEU")==0 || strcmp(tmp, "LYS")==0 || strcmp(tmp, "MET")==0 || strcmp(tmp, "PHE")==0 || strcmp(tmp, "PRO")==0 || strcmp(tmp, "SER")==0 || strcmp(tmp, "THR")==0 || strcmp(tmp, "TRP")==0 || strcmp(tmp, "TYR")==0 || strcmp(tmp, "VAL")==0 || strcmp(tmp, "MSE")==0)  && inp[12]!='H' && inp[13]!='H' )
			{
				strcpy(Rsname[cnt], tmp);
					for(i=13;i<16;i++)
						tmp[i-13]=inp[i];
					tmp[i-13]='\0';
				strcpy(at_name[cnt], tmp);

				
				if(inp[16]=='A' && inp[13]=='C' && inp[14]=='A')
				{
					MResid[Mcnt]=get_num(inp,22,26);
					strcpy(MRsname[Mcnt],Rsname[cnt]);
					Mcnt++;
				}

				if(inp[13]=='C' && inp[14]=='A'&& (inp[16]=='A' || inp[16]==' '))
				{
					CA_cnt[cnt2]=get_num(inp,22,26);
					cnt2++;

				}
                
				if(inp[16]=='A' || inp[16]==' ')
				{
					Resid[cnt]=get_num(inp,22,26); //Get the residue number

					X[cnt]=get_cord(inp,26,38);

					Y[cnt]=get_cord(inp,38,46);

					Z[cnt]=get_cord(inp,46,54);

					cnt++;
				}
			}
		
		}
	}


	for(i=0;i<cnt2-2;i++)
	{
		if(CA_cnt[i+1]!=CA_cnt[i]+1)
		{
			printf(" There is chain break at %d and next residue is %d \n", CA_cnt[i], CA_cnt[i+1]);
			printf(" RCO for file %s.pdb not generated\n", argv[1]);
			exit(0);
		}
	}

// Calculate distance and grab the residue info	

	float rid2[10000],jj, RRct[10000];
	int  k;

	char Nresname[10000][20];
	int Nresn[10000];



	for(k=0; k<1000;k++)
	{
		RRct[k]=0;
		rid2[k]=0;
	}

	cnt3=0;
	int n1, n2;

        n1=atoi(argv[3]);
        n2=atoi(argv[4]);

		
	for(i=0;i<cnt;i++)
	{
		strcpy(Nresname[cnt3],Rsname[i]);

		Nresn[cnt3]=Resid[i];


		for(j=0;j<cnt;j++)
		{
			//if( Resid[i]!=Resid[j] && Resid[j]!= Resid[i]+1 && Resid[j]!= Resid[i]-1)
			if( Resid[i]!=Resid[j] && Resid[i]>n1 && Resid[i]<n2 && Resid[j]>n1 &&Resid[j]<n2)
			{

				dist=sqrt((X[j]-X[i])*(X[j]-X[i])+(Y[j]-Y[i])*(Y[j]-Y[i])+(Z[j]-Z[i])*(Z[j]-Z[i]));

				//printf("%5d %5d dist %5.2f\n", Resid[i], Resid[j], dist);

				if( dist > 2.0 && dist <= 3.5)
				{			
					//printf("%5d %5d \n", Resid[i], Resid[j]);
					//printf("%5d %5d %s  %s  dist %5.2f\n", Resid[i], Resid[j], at_name[i], at_name[j], dist);

					if(Resid[j] > Resid[i])
				 		jj= Resid[j]-Resid[i];	
					if(Resid[i] > Resid[j])
				 		jj= Resid[i]-Resid[j];	
						
					rid2[cnt3]=rid2[cnt3]+jj;
				
					RRct[cnt3]++;

				}
			}
		}

		if(Resid[i]!=Resid[i+1])
			cnt3++;
	}


	if(Mcnt>0)
	{
		printf("\nWARNING: Multiple conformations present for following residues. Only first conformation is considered for calculations\n");
	for(i=0;i<Mcnt;i++)
		if(MResid[i] >n1 && MResid[i]<n2)
			printf("Res. No=%5d Res. Name=%s\n", MResid[i],MRsname[i]);
	}
	

	double rco;

	for(i=0;i<cnt3;i++)
	{
			if(Nresn[i]>n1 && Nresn[i]<n2)
			{
				rco = rid2[i]/ (RRct[i]);
				fprintf(fpout2," %5s %5d  %10.3f\n", Nresname[i], Nresn[i], rco);
				fprintf(fpout,"%d,%.3f\n", Nresn[i], rco);
			}
	}
	

	fclose(fpinp);
	fclose(fpout);
	fclose(fpout2);


}



float get_cord(char In[], int a, int b)
{
	int i;
	char tmp[50];
	for(i=a;i<b;i++)
		tmp[i-a]=In[i];
	tmp[i-a]='\0';

	return(atof(tmp));
}

int get_num(char In[], int a, int b)
{
	int i;
	char tmp[50];
	for(i=a;i<b;i++)
		tmp[i-a]=In[i];
	tmp[i-a]='\0';

	return(atoi(tmp));
}

