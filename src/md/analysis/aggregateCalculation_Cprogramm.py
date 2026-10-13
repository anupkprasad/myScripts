import sys 
import subprocess 
import os 
 
prog = r''' 
// line no 143


#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;



void count(int natom, int fr, int totfr)
{

	int xsig{}, i = {}, j = {}, count=0, cnt=0, fin_count=0, id = {};      ////natom
	int nhyd =natom;
	double x = {}, y = {}, z = {}, sumx = {}, sumy = {}, sumz = {}, sum=0;
	char dump[1000] = {};
	FILE *file, *output, *list, *Nc, *per_frame, *bin;
	file = fopen("coordinate", "r");
	output = fopen("number", "w");
	list = fopen("list", "w");
	per_frame = fopen("cluster_frame", "a");
	bin = fopen("frame_cluster_bin", "a");
	fprintf(per_frame, "##No.	 	#cluster, frameNo.%d\n", fr);
	Nc = fopen("Nc_VS_frame", "a");
	int min = 1;
	printf("F.Enter the minimum number of atoms for critical clusters :%d\n", min);

//*****************************************************************************//
		double r_cut = 20, deltar = 1;
		int Nbin = r_cut/deltar;
		int BIN[Nbin];
		for (int i=0; i<Nbin; i++)
			 BIN[i]=0;
//*****************************************************************************//



	while (fgets(dump, 1000, file)!=NULL)
	{
		if (strlen(dump)==1)
			count++;
	}
	fclose(file);
	printf("F.cluster count = %d\n",count);
	file = fopen("coordinate", "r");
	int clus_cnt[count];
	for(i=0; i<count; i++)
		clus_cnt[i] = 0;
	while (fgets(dump, 1000, file)!=NULL)
	{
		if (strlen(dump)==1)
		{
			fprintf(list, "%5d%5d\n", (cnt+1), clus_cnt[cnt]);
			fprintf(per_frame, "%5d%5d\n", (cnt+1), clus_cnt[cnt]);
			cnt++;
		}
		else
		{
			fprintf(list, "%s", dump);
			clus_cnt[cnt]++;
		}
	}
	printf("cluster count in clust_cnt[] array\n");
	for (int l=0; l<count; l++)
		{
			printf("%d", clus_cnt[l]);
			printf("\n");
			// xsig=floor(clus_cnt[l]);
			// cout << "xsig is: " << xsig << endl;
			BIN[(clus_cnt[l]-1)]++;
		}
			printf("\n");
	fprintf(bin, "%d\t\t", fr+1);
	for (int m=0; m<Nbin; m++)
	{
		printf("%d\t",BIN[m]);
		fprintf(bin, "%d\t", BIN[m]);
	}
		fprintf(bin, "\n");
		printf("\n");


	fclose(file);
	fclose(list);
	fclose(per_frame);
	fclose(bin);
	file = fopen("list", "r");
	for (i=0; i<count; i++)
	{
		if (clus_cnt[i] >= min)
		{
			sumx=0;
			sumy=0;
			sumz=0;
			for (j=0; j<clus_cnt[i]; j++)
			{
				fscanf(file, "%lf%lf%lf\n", &x, &y, &z);
				sumx += x;
				sumy += y;
				sumz += z;
			}
			//fprintf (output, "%8d%8d%12.3f%12.3f%12.3f\n", i+1, clus_cnt[i], sumx, sumy, sumz);
			fprintf (output, "%5d%5d%8.3f%8.3f%8.3f\n", i+1, clus_cnt[i], sumx/(double)clus_cnt[i], sumy/(double)clus_cnt[i], sumz/(double)clus_cnt[i]);   //clsuter com
			fgets(dump, 1000, file);
			fin_count++;
		}
		else
		{
			for (j=0; j<clus_cnt[i]+1; j++)
				fgets(dump, 1000, file);
		}		
	}
	printf("F.fin_count = %d\n", fin_count);					// no of cluster
	fprintf(Nc,"%d 	%d\n", fr+1, fin_count);
	fclose(file);
	fclose(output);
	file = fopen("number", "r");	
	for (i=0; i<fin_count; i++)
	{
		fscanf(file, "%d%d%lf%lf%lf\n", &id, &cnt, &x, &y, &z);
		sum += cnt;
	}
	if (fin_count != 0)
		sum /= (double)fin_count;
	printf ("F.The Average Number of Clusters :%8.3f\n", sum);
	fclose(file);
	fclose(Nc);
}





int main()
{

	// FILE *filebarrel, *output, *coords, *AVG, *meshsize;
	fstream fout;
	fout.open("com_eachpep.dat", ios::in);
	FILE  *output, *coords, *AVG, *meshsize;
		// filebarrel = fopen("com_eachpep.dat", "r");                     
		int k, m = 0, p = 0, totfr = 5;
		int i, natom=4, nstep=1200, fr, type, index=0, j, lk;         /////nmonomer=7200    natom=12240  ///nstep=610
		// char id[20];
		string id;
		double cutoff = 12, DIST[totfr]{}, D;
		int nhyd = natom;
		int l[nhyd], ind[nhyd], flag[nhyd];
		printf("M.Enter the cutoff distance :%f\n", cutoff);
		double xs, ys, zs, xi, yi, zi, x[nhyd], y[nhyd], z[nhyd], xj, yj, zj, vx, vy, vz, dist, rxk, ryk, rzk;
		
		for (fr = 0; fr < nstep; fr++)
		{
			output = fopen("cluster_list", "w");
			coords = fopen("coordinate", "w");
			// printf("filebarrel pointer: %ld\n", ftell(filebarrel));
			char dump[1000] = {};
			string line[nhyd];
			int index = 0, l[nhyd] = {}, ind[nhyd] = {}, flag[nhyd] = {};
			double xs = {}, ys = {}, zs = {}, xi = {}, yi = {}, zi = {}, x[nhyd] = {}, y[nhyd] = {}, z[nhyd] = {}, xj = {}, yj = {}, zj = {}, vx = {}, vy = {}, vz = {}, dist = {}, rxk = {}, ryk = {}, rzk = {};
			for (i=0; i<1; i++)
			{
				// fgets(dump, 1000, filebarrel);
				fout >> id;
				cout << id << endl << endl;
			}
			for (i=0; i<natom; i++)         ////////replaced natom
			{
				// fscanf(filebarrel, "%s%lf%lf%lf\n", &id, &xs, &ys, &zs);
				fout >> id >> xs >> ys >> zs;
				// if (type == 2)
				{
					x[index] = xs;
					y[index] = ys;
					z[index] = zs;
					// ind[index] = id;
					line[index] = id;
					index++;
				}
			}

			printf("M.array size index= %d\n", index);
			//-------------
			for (i=0; i<nhyd; i++)
				l[i] = i;
			for (i=0; i<nhyd-1; i++)
			{
				if (i == l[i])
				{
					j=i;
					xj = x[j];
					yj = y[j];
					zj = z[j];
					for (k=i+1; k<nhyd; k++)
					{
						lk = l[k];
						if (lk == k)
						{
							rxk = xj-x[k];
							ryk = yj-y[k];
							rzk = zj-z[k];
							dist = sqrt(pow(rxk,2)+pow(ryk,2)+pow(rzk,2));
							if (dist <= cutoff)
							{
								l[k] = l[j];
								l[j] = lk;
							}
						}
					}
					j = l[j];
					xj = x[j];
					yj = y[j];
					zj = z[j];
					while (j != i)
					{
						for (k=i+1; k<nhyd; k++)
						{
							lk = l[k];
							if (lk == k)
							{
								rxk = xj-x[k];
								ryk = yj-y[k];
								rzk = zj-z[k];
								dist = sqrt(pow(rxk,2)+pow(ryk,2)+pow(rzk,2));
								if (dist <= cutoff)
								{
									l[k] = l[j];
									l[j] = lk;
								}
							}
						}
						j = l[j];
						xj = x[j];
						yj = y[j];
						zj = z[j];
					}
				}
			}
			//-------------
			for (i=0; i<nhyd; i++)
			{
				flag[i] = 0;
				fprintf(output, "%8d%8d\n", i+1, l[i]+1);
			}
			
			int lit;
			for (i=0; i<nhyd; i++)
			{
				if (flag[i] == 0)
				{
					lit = l[i];
					if (lit != i)
					{
						fprintf(coords, "%10.3f %10.3f %10.3f\n", x[i], y[i], z[i]);
						while (lit != i)
						{
							flag[lit] = 1;
							fprintf(coords, "%10.3f %10.3f %10.3f\n", x[lit], y[lit], z[lit]);
							lit = l[lit];
						}
						fprintf(coords, "\n");
					}
					else
						fprintf(coords, "%10.3f %10.3f %10.3f\n\n", x[i], y[i], z[i]);			
				}
			}
			fclose(output);
			fclose(coords);

			count(natom, fr, totfr); 						///storing 4node-avg for each single frame in D and then in array DIST[p] for all frames 
		}
		remove("cluster_list");
		remove("coordinate");
		remove("cluster_frame");
		remove("Nc_VS_frame");
		remove("list");
		remove("number");
	// fclose(filebarrel);
	fout.close();
	return 0;
}

''' 

mainpath = "/media/anup/Backup Plus/simulation/Reeju/analysis_anup/4uperin/"

filepathlist = ['0_nosalt/set1/','10_nosalt/set1/', '10_nosalt/set2/', '10_nosalt/set3/','15_nosalt/set1/', '15_nosalt/set2/', '15_nosalt/set3/',
               '20_nosalt/set1/', '20_nosalt/set2/', '20_nosalt/set3/', '40_nosalt/set1/', '40_nosalt/set2/', '40_nosalt/set3/',"40_withsalt/set1","40_withsalt/set2","40_withsalt/set3"]

for filepath in filepathlist:
    curdir = mainpath + filepath
    print(curdir)
    os.chdir(curdir)
    f = open('eachpep_aggre.cpp', 'w') 
    f.write(prog)
    f.close() 
    #subprocess.call(["gcc", "eachpep_aggre.cpp", "-ofoo", "-std=c99", '-w', '-Ofast'])
    subprocess.call(["g++", "eachpep_aggre.cpp"])
    subprocess.call("./a.out")
    
print("done")

