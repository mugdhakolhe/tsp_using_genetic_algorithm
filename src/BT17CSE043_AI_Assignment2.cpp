/*
																				Assignment No. 2
																				
Use Genetic Algoritm to solve the Travelling Salesperson Problem (TSP on a large graph (greater than 50 nodes). 
=============================================================================================================================================================
Approach:

 In the following implementation, cities are taken as genes, string generated using these characters is called a chromosome, 
 while a fitness score which is equal to the path length of all the cities mentioned, is used to target a population.
 Fitness Score is defined as the length of the path described by the gene. Lesser the path length fitter is the gene. 
 The fittest of all the genes in the gene pool survive the population test and move to the next iteration. The number of iterations depends upon the value of a cooling variable. 
 The value of the cooling variable keeps on decreasing with each iteration and reaches a threshold after a certain number of iterations.

1. Initialize the population randomly.
2. Determine the fitness of the chromosome.
3. Until done repeat:
    1. Select parents.
    2. Perform crossover and mutation.
    3. Calculate the fitness of the new population.
    4. Append it to the gene pool
    
Assumptions:
1. In case graph is not fully connected put INT_MAX
2. Optimal values of variables
===============================================================================================================================================================

Made by

Author: Mugdha Satish Kolhe
Enrollment No: BT17CSE043
Course Code: CSL-412
Course: Artificial Intelligence
================================================================================================================================================================

*/


#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#define SIZE 75
#define START 1
#define NI 100
#define K 50
#define MR 0.5
#define CR 50

using namespace std;

vector<int> cost;

//function to print fitness
void print_fitness(vector <int> fit)
{
	for(int i=0; i<fit.size(); i++)
	{
		cout<<" fit "<<i<<": "<<fit[i]<<" ";
	}
	cout<<endl;
}

//function to print from a vector of vector
void print_pop(vector <vector < int > > pop)
{
	for(int i=0; i<K; i++)
	{
		for(int j=0; j<pop[i].size(); j++)
		{
			cout<<pop[i][j]<<" ";
		}
		cout<<endl;
	}
}

//function to print from vector of integers
void print(vector <int> vect)
{
	for(int i=0; i<vect.size(); i++)
	{
		cout<<vect[i]<<" ";
	}
	cout<<endl;
}

//function to print from vector of floats
void printf(vector <float> vect)
{
	for(int i=0; i<vect.size(); i++)
	{
		cout<<vect[i]<<" ";
	}
	cout<<endl;
}

//function to calculate fitness of each member from the population
int fitness_calc(vector<int> sol, vector< vector < int > > graph)
{
	int no, sum, c1, c2;
	sum=0;
	no=sol.size();
	
	//sum the distances between cities
	for(int i=0; i<no-1; i++)
	{
		c1=sol[i];
		c2=sol[i+1];
		sum+=graph[c1][c2];
	}
	return sum;
}

//function to create initial randomize population
vector<int>permute(int n)
{
	vector <int>result, temp;
	for(int i=0; i<n; i++)
	{
		if(i!=START)	//shuffle all the nodes except ths start
		{
			temp.push_back(i);
		}
	}
	int size=temp.size();
	for (int i = 0; i < size - 1; i++) //swaping
	{
      int j = i + rand() % (size - i);
      swap(temp[i], temp[j]);
   }
   
   
	result.push_back(START);	//push the start
	for(int i=0; i<size; i++)	//push the remaining nodes
	{
		result.push_back(temp[i]);
	}
	result.push_back(START);	//push the start node again for a round trip
	
	return result;
}

//function to create intial population
vector <vector < int> > create_pop(int n)
{
	vector <vector< int > > pop;
	vector <int > p;
	for(int i=0; i<K; i++)	//create K random permutations 
	{
		p=permute(n);
		pop.push_back(p);
	}
	return pop;
}

//select parent according to its fitnesss value
vector<int> random_selection (vector <vector <int > > population, vector <int> fitness)
{
 	//Roulette Wheel Selection
	vector <float> prob;
	float pop=0, t, sum=0;
	int no;
	for(int i=0; i<K; i++)
	{
		sum+=fitness[i];
	}
	//cout<<"sum: "<<sum<<endl;
	for(int i=0; i<K; i++)
	{
		t=(fitness[i]/sum)*100+pop;
		pop+=t;
		prob.push_back(t);
	}
	no=rand()%100;
	//printf(prob);
	//cout<<"n: "<<no;
	if(no<prob[0])
	{
		return population[0];
	}
	for(int i=1; i<K; i++)
	{
		if(no<prob[i] && no>=prob[i-1])
		{
			return population[i];
		}
	}
}

//Stochastic Universal Sampling as no negative values of fitness function
vector<int> sus (vector <int> x, vector <int> y)
{
	/*
	SUS(Population, N)
    F := total fitness of Population
    N := number of offspring to keep
    P := distance between the pointers (F/N)
    Start := random number between 0 and P
    Pointers := [Start + i*P | i in [0..(N-1)]]
    return RWS(Population,Pointers)

	RWS(Population, Points)
    Keep = []
    for P in Points
        i := 0
        while fitness sum of Population[0..i] < P
            i++
        add Population[i] to Keep
    return Keep
	*/
}

//create a child from two parents
vector<int> cross_over(vector<int>x, vector<int> y)
{
	//ordered crossover

	vector<int> result, temp, xtemp;
	int pt1, pt2, n, t, ptr;
	
	result.push_back(START);	//start remains same
	n=x.size()-2;
	pt1=(rand()%n)+1;	//create two pointers for making a subset
	pt2=(rand()%n)+1;
	
	if(pt2<pt1)
	{
		t=pt1;
		pt1=pt2;
		pt2=t;
	}
	
	for(int i=pt1; i<=pt2; i++)	//pushing the data from subset in y to temp
	{
		temp.push_back(y[i]);
	}
	
	for(int i=1; i<x.size()-1; i++)		//pushing the data from x not in subset in temp vector
	{
		if(count(temp.begin(), temp.end(), x[i])==0)
		{
			xtemp.push_back(x[i]);
		}
	}
	ptr=0;
	
	//creating a crossover child
	for(int i=1; i<pt1; i++)	
	{
		result.push_back(xtemp[ptr]);
		ptr++;
	}
	for(int i=pt1; i<=pt2; i++)	
	{
		result.push_back(y[i]);
	}
	for(int i=pt2+1; i<=n; i++)	
	{
		result.push_back(xtemp[ptr]);
		ptr++;
	}
	result.push_back(START);	//start remains same
	return result;
}

//function to swap a vector given two indexes
vector<int> swap(vector <int> s, int i, int j)
{
	vector<int> result;
	int k;
	for(int l=0; l<s.size(); l++)
	{
		k=l;
		if(l==i)
		{
			k=j;
		}
		if(l==j)
		{
			k=i;
		}
		result.push_back(s[k]);
	}
	return result;
}

//.function to create mutation according to given mutation rate
vector<int> mutation(vector<int> vect)
{
	vector<int> result;
	float rate;
	int r, n;
	n=vect.size()-1;
	
	result=vect;
	//swap
	for(int i=1; i<vect.size()-1; i++)
	{
		rate=((float)(rand()%100))/100;
		//cout<<"rate: "<<rate<<endl;
		if(rate<MR)
		{
			r=(rand()%(n-1))+1;
			//cout<<"i: "<<i<<" r: "<<r<<endl;
			result=swap(vect, i, r);
		}
	}
	return result;
}

//func in case you use temp in place of iterations
int cooldown(int temp) 
{ 
    return (90 * temp) / 100; 
} 

//Genetic Algorithm for TSP
vector<int> GeneticAlgoTSP(vector <vector < int > > graph)
{
	vector < vector < int > > population, new_population, best_population;
	vector <int> fitness;
	
	vector<int> child, x, y, best_child;
	int fit, ptr, min, n, temp, improvement;
	float r;
	
	n=graph.size();
	
	//creating the population
	population=create_pop(n);
	cout<<"Initial Population: "<<endl;
	print_pop(population);
	
	//applying GA
	int temperature=10000;
	int cutoff=10;
	for(int itr=0; itr<NI; itr++)	//while(temperature>1000)	//while(improvement<cutoff)
	{
		temperature=cooldown(temperature);
		for(int i=0; i<population.size(); i++)
		{
			fit=fitness_calc(population[i], graph);
			fitness.push_back(fit);
		}
		//print_fitness(fitness);
		new_population.clear();
		for(int i=0; i<K; i++)
		{
			//roulette wheel
			x=random_selection(population, fitness);
			y=random_selection(population, fitness);
			
			//stocastic universal sampling
			//child=(x, y);
			
			//evaluating the child
			//child=cross_over(x, y);
			
			//best child selection
			min=INT_MAX;
			for(int k=0; k<CR; k++)	
			{
				child=cross_over(x, y);
				temp=fitness_calc(child, graph);
				if(temp<min)
				{
					min=fitness_calc(child, graph);
					best_child=child;
				}
				/*else
				{
					float prob = pow(2.7, -1 * ((float)(temp - fitness_calc(x, graph)) / temperature)); 
                    if (prob > 0.5) 
					{ 
                        min=fitness_calc(child, graph);
                        best_child=child;
                	}
				}*/
			}
			child=best_child;
			/*cout<<"Parent 1: ";
			print(x); 
			cout<<"Parent 2: ";
			print(y); 
			cout<<"Child: ";
			print(child);*/
			
			//mutation
			r=((float)(rand()%100))/100;
			if(r<MR)
			{
				child=mutation(child);
			}
			new_population.push_back(child);
		}
		population=new_population;
		cout<<endl;
		cout<<"==============================================================================="<<endl;
		cout<<"Generation: "<<itr+1<<endl<<endl;
		//print_pop(population);
		
		//select best from this generation
		int prev=min;
		min=INT_MAX;
		for(int i=0; i<K; i++)
		{
			fit=fitness_calc(population[i], graph);
			if(fit<min)
			{
				min=fit;
				ptr=i;
			}

		}
		//improvement=prev-min;
		
		best_population.push_back(population[ptr]);
		cout<<"Best fit: "<<endl;
		print(population[ptr]);
		cout<<"Cost: "<<min<<endl;
		cost.push_back(min);
	}
	
	//return best from the population
	min=INT_MAX, ptr;
	for(int i=0; i<best_population.size(); i++)
	{
		fit=fitness_calc(best_population[i], graph);
		if(fit<min)
		{
			min=fit;
			ptr=i;
		}
	}
	vector<int> result=best_population[ptr];	//best of best population
	return result;
}

int main()
{
	
	int no, i, j, x, flag=0, min;	//no = vertices, START=Starting node
	string str, city, dist;
	ofstream fout;

//===============================================================================================================================================	
	//Uncomment this part of graph for  random graph generator
	
/*	
	int a[SIZE][SIZE];
	for(int i=0;i<SIZE;i++)	//randomized graph generator
	{
		for(int j=0;j<SIZE;j++)
		{
			int b=rand()%100+1;	//1 is added to prevent 0 from getting add as it would be meaningless
			if(i==j)
			{
				a[i][j]=0;	//same city to same city distance is 0
			}
			else if(b>0)	//symmetric graph
			{
				a[i][j]=a[j][i]=b;
			}
			else
			{
				j--;
			}
		}
	}
	string line; 
	fout.open("input3.txt");	//code to write graph generated to file
	fout<<SIZE<<endl;
	for(int i=0; i<SIZE; i++)
	{
		fout<<"city"<<i<<" ";
	} 
	fout<<endl;
	for(int i=0; i<SIZE; i++)
		{
			for(int j=0; j<SIZE; j++)
			{
				fout<<a[i][j]<<" ";
			}
			fout<<endl;
		}
	fout.close();*/
	
//===========================================================================================================================================================	
	
	 
	ifstream file("input3.txt");	//input file reading
	getline (file, str);
	stringstream ss(str); 
  	ss>>no;									//no stores the no of cities
  	cout<<"No of vertices: "<<no<<endl;;
	getline (file, str);
	string cities[no];
	stringstream iss(str);
	i=0; 
  	while (iss >> city) 
    {
    	cities[i]=city;
		i++;	
    }
    for(i=0; i<no; i++)	//cities array stores the cities
    {
    	cout<<cities[i]<<" ";
	}
	cout<<endl;
    i=0;
    vector < vector <int > > graph;	//graph stores the adjancy matrix of input distances
	while (getline (file, str)) 
	{
    	j=0; 
    	stringstream s(str);
    	vector<int> v;
  		while (s>>dist) 
    	{
    		stringstream geek(dist); 
  			geek>>x;   
    		//graph[i][j]=x;	
    		v.push_back(x);
    		//j++;
    	}
    	graph.push_back(v);
    	//i++;
  	}
  	for(i=0; i<no; i++)	//display graph
  	{
  		cout<<cities[i]<<": ";
  		for(j=0; j<no; j++)
  		{
  			cout<<graph[i][j]<<" ";
		}
		cout<<endl;
	}
	file.close();
	
	
	
	//genetic algorithm starts here

	vector<int> sol;
	int sum=0;
	fout.open("output3.txt");			//file in which output is stored
	
	fout<<"Population size = "<<K<<endl;
	fout<<"Mutation Rate = "<<MR<<endl;
	fout<<"No of children: = "<<CR<<endl;
	fout<<"No. of generations = "<<NI<<endl;
	
	fout<<endl;
	fout<<"========================================================================"<<endl<<endl;
	cout<<"Applying Genetic Algorithm: "<<endl;
	fout<<"Applying Genetic Algorithm: "<<endl;
	
	sol=GeneticAlgoTSP(graph);
	
	cout<<"========================================================================"<<endl;
	cout<<"The solution by GA is: "<<endl;
	fout<<"========================================================================"<<endl;
	fout<<"The solution by GA is: "<<endl;
	for(int i=0; i<sol.size()-1; i++)
	{
		cout<<sol[i]<<" -> ";
		fout<<sol[i]<<" -> ";
	}
	cout<<sol[i]<<endl;
	cout <<endl;
	cout<<"Cost of the solution found by Genetic Algorithm is: "<< fitness_calc(sol, graph)<<endl;
	fout<<sol[i]<<endl;
	fout <<endl;
	fout<<"Cost of the solution found by Genetic Algorithm is: "<< fitness_calc(sol, graph)<<endl;
	fout<<endl<<endl;
	for(int i=0; i<cost.size(); i++)
	{
		fout<<"Best cost in generation "<<i+1<<": "<<cost[i]<<endl;
	}
	fout.close();
	return 0;
}
