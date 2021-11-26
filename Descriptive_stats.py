import numpy as np
from pandas.core import arrays
from TumourCA import run
from main import get_kwargs, parallel_run
import csv
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, MultipleLocator
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split,cross_val_score
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.svm import SVR

n_iter=3
mort=0.5
agent_ratio=0

def write_csv(file,results):
    with open(file,'a') as f:
        write=csv.writer(f)
        write.writerow(results)

def evaluate_data(output_file):
    # write_csv(output_file,["mot","mul","mort","agent_ratio","mean_tot","std_tot","mean_half","std_half"])
    i=0
    agent_ratio=0
    mul=0.5
    for mot in np.linspace(0,1,5):
        for mort in np.linspace(0,1,5):
            # for mort in np.linspace(0,1,5):
            if i%5==0:
                print("%s/%s complete"%(str(i),25))
            kwargs,_ = get_kwargs(mot,mul,mort,agent_ratio)
            times,_,_,cellcounts = parallel_run(run,kwargs,n_iter)
            mean_counts = np.mean(cellcounts,axis=0)
            std = np.std(cellcounts,axis=0)
            t_half = int(len(times)/2)
            mean_half = mean_counts[t_half]
            std_half = std[t_half]
            mean_tot=mean_counts[-1]
            std_tot=std[-1]
            results=[mot,mul,mort,agent_ratio,mean_tot,std_tot,mean_half,std_half]
            write_csv('results/results.csv',results)
            i+=1

def load_csv(file):
    out=[]
    with open(file) as csvfile:
        f = csv.reader(csvfile, delimiter=',')
        for row in f:

            out.append(row)
    return out

def visualise(resultsfile):
    results = load_csv(resultsfile)
    df = pd.DataFrame(results[1:],columns=results[0])
    mot,mul,mort,z ='mot','mul','mort','std_tot'
    
    mot_vals,mul_vals,mort_vals=[set(df[x].values) for x in ['mot','mul','mort']] # Get unique values
    
    dep_var1=mul
    dep_var2=mort
    indep_var = z

    if (dep_var1==mot) and (dep_var2==mul):
        df = df.iloc[:25,:] # Mot and mul vary, mort constant
        series = mul_vals
        series_name = "Growth"
        xlabel='Motility'
    elif (dep_var1==mul) and (dep_var2==mort):
        df = df.iloc[25:50,:] # Mort and mul vary, mot constant
        series=mort_vals
        xlabel='Growth'
        series_name = "Mortality"
    elif (dep_var1==mort) and (dep_var2==mot):
        df = df.iloc[50:,:] # Mot and mort vary, mul constant
        series=mot_vals
        xlabel="Mortality"
        series_name = "Motility"
    else:
        raise ValueError("Incorrect inputs")

    for val in sorted(series):
        df_sub = df[(df[dep_var2]==val)]
        x = df_sub[dep_var1].values
        y = df_sub[indep_var].values
        plt.plot(x,y,label='%s: %s'%(series_name,val))
    
    plt.gca().yaxis.set_major_locator(MultipleLocator(5))
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) 
    plt.xlabel(xlabel)
    plt.ylabel(z)
    plt.legend()
    plt.show()


if __name__=="__main__":
    resultsfile='results/results.csv'

    # evaluate_data(resultsfile)       # Generate data with gridsearch
    visualise(resultsfile)















