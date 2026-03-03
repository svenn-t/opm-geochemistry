import pandas as pd
import numpy as np
import os
import glob
import sys
import matplotlib.pyplot as plt
import itertools
import pathlib as pt
from scipy import stats

def get_col_data(col_name,unique_name,df):
    ''' returns data frame with column that matches unique_name'''
    try:
        return df.loc[df[col_name] == unique_name]
    except:
        print('Could not find data for '+ unique_name)
        return pd.DataFrame()

def replace_chars(name, chars=[" ", "/","Å", "Ø", "Æ"], new_chars=["_","_","AA","O","AE"]):
    ''' replace Norwegian characters and space in names'''
    new_name = name
    for ch,nch in zip(chars,new_chars):
        new_name = new_name.replace(ch, nch)
    return new_name

def get_columns_with(name,col_names):
    ''' returns list from col_names that contais name''' 
    l=[]
    for vi in col_names:
        if name in vi:
            l.append(vi)
    return l


def combine_images(list_names,no_col=3,no_row=5,dir_i='',dir_o = '../../fig/',pre_n='combined',plot_type='png'):
    if sys.platform == 'linux':
        exe='doconce'
    else:
        exe='ubuntu run /home/ah/anaconda3/bin/doconce'
    command = exe+' combine_images '+plot_type + ' -' +str(no_col) +str(' ')
    dd = list_names.copy()
    iter=0
    while(len(dd)>=no_row*no_col):
        files =""
        outfname= pre_n+str(iter)+'.'+plot_type
        
        outf2 = dir_o+outfname
        for i in range(no_row*no_col):
            files += str(dir_i)+str(dd.pop(0)) + " "
        print('11',command + files + str(outf2))
        
        os.system(command + files + str(outf2))
        iter += 1
    if(len(dd)>0):
        outfname= pre_n+str(iter)+'.'+plot_type
        outf = dir_o+outfname
        files = [str(dir_i+f) for f in dd]
        files = " ".join(files)
        print('22'+command + str(files) + ' ' + str(outf))
        os.system(command + str(files) + ' ' + str(outf))

def scatter_plot(x_data,y_data,labels,error_bar=[],xlim=None,ylim=None,fig_dir='../../fig/',
savefig='',dpi=200,trendline=False, xlabel='',ylabel='',skip=[],xscale=None,yscale=None,include_in_plot=True):
    '''
    plots x_data and y_data, possible to set limits, add error bars
    and add a regression (trend) line
    '''
    if xlim is not None:
        plt.xlim(xlim)
        xmin=min(xlim)
        xmax=max(xlim)
    else:
        xmin = -np.inf
        xmax = np.inf
    if ylim is not None:
        plt.ylim(ylim)
        ymin=min(ylim)
        ymax=max(ylim)
    else:
        ymin=-np.inf
        ymax=np.inf
    marker = itertools.cycle(('.',',','o','v','^','<','>','1','2','3','4','8',
                                's','p','P','*','h','H','+','x','X','D','d','|',
                                '_')) 
    color = itertools.cycle(('b','g','r','c','m','y','k'))
    x_new=[]
    y_new=[]
    x_plot=[]
    id=0
    for vi,vm,l in zip(x_data,y_data,labels):
        if l not in skip or include_in_plot:
            if xmin < vi < xmax and ymin < vm < ymax:
                x_plot.append(vi)
                if len(error_bar)==0:
                    plt.scatter(vi,vm,label=l,marker=next(marker),color=next(color))
                else:
                    plt.errorbar(vi,vm,error_bar[id],label=l,marker=next(marker),color=next(color))
                if l not in skip:
                    x_new.append(vi)
                    y_new.append(vm)
                else:
                    print('Skipping ', l, ' from linear fit')
            else:
                print('Skippng ', l,', x= ',vi,', y= ', vm)
            id+=1
    x_new=np.array(x_new)
    y_new=np.array(y_new)
    x_plot=np.array(x_plot)
    plt.grid()
    if trendline:
        slope, intercept, r_value, p_value, std_err= stats.linregress(x_new,y_new)
        line = slope*x_plot+intercept
        l1,=plt.plot(x_plot,line,'k:')
        if(intercept<0):
            l1_t2=str(format(intercept,'.3f'))
        else:
            l1_t2=r' + '+str(format(intercept,'.3f'))
        l1_t= r'$y$ = ' + str(format(slope,'.3f'))+ r'x '+l1_t2 +r' $, R^2$ = ' + str(format(r_value,'.3f'))
        legend1=plt.legend([l1], [l1_t], loc=2)
        plt.gca().add_artist(legend1)
        print('R-value=',r_value)

    plt.legend(loc='upper center', ncol=6, bbox_to_anchor=(0.5,-0.2),prop={'size': 7})
    plt.minorticks_on()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if xscale:
        plt.xscale(xscale)
    if yscale:
        plt.yscale(yscale)
    if not savefig =='':
        plt.savefig(fig_dir+savefig, bbox_inches='tight',transparent=True,dpi=dpi)
    plt.show()


def path_exists(name):
    f=pt.Path(name)
    if f.exists():
        return True
    else:
        print('file or folder ' + name + ' does not exists')
        return False

def part_of_list(list,elem):
    '''
    checks if elem contains any elements in list
    '''
    for li in list:
        if li in elem:
            return True
    return False


def insensitive_glob(pattern):
    def either(c):
        return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c
    return glob.glob(''.join(map(either, pattern)))

def find_files_of_type(type, dir='',name='',ignore=['~'],ignore_case=True):
    path=dir+'/'+'*'+name+'*'+type
    if ignore_case:
        list=insensitive_glob(path)
    else:
        list=glob.glob(path)
    new_list = [x for x in list if not part_of_list(ignore,x)]
    return new_list

def sort_lists(list1,list2,reverse=False):
    '''
    sort list1 low to high and list2 correspondingly
    '''
    l2=[x for _, x in sorted(zip(list1, list2),reverse=reverse)]
    l1=sorted(list1,reverse=reverse)
    return l1,l2

