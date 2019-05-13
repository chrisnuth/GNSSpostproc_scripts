import os
import glob
import datetime as dt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from pyproj import Proj
from scipy import stats

def get_file_info(in_filestring):
    in_filename = os.path.basename(in_filestring)
    in_dir = os.path.dirname(in_filestring)
    if in_dir == '':
        in_dir = '.'
    return in_filename, in_dir

def read_rtklib_file(in_filename):
    """ Read in .pos file from RTKLib
    """
    print("reading file: ",in_filename)
    
    myProj = Proj("+proj=utm +zone=33X, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    
    mynames = ('UTCdate','UTChour','latitude','longitude','height','Q','ns','sdn','sde','sdu','sdne','sdeu','sdun','age','ratio')
    df = pd.read_csv(in_filename,  comment='%', delim_whitespace=True, names=mynames)
    df['timestamp']=df.UTCdate + ' ' + df.UTChour
    df.timestamp=pd.to_datetime(df.timestamp, format='%Y/%m/%d %H:%M:%S')
    df.set_index(df.timestamp, inplace=True)
    df.drop(columns=["UTCdate","UTChour"],inplace=True)
    
    df.rename(columns={"latitude": "Y","longitude": "X","height": "Z"}, inplace=True)
    df.X, df.Y = myProj(df.X.values,df.Y.values)
    
    
    return df
   

class GNSS(object):
    """Create a GNSS object from a RTKLib processed dataset.
   
    Parameters
    ----------
    in_filename : str or gdal object.
        If in_filename is a string, the GNSS is created by reading the file 
        corresponding to that filename. If in_filename is a gdal object, the 
        GeoImg is created by operating on the corresponding object.
    in_dir : str, optional
        Directory where in_filename is located. If not given, the directory
        will be determined from the input filename.
    dtype : numpy datatype
        numpy datatype to read input data as. Default is np.float32. See numpy docs for more details.
        
    """
    def __init__(self, in_file_dir):
        
        if os.path.isfile(in_file_dir):
            in_filename, in_dir = get_file_info(in_file_dir)
            self.filename = in_filename
            self.in_dir_path = in_dir
            self.df = read_rtklib_file(in_file_dir)
        elif os.path.isdir(in_file_dir):
            myfiles = glob.glob(os.path.join(in_file_dir,"*.pos"))
            self.filename = myfiles
            self.in_dir_path = in_file_dir
#            dict={}
            mynames = ('Y','X','Z','Q','ns','sdn','sde','sdu','sdne','sdeu','sdun','age','ratio')
#            for name in mynames:
                
            self.df = pd.DataFrame(columns=mynames)
            start_tot = time.time()
#            i=1
            data=list()
            for filename in myfiles:
                if os.stat(filename).st_size>5000:
                    #start = time.time()
#                    print(i + "fuck")
                    data.append(read_rtklib_file(os.path.join(self.in_dir_path,filename)))
#                    i=i+1
#                    self.df = self.df.append(data)
                    #end = time.time()
                    #print("Finished reading in :",end-start)
#            self.df = pd.DataFrame(data,columns=mynames)
            self.df = pd.concat(data,axis=0)                    
            end = time.time()
            print("TOTAL IMPORT TIME: ", end-start_tot)
        else:
            raise Exception('must be a string path to file or directory of .pos files')
       
        
    def filter_Q(self):
        '''
        Filter all solutions where ambiguities are not resolved. 
        '''
        self.df_raw = self.df
        
        n, _ = self.df.shape
        nf = sum(i != 1 for i in self.df.Q)
        psolve = 100-(nf/n)*100
        print ("Percent Amb Resolved: ", psolve," %")

        self.df = self.df[self.df.Q==1]
        
        return self
    
    def filter_st(self,thr_sdn=0.02,thr_sde=0.02,thr_sdu=0.04):
        '''
        Filter data based on std deviations of raw positions
        '''
        self.df = self.df[(self.df.sdn<thr_sdn) & (self.df.sde<thr_sde) & (self.df.sdu<thr_sdu)]
        return self

    def filter_pos_3std(self,grpcode='H'):
        '''
        
        '''
        def myfilter(df):
            df = df[(np.abs(stats.zscore(df.X)) < 3)]
            df = df[(np.abs(stats.zscore(df.Y)) < 3)]
            df = df[(np.abs(stats.zscore(df.Z)) < 3)]
            return df
        
#        self.df = self.df.groupby(pd.Grouper(freq='M')).apply(myfilter)
        self.df = self.df.groupby(pd.Grouper(freq=grpcode)).apply(myfilter)
        self.df.reset_index(drop=True,inplace=True)
        self.df.set_index(self.df.timestamp,inplace=True)
        return self
    
    def pos_filter(self,grpcode='H'):
        '''
        Positional filter: Checks horizontal and vertical distance
        to the hourly [default] average. horizontal distances outside 
        the 1 std dev are removed, then euclidean distances greater
        than 2 standard deviations are removed. 
        '''
        data=self.df
        def myfilter(df):
            dx = df.X - df.X.mean()
            dy = df.Y - df.Y.mean()
            dz = df.Z - df.Z.mean()

            euc_d = np.sqrt( dx**2 + dy**2 + dz**2) 
            hor_d = np.sqrt( dx**2 + dy**2 ) 
            zsc_e = stats.zscore(euc_d)
            zsc_d = stats.zscore(hor_d)
            
#            plt.figure()
#            plt.scatter(df.X,df.Y,c=euc_d)
#            plt.figure()
#            plt.scatter(df.X,df.Y,c=zsc_d)
#            plt.figure()
#            plt.hist(zsc_d,bins=100)
#            plt.scatter(df.X,df.Y,c=zsc)
            
            df = df[((np.abs(zsc_d) < 1.5) | (np.abs(zsc_e) < 2))]
            return df
        
#        self.df = self.df.groupby(pd.Grouper(freq='M')).apply(myfilter)
        data1 = data.groupby(pd.Grouper(freq=grpcode)).apply(myfilter)
        data1.reset_index(drop=True,inplace=True)
        data1.set_index(data1.timestamp,inplace=True)
        
        plt.figure()
        plt.plot(data.X,data.Y,'r.')
        plt.plot(data1.X,data1.Y,'k.')

        self.df = data1        
        return self
             
    
    def hourly_ave(self):
        '''
        Replace dataframe with hourly averages
        '''
        self.df_hour = self.df.resample('H').mean()
        return(self)

    def daily_ave(self):
        '''
        Replace dataframe with hourly averages
        '''
        self.df_day = self.df.resample('D').mean()
        return(self)
    
    def velocity(self):
        '''
        calculate velocity in meters per day
        '''
        df = self.df_day
        df['Time'] = df.index.asi8
        dist = df.diff().fillna(0.)
        dist['Dist'] = np.sqrt(dist.X**2 + dist.Y**2)
        dist['Speed'] = dist.Dist / ((dist.Time /1E9) / 86400)

        plt.figure()
        dist['Speed'].plot()
        
        
        self.vel = dist
        return self
    
    def plot_xyztime(self,ptype='original',daterange=None):
        '''
        Plot time series of positions and statistics
        '''

        if daterange is None:
            df = self.df
            if hasattr(self,"df_hour"): df_hour = self.df_hour
            if hasattr(self,"df_day"): df_day = self.df_day
        else:
            df = self.df[daterange[0]:daterange[1]]
            if hasattr(self,"df_hour"): df_hour = self.df_hour[daterange[0]:daterange[1]]
            if hasattr(self,"df_day"): df_day = self.df_day[daterange[0]:daterange[1]]
        
        if hasattr(self,"df_hour"):
            pos_df = df_hour.iloc[:,:3]
            stats_df = df_hour.iloc[:,5:8]
            pos_df.plot(subplots=True)
            stats_df.plot(subplots=True)
        else:
            pos_df = df.iloc[:,:3]
            stats_df = df.iloc[:,5:8]
#            plt.figure()
            pos_df.plot(subplots=True)
#            plt.figure()
            stats_df.plot(subplots=True)
                

    def plot_xy(self,daterange=None):
        '''
        Plot the cartographic displacements through time
        daterange: is a tuple with string dates e.g. YYYY-MM-DD
        for the start and stop dates to plot
        '''
        
        if daterange is None:
            df = self.df
            if hasattr(self,"df_hour"): df_hour = self.df_hour
            if hasattr(self,"df_day"): df_day = self.df_day
        else:
            df = self.df[daterange[0]:daterange[1]]
            if hasattr(self,"df_hour"): df_hour = self.df_hour[daterange[0]:daterange[1]]
            if hasattr(self,"df_day"): df_day = self.df_day[daterange[0]:daterange[1]]
        
            
        
        fig = plt.figure()
        plt.title('Positions through time')
        if hasattr(self,"df_hour") & hasattr(self,"df_day"):
            plt.plot(df.X,df.Y,'.',color='k',ms=1)
            plt.plot(df_hour.X,df_hour.Y,'o:',color='b',ms=5)
            plt.plot(df_day.X,df_day.Y,'s-',color='r',ms=8)
            plt.legend(['Raw','Hour','Day'])
        elif hasattr(self,"df") and not hasattr(self,"df_hour") and hasattr(self,"df_day"):
            plt.plot(df.X,df.Y,'.',color='k')
            plt.plot(df_day.X,df_day.Y,'s-',color='r',ms=3)
            plt.legend(['Raw','Hour'])
        elif hasattr(self,"df") and not hasattr(self,"df_hour") and not hasattr(self,"df_day"):
            plt.plot(df.X,df.Y,'.',color='k')
            plt.legend(['Raw'])
        else:
            print("dataframes not available")
#            raise TypeError('input data must hourly and daily dataframes.')
        
        
        