import numpy as np
import pandas as pd
import time
import sys
import zipfile
import os
from pathlib import Path
from os import path
import _pickle as cPickle
import datetime
from datetime import timedelta



class Base():
    @staticmethod
    def  file_date(path_to_file):
        '''
        returns a files time stamp as a datetime object
        '''
        stat = os.stat(path_to_file)
        date = datetime.datetime.fromtimestamp(stat.st_mtime)
        return date, stat

    @staticmethod
    def  file_parts(file):
        '''
        splits a string file to path, file_name, file_type
        '''
        head_tail = os.path.split(file)
        path = head_tail[0]
        file_name, file_type = head_tail[1].split('.')
        return path, file_name, file_type

    @staticmethod
    def unfold_summary(df):
        for date in df.datetimeObserved.unique():
            idx = df.datetimeObserved==date
            vq = df.iloc[idx].valueQuantity.values
            df.valueDatetimeStart.iloc[idx] = df.loci[idx].datetimeObserved -(np.cumsum(vq[::-1])[::-1]).astype("timedelta64[s]")
            df.valueDatetimeEnd.iloc[idx] = df.valueDatetimeStart.iloc[idx] + vq.astype("timedelta64[s]")           
        return df

    @staticmethod
    def  remap_cat(_cat,_mapping,_df):
        _df[_cat] = pd.Categorical(_df[_cat])
        _df[_cat] = _df[_cat].cat.rename_categories(_mapping)
        return _df

    @staticmethod
    def  pickle_exist(pickle_file):
        return os.path.exists(pickle_file)

    @staticmethod
    def  map(value, key):
        return pd.Series(value, index=key).to_dict()

    
    def set_attributes(self):
        # self.clinical = []    # TODO: ask magda or helen to get the clinical data  
        self.domains  =   ["flags",
                            "sleep",
                            "physiology",
                            "wellbeing",
                            "location",
                            "demographics",
                            "movement",
                            "activity",
                            "clinical",
                            "doors",
                            "appliances",
                            "temperature",
                            "light",
                            "sleep_disturbance"]
        
        for domain in self.domains:
            setattr(self,domain,pd.DataFrame())
            
    def  find_duplicate(self,L):
        '''
        identifies duplicates in a list and returns their index
        '''
        seen, duplicate = set(), set()
        index = np.zeros(len(L), dtype=bool)
        seen_add, duplicate_add = seen.add, duplicate.add
        for idx, item in enumerate(L):
            if item in seen:
                duplicate_add(item)
                index[idx] = True
            else:
                seen_add(item)

        return self. ismember(L, list(duplicate)) != None, duplicate

    @staticmethod
    def  ismember(a, b):
        '''
        mimic's Matlabs ismemeber function (should be removed in later versions)
        '''
        bind = {}
        for i, elt in enumerate(b):
            if elt not in bind:
                bind[elt] = i
        return np.array([bind.get(itm, None) for itm in a])


    def  save_pickle(self):
        if not os.path.exists(self.pickle_file):
            filepath, _, _ = self. file_parts(self.pickle_file)
            Path(filepath).mkdir(parents=True, exist_ok=True)
        output_pickle = open(self.pickle_file, "wb")
        cPickle.dump(self, output_pickle)
        output_pickle.close()

    
    def  load_pickle(self):
        input_pickle = open(self.pickle_file, 'rb')
        data = cPickle.load(input_pickle)
        for domain in self.domains:
            setattr(self,domain,getattr(data,domain))
        input_pickle.close() 

class LoadData(Base):
    
    def __init__(self, input_path, output_path, datasets=['10','15','dri'], reload_data=True,reload_subset=True,verbose=True,domains=True):
        if verbose:tic = time.perf_counter()
        self.set_attributes()
        self.pickle_file = f"{output_path}/pkl/merged_tihm.pkl"                
        if self.pickle_exist(self.pickle_file) & reload_data:
            self.load_pickle()
            if verbose:
                print(f"Loading previously parsed merged file in {time.perf_counter()-tic:0.2f} seconds")
        else:        
            data_tmp = {}
            if len(datasets)>1:
                for dataset in datasets:
                    data_tmp[dataset] = PreProcess(dataset, input_path,output_path, reload_subset, verbose)
                for domain in  self.domains:  
                    _tmp = pd.concat([getattr(data_tmp[dataset],domain) for dataset in datasets])
                    _tmp = _tmp.drop_duplicates().reset_index(drop=True)
                    setattr(self, domain, _tmp)
            else:
                _tmp = PreProcess(datasets[0], input_path, output_path, reload_subset, verbose)
                for domain in  self.domains:  
                    setattr(self, domain, getattr(_tmp,domain))
            self.save_pickle()



class PreProcess(Base):
    tic = time.perf_counter()
    def __init__(self, dataset, input_path,output_path, reload_data=False,verbose=True):
        
        self.set_attributes()
        self.name = dataset
        if verbose:tic = time.perf_counter()
        file = str(list(Path(input_path).rglob(f'tihm{self.name}.zip'))[0])   
        date, _ = self.file_date(file)
        self.date = date.strftime("%Y%m%d")
        self.pickle_file = f"{output_path}/pkl/tihm_{dataset}_{self.date}.pkl"

        if self.pickle_exist(self.pickle_file) & reload_data:
            self.load_pickle()
            if verbose:
                print(f"Loading previously parsed TIHM {self.name} project in {time.perf_counter()-tic:0.2f} seconds")
        else:
            self.__parse_domains(file, verbose)
            if verbose:
                print(f"Total elapsed parsing time of TIHM {self.name} project: {time.perf_counter()-self.tic:0.2f} seconds")

      
    def __parse_domains(self, file, verbose):
        if verbose:
            tic = time.perf_counter()
        _zip = zipfile.ZipFile(file)
        _pid = pd.read_csv(_zip.open('Patients.csv'))
        _mapping = self.map(_pid.sabpId.values, _pid.subjectId.values)
        _pid['project'] = self.name
        self.demographics = _pid
        self.__parse_observations(_mapping, _zip, verbose)
        self.__parse_flags(_mapping, _zip, verbose)
        self.__parse_wellbeing(_mapping, _zip, verbose)

    def __parse_observations(self, _mapping, _zip, verbose):
        if verbose:tic = time.perf_counter()
        _df = pd.read_csv(_zip.open('Observations.csv'),
                        encoding='unicode_escape', low_memory=False)
        _type = pd.read_csv(_zip.open('Observation-type.csv'))
        _df.type = pd.Categorical(_df.type)
        if verbose:
            print(f"Loading TIHM {self.name} csv file took : {time.perf_counter()-tic:0.2f} seconds")
            tic = time.perf_counter()
            # subset any event that is logged but contains no actual values to a Dataframe under null key
        idx_null = _df[["valueBoolean", "valueState", "valueQuantity",
                    "valueDatetimeStart", "valueDatetimeEnd"]].isnull().values.all(axis=1)
        _df = _df[idx_null == False]
        _df['datetimeObserved'] = pd.to_datetime(_df['datetimeObserved'])
        _df = _df[_df['datetimeObserved'].dt.year > 2014]
        _df = _df[_df.type != "724061007"] # filter out device status
        _df = _df.sort_values(by=['datetimeObserved']).reset_index(drop=True)
        _df['activity'] = 1
        _df['project_id'] = _df.subject
        _df['display'] = _df.type
        _df = self.remap_cat('display', self.map(_type.display.values, _type.code.values), _df)
        _df = self.remap_cat('subject',_mapping, _df)
        _df['project'] = self.name
        _df['project'] = pd.Categorical(_df['project'])
        _df['subject'] = pd.Categorical(_df['subject'])
        self.__parse_location(_df, verbose)
        self.__parse_activity(_df, verbose)
        self.__parse_sleep(_df, verbose)
        self.__parse_physiology(_df, verbose)
        if verbose:
            print(f"Processing TIHM {self.name} observations domain took : {time.perf_counter()-tic:0.2f} seconds")
            tic = time.perf_counter()
        self.save_pickle()
        if verbose:
            print(f"Saving pickle: {time.perf_counter()-tic:0.2f} seconds")

    def __parse_activity(self,  _df, verbose):
        # convert movement, doors activity and appliences activity to a cleaned dataframe in day, hour and raw frequencies
        doors = []
        for k,subset in _df[_df.display=='Door'].groupby(['subject','location','project']):
            subset = subset[["datetimeObserved","valueState"]].pivot(columns="valueState",values="datetimeObserved").reset_index()
            if subset.shape[0] > 1:
                idx = np.where(subset.Open.isnull())
                open = subset.iloc[idx[0]-1].Open.values
                close = subset.iloc[idx[0]].Close.values
                delta = (close - open).astype('timedelta64[s]')
                m=open.shape[0]
                doors.append(pd.DataFrame({'project':[k[2]]*m,'subject':[k[0]]*m,'location':[k[1]]*m,'datetimeObserved':open
                                           ,"Close":close,"delta":delta,"activity":[1]*m}))

        doors = pd.concat(doors)
        doors = doors[~doors.location.isin(['B','Bathroom','C','Dining Room'])]
        doors = doors[doors.delta < timedelta(minutes=15)].reset_index(drop=True)
        self.doors = doors
        self.appliances = self.appliances[~self.appliances.location.isin(['A','B'])]
        self.appliances.location[self.appliances.location.isin(['Microwave','Toaster'])] = 'Oven'
        self.movement = self.movement[~self.movement.location.isin(['D','Study','Living Room',"Front Door"])]
        self.activity = pd.concat([self.doors,self.movement,self.appliances])[['project','subject','datetimeObserved','location']]
        self.activity = pd.get_dummies(self.activity, columns=['location'], prefix='', prefix_sep='').copy()

    def __parse_location(self,  _df, verbose):
        if verbose:tic = time.perf_counter()
        _df = _df[_df.location.notnull()].drop(columns=["datetimeReceived", "provider","valueUnit", "valueDatetimeStart", "valueDatetimeEnd"])
        if verbose:
            print(f"Processing TIHM {self.name} movement sensors : {time.perf_counter()-tic:0.2f} seconds")  
        self.location = _df
        self.appliances =_df[_df.display=='Does turn on domestic appliance']\
                            [["project","subject","datetimeObserved","location","activity"]].copy()
        self.movement = _df[_df.display =='Movement']\
                           [["project","subject","datetimeObserved","location","activity"]].copy()
        self.temperature = _df[(_df.valueQuantity.notnull())&(_df.display == "Room temperature")]\
                          [["project","subject","datetimeObserved","location","valueQuantity"]].copy()       
        self.temperature = self.temperature[~self.temperature.location.isin(["Living Room","Study"])].copy()                  
        self.light = _df[(_df.valueQuantity.notnull())&(_df.display == "Light")]\
                          [["project","subject","datetimeObserved","location","valueQuantity"]].copy()
        self.light = self.light[~self.light.location.isin(["Living Room"])]                  
        if verbose:
            print(f"Processing TIHM {self.name} location sensors : {time.perf_counter()-tic:0.2f} seconds")                                      
    
    def __parse_sleep(self,  _df, verbose):
        if verbose:tic = time.perf_counter()
        # TODO: compare the DRI and 1.5 data 
        self.sleep_disturbance = _df[_df.type.isin(['67233009'])].drop(columns=["datetimeReceived", "provider","type","device",
                                                                                "location", "valueBoolean", "valueState","provider",
                                                                                "device","valueDatetimeStart",
                                                                                "valueDatetimeEnd","activity"]).copy()
        _df = _df[_df.type.isin(['258158006', '29373008', '248218005', '60984000', '89129007', '421355008','307155000'] )].copy()
        _df['valueDatetimeStart'] = pd.to_datetime(_df['valueDatetimeStart'])
        _df['valueDatetimeEnd'] = pd.to_datetime(_df['valueDatetimeEnd'])
        
        _df['Start_End_logged'] = True
        _df['Start_End_logged'][_df.valueDatetimeStart.isnull()]= False
        _df.reset_index()
        if verbose:
            print(f"Processing TIHM {self.name} sleep_disturbances domain : {time.perf_counter()-tic:0.2f} seconds")         
        idx = _df['Start_End_logged']==True
        _df.valueQuantity.loc[idx]  = (_df[idx].valueDatetimeEnd - _df[idx].valueDatetimeStart).dt.seconds

        self.sleep = _df.drop(columns=["datetimeReceived","device", "provider", "location", "valueBoolean", "valueState","valueUnit","activity"])
        if verbose:
            print(f"Processing TIHM {self.name} sleep domain : {time.perf_counter()-tic:0.2f} seconds")         

    def __parse_physiology(self,  _df, verbose):
        if verbose:tic = time.perf_counter()
        _df = _df[_df.type.isin(
            ['8310-5', '55284-4', '29463-7', '251837008', '163636005', '248362003', '8462-4', '8480-6', '150456'])]
        self.physiology = _df.drop(columns=["datetimeReceived", "provider", "location", "valueBoolean",
                                            "valueState", "valueDatetimeStart", "valueDatetimeEnd"])
        if verbose:
            print(f"Processing TIHM {self.name} physiology domain : {time.perf_counter()-tic:0.2f} seconds") 

    def __parse_flags(self, _mapping, _zip, verbose):
        if verbose:tic = time.perf_counter()
        _df = pd.read_csv(_zip.open('Flags.csv'))
        _type = pd.read_csv(_zip.open('Flag-type.csv'))
        _cat = pd.read_csv(_zip.open('Flag-category.csv'))
        _val = pd.read_csv(_zip.open('FlagValidations.csv'))
        _df = pd.merge(_df, _val, how='outer', on=None,
                    left_on="flagId", right_on="flag",
                    suffixes=('_df', '_val'), copy=True)
        _df.category = pd.Categorical(_df.category)
        _df.rename(columns={'subject_df': 'subject'}, inplace=True)
        _df['project_id'] = _df.subject
        _df = self.remap_cat('subject', _mapping, _df)
        _df = self.remap_cat('category', self.map(_cat.display.values, _cat.code.values), _df)
        idx = self.find_duplicate(_type.display.values)[0]
        if any(idx):
            values = list(_type.code.values[idx])
            key = list(np.where(idx)[0])
            for key, val in dict(zip(key[0:-1], values[0:-1])).items():
                _df.type[_df.type == val] = values[-1]
                _type = _type.drop(key)
        _df = self.remap_cat('type', self.map(_type.display.values, _type.code.values), _df)
        _df['project'] = self.name
        _df['project'] = pd.Categorical(_df.project)
        self.flags = _df
        self.save_pickle()
        if verbose:
            print(f"Processing TIHM {self.name} flags domain : {time.perf_counter()-tic:0.2f} seconds") 

    def __parse_wellbeing(self, _mapping, _zip, verbose):
        if verbose:tic = time.perf_counter()
        _df = pd.read_csv(_zip.open('QuestionnaireResponses.csv'))
        _df['datetimeAnswered'] = pd.to_datetime(_df['datetimeAnswered'])
        _df = _df.sort_values(by=['datetimeAnswered'])
        _df = _df.drop(columns=["questionnaire", "datetimeReceived"])
        _df.question, questions = pd.factorize(_df.question)
        _df = _df.drop_duplicates()
        _df['project_id'] = _df.subject
        _df = self.remap_cat('subject', _mapping, _df)
        _df = _df.dropna().reset_index(drop=True)
        index = pd.MultiIndex.from_tuples(zip(_df.subject, _df.datetimeAnswered,
                                            _df.question), names=['subject', 'datetimeAnswered', 'question'])
        _df = pd.DataFrame(_df.answer.values, index=index,
                        columns=['answer']).unstack()
        _df.columns = _df.columns.droplevel()
        _df.columns = questions
        _df = _df.reset_index()
        _df['project'] = self.name
        _df['project'] = pd.Categorical(_df['project'])
        for col in _df.columns:
            _df[col] = pd.Categorical(_df[col])
        self.wellbeing = _df.copy()
        self.save_pickle()
        if verbose:
            print(f"Processing TIHM {self.name} wellbeing domain : {time.perf_counter()-tic:0.2f} seconds")   