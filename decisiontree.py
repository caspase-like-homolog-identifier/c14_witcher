#!/usr/bin/env python
import pandas as pd
import pickle


class DecisionTree(object):

    def __init__(self, query_df,  pkl_model = "c14classifier.pickle"):
        
        with open(pkl_model, 'rb') as file_obj:
            self.c14classifier = pickle.load(file_obj)
            
        self.query_df = query_df
        self.query_data = query_df.values
        

        
    def classify(self):
        
        data_pred = self.c14classifier.predict(self.query_data)
        self.query_df["classification"] = data_pred

        return self.query_df





if __name__ ==  '__main__':
    query_df = pd.read_csv("real_data.tsv",
                           index_col=0,
                           header=None,
                           sep = "\t")
    dt = DecisionTree(query_df)
    df = dt.classify()
    print(df)
