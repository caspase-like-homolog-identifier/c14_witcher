#!/usr/bin/env python
from sklearn.model_selection import train_test_split
from sklearn.tree import export_graphviz
from sklearn.tree import DecisionTreeClassifier
from IPython.display import Image  
from sklearn import metrics
from six import StringIO
import pandas as pd
import pydotplus
import argparse

c14reference = pd.read_csv("c14reference.tsv", delimiter = "\t")
c14reference.shape

c14_ref = c14reference.dropna()
c14_ref.shape

feature_cols = c14_ref.columns[:-1]

c14_ref = c14_ref[['p20', 'linker', 'p10','Classification']]

# +
#attributes 
y = c14_ref.loc[:,"Classification"].values

# #labels
X = c14_ref.drop(["Classification"], axis = 1).values
# -

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

classifier =  DecisionTreeClassifier(random_state=0)

classifier.fit(X_train, y_train)

y_pred = classifier.predict(X_test)

y_pred          


y_test

print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

dot_data = StringIO()
export_graphviz(classifier, 
                out_file=dot_data,  
                filled=True, 
                rounded=True,
                special_characters=True,
                feature_names = feature_cols)

graph = pydotplus.graph_from_dot_data(dot_data.getvalue())  
graph.write_png('c14classifier.png')
Image(graph.create_png())


