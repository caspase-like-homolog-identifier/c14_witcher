#!/usr/bin/env python
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
import pandas as pd
import argparse

c14reference = pd.read_csv("c14reference.tsv", delimiter = "\t")
c14reference.shape

c14_ref = c14reference.dropna()
c14_ref.shape

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



# +
n_entries, n_features  = X_train.shape

print("n_feat\tn_trees\tmae\tmse\trmse\taccuracy")

for n_feat  in range(1, n_features + 1):
      print()
      for n_trees in range(min_trees,  max_trees+trees_step,  trees_step):
          # y_pred = regressor.predict(X_test)
          # mae = metrics.mean_absolute_error(y_test, y_pred)
          # mse = metrics.mean_squared_error(y_test, y_pred)
          # rmse = np.sqrt(metrics.mean_squared_error(y_test, y_pred))
          # errors = abs(y_pred - y_test)
          # mape = 100 * (errors / y_test)
          # accuracy = 100 - np.mean(mape)
          # print("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(n_feat,n_trees, mae, mse, rmse, accuracy))
# -

             


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('c14reference', type=argparse.FileType('r')
    parser.add_argument('-n','--n_jobs', default = 4)
    parser.add_argument('-z','--test_size', default = 0.2)

    
    params = parser.parse_args()
    get_data(params)
