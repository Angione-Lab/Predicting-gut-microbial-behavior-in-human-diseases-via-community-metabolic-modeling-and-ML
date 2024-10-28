
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


from sklearn.metrics import average_precision_score, roc_auc_score, accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict, cross_val_score
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB


from lightgbm import LGBMClassifier
import shap



from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

rc = {'figure.figsize':(10,5),
'axes.facecolor':'white',
'axes.grid' : True,
'grid.color': '.8',
'font.family':'Times New Roman',
'font.size' : 15}
plt.rcParams.update(rc)

#' GRASolo ', ' GRBSolo '

df = pd.read_csv('ExchangeReactions.csv', index_col = 0)
df.drop(['label'], axis = 1, inplace = True)


df.columns = [col.replace('modelA_', '') for col in df.columns]
df.columns = [col.replace('modelB_', '') for col in df.columns]


df_t = df.T
df_t['reactions'] = df_t.index

df = df_t.groupby(['reactions']).mean().T




flux = pd.read_csv("covid_healthy.csv", index_col=0)

flux['CommunityGrowthRate'] = flux[[' GRSpeciesAFull ', ' GRSpeciesBFull ']].mean(axis = 1)
flux['SoloGrowthRate'] = flux[[' GRASolo ', ' GRBSolo ']].mean(axis = 1)

flux = flux[[' Neutralism', ' Amensalism', ' Competition', 'CommunityGrowthRate', 'SoloGrowthRate', 'label']]



df = pd.merge(df, flux, left_on=df.index, right_on=flux.index)
df.index = df['key_0']
df.drop('key_0', axis=1, inplace=True)
# df.to_csv("aggregated_patient_level_features.csv")

df = pd.get_dummies(df, columns=['label'], drop_first=True)

x = df.drop(['label_H'], axis = 1)



y = df['label_H'].values


RANDOM_STATE = 1234

X_train, X_test, Y_train, Y_test = train_test_split(x,y, test_size=0.2,
                                     stratify=y,
                                    random_state=RANDOM_STATE,
                                    shuffle=True)

mm = StandardScaler(with_mean=False)
X_train = pd.DataFrame(mm.fit_transform(X_train), columns=X_train.columns)
X_test = pd.DataFrame(mm.transform(X_test), columns=X_test.columns)



import time

# Use GridSearchCV to find the best parameters:
def model_best_estimator(x_train, y_train, RANDOM_STATE=12, cv=5):
    
    # SVC 
    t0 = time.time()
    SVC_grid = {"kernel": ["sigmoid", 'poly'],
                "gamma" : ["scale", "auto"], 
                'C': [0.001, 0.01, 0.1, 1, 10, 100]}

    grid_log_reg = GridSearchCV(SVC(random_state=RANDOM_STATE, max_iter=1000),
                                SVC_grid, cv=cv, n_jobs=-1)
    grid_log_reg.fit(x_train, y_train)

    # get the SVC with the best parameters.
    svc_model = grid_log_reg.best_estimator_
    t1 = time.time()

    print("Best fit parameter for SVC", svc_model)
    print("Elapsed time {:.2f} s".format(t1 - t0))

    
    # KNN
    t2 = time.time()
    knears_params_grid = {"n_neighbors": list(range(5,8)), 
                          "metric": ('minkowski', 'euclidean')}
    
    grid_knears = GridSearchCV(KNeighborsClassifier(), knears_params_grid, cv=cv)
    grid_knears.fit(x_train, y_train)
   
    # KNN best estimator
    knn = grid_knears.best_estimator_
    t3 = time.time()
    print("\nBest fit parameter for KNN", knn)
    print("Effective metric knn:", knn.effective_metric_)
    print("Elapsed time {:.2f} s".format(t3 - t2))
    
    # LGBMClassifier
    t10 = time.time()
    nb_params_grid = {
                        "max_depth": list(range(4,12,1)),
                        "num_leaves": [4,6],
                        }
    
    grid_nb = GridSearchCV(LGBMClassifier(), nb_params_grid, cv=cv)
    grid_nb.fit(x_train, y_train)
   
    # NB best estimator
    nb = grid_nb.best_estimator_
    t11 = time.time()
    print("\nBest fit parameter for NB", nb)
    print("Elapsed time {:.2f} s".format(t11 - t10))


    # DecisionTree Classifier:
    t4 = time.time()
    tree_params_grid = {"criterion": ["gini", "entropy"], 
                        #"splitter":["best"],
                        #"min_samples_split":range(3,10), 
                        "max_depth": list(range(8,12,1)),
                        #"min_samples_leaf": [3,4,6]
                        }
    grid_tree = GridSearchCV(DecisionTreeClassifier(random_state=RANDOM_STATE),
                             tree_params_grid, cv=cv)
    grid_tree.fit(x_train, y_train)
    
    # tree best estimator
    tree_clf = grid_tree.best_estimator_
    t5 = time.time()
    
    print("\nBest fit parameter for Decision Tree:", tree_clf)
    print("Elapsed time {:.2f} s".format(t5 - t4))

    # Random Forest Classifier
    t6 = time.time()
    rf_params_grid = {"criterion": ["gini", "entropy"], "max_depth": list(range(2,10,1)), 
                "min_samples_leaf": [1,2,3,4,6],  "n_estimators":[4,6,8]}

    grid_rf = GridSearchCV(RandomForestClassifier(random_state=RANDOM_STATE), 
                           rf_params_grid, cv=cv)
    grid_rf.fit(x_train, y_train)

    # random forest best estimator
    rf = grid_rf.best_estimator_
    t7 = time.time()

    print("\nBest fit parameter for Random Forest:", rf)
    print("Elapsed time {:.2f} s".format(t7 - t6))


    # GBoost Classifier
    t8 = time.time()
    gb_para = {'learning_rate': [0.01, 0.1, 1], "n_estimators":[4,6,8]}

    grid_gb = GridSearchCV(GradientBoostingClassifier(random_state=RANDOM_STATE), 
                           gb_para, cv=cv)
    grid_gb.fit(x_train, y_train)

    # GB best estimator
    gb = grid_gb.best_estimator_
    t9 = time.time()

    print("\nBest fit parameter for Random Forest:", gb)
    print("Elapsed time {:.2f} s".format(t9 - t8))
    
    
    return [svc_model, knn, nb, tree_clf, rf, gb]   



# Get testing model results
def predict_model(classifier, x_test, y_test):
    y_pre = classifier.predict(x_test)
    print(classification_report(y_test, y_pre, labels=[1,0]))
    
    # Confusion Matrix
    print('Confusion matrix:', classifier)
    cf_matrix = confusion_matrix(y_test, y_pre, labels=[1,0])
    ax =sns.heatmap(cf_matrix, annot=True, fmt="d", cmap="Blues",
                xticklabels=['Healthy', 'Covid'],
                yticklabels=['Healthy', 'Covid'])
    ax.set(xlabel="Predicted outputs", ylabel = "Actual outputs")
    plt.show()
  
    return y_pre

def train_test(X_train, y_train, X_test, y_test, RANDOM_STATE=RANDOM_STATE, cv=5):
    # Find best parameter for model
    model_select_result = model_best_estimator(X_train, y_train, RANDOM_STATE=RANDOM_STATE)
    
    log_reg, knn, lgbm, tree_clf, rf, gb = model_select_result
    
    results =[]
    # Train and get result
    for classifier in model_select_result:
        print("\nPredict model:", classifier)
        #print("Train result:")
        #y_pre = predict_model(classifier, X_train, y_train)
        
    
        print("Testing result:")
        y_pre = predict_model(classifier, X_test, y_test)
        accuracy = accuracy_score(Y_test, y_pre)
        precision = precision_score(Y_test, y_pre, average='weighted')
        recall = recall_score(Y_test, y_pre, average='weighted')
        f1 = f1_score(Y_test, y_pre, average='weighted')
        print(accuracy, precision, recall, f1)
        results.append({"classifier": type(classifier).__name__, "accuracy": accuracy, "precision":precision, "recall":recall, "f1": f1})
    
    
    #can add feature importance
    return [log_reg, knn, lgbm, tree_clf, rf, gb, results]

#%%

svc_model, knn, lgbm, tree_clf, rf, gb, results = train_test(X_train, Y_train, X_test, Y_test)
print(pd.DataFrame(results))

#%%
def visualise_shap(shap_values, feature_names, file_name = "svm_shap"):
    
    mean_abs_shap_values = pd.DataFrame(shap_values.values, columns=feature_names).abs().mean().sort_values(ascending=False)
    
    # Get the top 20 features
    top_features = mean_abs_shap_values.head(20).index
    
    # Filter SHAP values for the top 20 features
    shap_values_top20_df = pd.DataFrame(shap_values.values, columns=X_train.columns)[top_features]
    
    # Melt the DataFrame for Seaborn
    shap_values_melted_top20 = pd.melt(shap_values_top20_df, var_name='Feature', value_name='SHAP Value')
    
    # Determine the average SHAP value for categorization
    mean_shap_value = shap_values_melted_top20['SHAP Value'].mean()
    
    # Categorize SHAP values as 'High' or 'Low'
    shap_values_melted_top20['Importance'] = np.where(shap_values_melted_top20['SHAP Value'] >= mean_shap_value, 'High', 'Low')
    
    # Set the plot size and style
    plt.figure(figsize=(30, 30))
    sns.set(style="white")
    
    # Create the swarm plot with hue based on importance
    sns.catplot(y="Feature", x="SHAP Value", data=shap_values_melted_top20, kind="violin", color=".9", inner=None)
    sns.swarmplot(y="Feature", x="SHAP Value", hue="Importance", data=shap_values_melted_top20, size=3 , palette={"High": "red", "Low": "blue"})
    # sns.despine(bottom = True, left = True)
    # Rotate feature names for better readability
    plt.xticks(rotation=0)
    plt.title('SHAP Value Distribution by Top 20 Features')
    # plt.legend(title='SHAP Importance', loc='lower right')
    plt.savefig(file_name +'.pdf', format='pdf')
    plt.show()


def shap_interpret(model, X_train, X_test, feature_names):
    explainer = shap.Explainer(model.predict, X_train, feature_names=feature_names)
    # Calculates the SHAP values - It takes some time
    shap_values = explainer(X_test)
    
    return shap_values

#%%
for model in [svc_model, knn, lgbm, tree_clf, rf, gb]:
    shap_val = shap_interpret(model, X_train, X_test, X_train.columns)
    visualise_shap(shap_val, X_train.columns, file_name = "Shap_images/"+ type(model).__name__)

tree_clf = DecisionTreeClassifier(random_state=RANDOM_STATE)
tree_clf.fit(X_train, Y_train)
predict_model(tree_clf,  X_test, Y_test)
shap_val = shap_interpret(tree_clf, X_train, X_test, X_train.columns)
shap.summary_plot(shap_val,  X_test, class_names=Y_test)

#%%
shap.plots.heatmap(shap_val, instance_order=shap_val.sum(1))

shap_df = pd.DataFrame(shap_val.values, columns=X_train.columns)
shap_df['target'] = Y_test

# Generate the SHAP heatmap
plt.figure(figsize=(12, 8))

# Define colors for the target classes
unique_classes = np.unique(Y_test)
colors = ["#F79B9B", "#9BB9F4"] #sns.color_palette("hsv", len(unique_classes))
class_colors = {cls: colors[i] for i, cls in enumerate(unique_classes)}

# Map target classes to colors
row_colors = shap_df['target'].map(class_colors)

# Generate the SHAP heatmap with custom row colors
plt.figure(figsize=(12, 8))
sns.clustermap(shap_df.drop(columns=['target']), 
               method='ward', 
               row_colors=row_colors, 
               row_cluster=True, 
               yticklabels = False,
               z_score=0, 
               cmap="coolwarm", 
               center=0)

plt.title('SHAP Heatmap clustered by Instances with Custom Row Colors')
plt.show()

#%%

df1 = df[[' Amensalism', ' Competition', ' Neutralism', 'label_H']]
# Generate the SHAP heatmap
df1.sort_values(by='label_H', inplace= True)
plt.figure(figsize=(12, 8))
#df1.sort_values(by=' Neutralism' , ascending=[False], inplace= True)

# Define colors for the target classes
unique_classes = np.unique(Y_test)
colors = ["#F79B9B", "#9BB9F4"] #sns.color_palette("hsv", len(unique_classes))
class_colors = {cls: colors[i] for i, cls in enumerate(unique_classes)}

# Map target classes to colors
row_colors = df1['label_H'].map(class_colors)

# Generate the SHAP heatmap with custom row colors
plt.figure(figsize=(12, 8))
sns.clustermap(df1.drop(columns=['label_H']), 
               method='ward', 
               row_colors=row_colors, 
               row_cluster=False, 
               yticklabels = False,
               #z_score=0, 
               cmap="coolwarm", 
               center=0)

plt.title('Interaction_type')
plt.show()

#%%
df1 = df[['EX_ca2_e', 'EX_cl_e', 'EX_co2_e', 'EX_cobalt2_e', 'EX_fe3_e',
       'EX_fe3pyovd_kt_e', 'EX_h2o_e', 'EX_h_e', 'EX_k_e', 'EX_mg2_e',
       'EX_mn2_e', 'EX_nh4_e', 'EX_o2_e', 'EX_pyovd_kt_e', 'EX_zn2_e',
       'label_H']]
# Generate the SHAP heatmap
plt.figure(figsize=(12, 8))
df1.sort_values(by='label_H' , ascending=[False], inplace= True)

# Define colors for the target classes
unique_classes = np.unique(Y_test)
colors = ["#F79B9B", "#9BB9F4"] #sns.color_palette("hsv", len(unique_classes))
class_colors = {cls: colors[i] for i, cls in enumerate(unique_classes)}

# Map target classes to colors
row_colors = df1['label_H'].map(class_colors)

# Generate the SHAP heatmap with custom row colors
plt.figure(figsize=(12, 8))
sns.clustermap(df1.drop(columns=['label_H']), 
               method='ward', 
               row_colors=row_colors, 
               row_cluster=False, 
               yticklabels = False,
               z_score=0, 
               cmap="coolwarm", 
               center=0)

plt.title('SHAP Heatmap clustered by Instances with Custom Row Colors')
plt.show()


