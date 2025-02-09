
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier

from sklearn.ensemble import VotingClassifier
import shap

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

flux = pd.get_dummies(flux, columns=['label'], drop_first=True)

df = pd.merge(df, flux, left_on=df.index, right_on=flux.index)
df.index = df['key_0']
df.drop('key_0', axis=1, inplace=True)


x = df.drop(['label_H'], axis = 1)



y = df['label_H'].values
RANDOM_STATE = 123

X_train, X_test, Y_train, Y_test = train_test_split(x,y, test_size=0.2,
                                     stratify=y,
                                    random_state=RANDOM_STATE,
                                    shuffle=True)

mm = StandardScaler(with_mean=False)
X_train = pd.DataFrame(mm.fit_transform(X_train), columns=X_train.columns)
X_test = pd.DataFrame(mm.transform(X_test), columns=X_test.columns)




#%%

# Use GridSearchCV to find the best parameters:
def model_best_estimator(x_train, y_train, RANDOM_STATE=12, cv=5):
    
    # define models
    svc = SVC(random_state=RANDOM_STATE, max_iter=3000)
    knn = KNeighborsClassifier()
    gb = GradientBoostingClassifier(random_state=RANDOM_STATE)
    
    # Define parameters
    params = {"svc__kernel": ["sigmoid", 'poly'],
                "svc__gamma" : ["scale", "auto"], 
                'svc__C': [0.001, 0.01, 0.1, 1],
                "knn__n_neighbors": list(range(5,7)), 
                "knn__metric": ('minkowski', 'euclidean'),
                "gb__learning_rate": [0.01, 0.1, 1], 
                "gb__n_estimators":[4,6,8] 
                }
       
    
    ensemble = VotingClassifier(
        estimators=[('svc', svc), 
                    ('knn', knn), 
                    ('gb', gb)
                    ],
        voting='hard')
    grid_ensemble = GridSearchCV(estimator=ensemble, param_grid=params, cv=cv)
    grid_ensemble.fit(x_train, y_train)

    # grid_ensemble.best_estimator_
    print(grid_ensemble.best_params_)
    
    return grid_ensemble  



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


#%%

from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score


ensemble_model = model_best_estimator(X_train, Y_train, RANDOM_STATE=RANDOM_STATE)

best_model = ensemble_model.best_estimator_

y_pre = predict_model(best_model, X_test, Y_test)

accuracy = accuracy_score(Y_test, y_pre)
precision = precision_score(Y_test, y_pre, average='weighted')
recall = recall_score(Y_test, y_pre, average='weighted')
f1 = f1_score(Y_test, y_pre, average='weighted')
print(accuracy, precision, recall, f1)


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
    plt.savefig(file_name +'_merged.pdf', format='pdf')
    plt.show()


def shap_interpret(model, X_train, X_test, feature_names):
    explainer = shap.Explainer(model.predict, X_train, feature_names=feature_names)
    # Calculates the SHAP values - It takes some time
    shap_values = explainer(X_test)
    
    return shap_values

#%%

shap_val = shap_interpret(ensemble_model, X_train,  X_train.loc[:19,:], X_train.columns)
visualise_shap(shap_val, X_train.columns, file_name = "Ensemble_model")

#%%
shap.plots.heatmap(shap_val, instance_order=shap_val.sum(1))

shap_df = pd.DataFrame(shap_val.values, columns=X_train.columns)
mean_abs_shap_values = shap_df.abs().mean().sort_values(ascending=False)

# Get the top 20 features
top_features = mean_abs_shap_values.head(20).index

# Filter SHAP values for the top 20 features
shap_values_top20_df = shap_df[top_features]

# shap_df = pd.DataFrame(shap_val.values, columns=X_train.columns)
shap_values_top20_df['Label'] = Y_test

# Generate the SHAP heatmap
plt.figure(figsize=(12, 8))

# Define colors for the target classes
unique_classes = np.unique(Y_test)
label_colors = ["#F79B9B", "#9BB9F4"] #sns.color_palette("hsv", len(unique_classes))
class_colors = {cls: label_colors[i] for i, cls in enumerate(unique_classes)}

# Map target classes to colors
row_colors = shap_values_top20_df['Label'].map(class_colors)

# Generate the SHAP heatmap with custom row colors
plt.figure(figsize=(12, 8))
g = sns.clustermap(shap_values_top20_df.drop(columns=['Label']).T, 
               method='ward', 
               col_colors=row_colors, 
               col_cluster=True, 
               row_cluster=False, 
               xticklabels = False,
               z_score=0, 
               cmap="PRGn", 
               center=0)
g.ax_heatmap.yaxis.set_ticks_position("left")
g.savefig('Ensemble_SHAP_heatmap.pdf', format='pdf')
# g.ax.title('SHAP Heatmap')
g.plt.show()


