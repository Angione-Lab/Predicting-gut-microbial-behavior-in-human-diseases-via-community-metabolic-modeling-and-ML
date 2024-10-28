<h1> Predicting Patient-specific Disease States </h1>
This repository contains the code and data to reproduce the results presented in the paper: Chaimaa Tarzi, Guido Zampieri, Suraj Verma, Stefano Campanaro, Neil Sullivan, Claudio Angione, "Predicting gut microbial behavior in human diseases via community metabolic modeling and machine learning". </br></br>

<p>This project focuses on developing a pipeline involving pairwise genome-scale metabolic modeling, statistical analysis, and machine learning to predict patient-specific disease states. Initially, genome-scale metabolic models are constructed and curated using CarveMe and ModelPolisher, followed by pairwise reconstruction to analyze species interactions within microbial communities. </br></br>
The interactions are evaluated under defined metabolic conditions, classifying them into positive (e.g., mutualism) or negative (e.g., competition) based on growth rate impacts. The data is further analyzed using R, where patient-level features are visualized through UMAP and statistical tests, highlighting significant differences in exchange reactions. </br></br>
A comprehensive machine learning framework is then employed, leveraging multiple classifiers including SVC, Random Forest, and Gradient Boosting to predict disease states. SHAP values are used to interpret the model, providing insights into feature importance and the underlying biological interactions. The ensemble model, combining the strengths of individual classifiers, enhances prediction accuracy and robustness.
</p>
<h1> Overview </h1> 
The Predicting Patient-specific Disease States project is structured into five main parts:
    <ul>
    <li><p>Genome-scale metabolic model construction: a Jupyter Notebook that guides users through the process of building and refining GEMs.</p></li>
    <li><p>Pairwise reconstruction and analysis: a pipeline for analyzing interactions between different metabolic models.</p></li>
    <li><p>Growth rate calculation: a method to evaluate the growth rates of microbial species under defined conditions.</p></li>
    <li><p>R analysis: statistical analysis and visualization of patient-level features using R.</p></li>
    <li><p>Machine learning: this part involves applying machine learning techniques to analyze patient data and predict disease states. By utilizing various classifiers and SHAP values, we can interpret model predictions and identify key features influencing patient outcomes. Additionally, ensemble modeling techniques are employed to enhance prediction accuracy, leveraging the strengths of multiple classifiers for more robust results.</p></li>
    </ul>
<h2>1. Genome-Scale Metabolic Model Construction</h2>
<p>
To get started, ensure you have the necessary input sequence files prepared and saved in the designated folders.
CarveMe serves as an automated reconstruction tool that generates a simulation-ready metabolic model based on MAG sequence data. </br> </br>
Begin by executing the command line in <strong>Tutorial GEM reconstruction.ipynb</strong> script, which will create the initial draft GEM. Following this, the use of the MEMOTE test suite is important to run a series of standardized tests that assess various aspects of the model, generating a report that highlights its strengths and areas for improvement. </br></br>
This automated testing process facilitates tracking incremental changes and ensures that the model meets community standards. Then to enhance and improve the quality of the model through improved annotation and curation. </br></br>
Note: Ensure that all prerequisite software and dependencies for CarveMe, MEMOTE, and ModelPolisher are installed prior to running the scripts. 
</p>
<h2>2. Pairwise Reconstruction and Analysis</h2>
<p>
In this step, users will place all polished genome-scale metabolic models in the <strong>single_model</strong> folder. The main script <strong>main.py</strong> will be used to run the pairwise modeling analysis.
Steps to Create Pairwise Models involves the widget 4  Metabolic models of 2-species communities created using COBRA Tools for the Python computational framework (COBRApy). This allows for the integration of metabolic models to analyze interactions within microbial communities.</br></br>
Under defined metabolic conditions in the <strong>data</strong> folder, the growth rates of each species are estimated both in isolation and in the presence of another species within the community in widget 5. This analysis helps to understand how species interact metabolically, and then the flux of each community exchange reaction is stored in the '_rawExRxns.csv'.</br></br>
The interactions occurring between pairs of species under the metabolic conditions defined in widget 5 are predicted. In widget 6 the interactions can be classified as; positive interactions: commensalism and mutualism, where at least one species benefits and no species suffers. Negative Interactions: parasitism, amensalism, and competition, where at least one species suffers. Neutralism: No interactions between the species. Interactions are determined based on their effect on the predicted growth rates compared to growth in isolation. An effect greater than 10% indicates an interaction, with the direction classified as negative or positive based on whether the species grows slower or faster in the community.
</p>
<h2>4. R Analysis</h2>
<p>
This section details the R analysis performed to visualize and analyze patient-level features.
Run R Code <strong>'Aggregated_Patient_level_analysis.R'</strong> in the R folder to aggregate the pairwise model to the patient-level feature. Following this, UMAP and violin plots are performed on a community level and patient level. The exchange reactions are read from the _rawExRxns.csv file, and UMAP is used to visualize the data. It is also presented in the manuscript and illustrated in Figures 3 and 4.</br></br>
In the next step statistical analysis of flux rates is performed the analysis includes pairwise t-tests to compare flux rates across different labels <strong> 'Community_reaction_enrichment.R'</strong>.
</p>
<h2>5. Machine learning</h2>
<p>
To prepare the data for machine learning analysis, we begin by organizing the samples into distinct categories. The <strong>'community_model_extract.ipynb'</strong> script is utilized to label the samples accordingly. Healthy samples are assigned the label "Healthy," while COVID samples are labeled "COVID." This is achieved by concatenating the two datasets into a single DataFrame, combined_samples. This step ensures that any columns with missing values are removed, resulting in a cleaner dataset for analysis. Next, we run the <strong>DataPreprocessing.py</strong> script, which includes a function extract_interaction_type(cohort, label). This function processes each sample in the cohort, reading the data from CSV files and calculating the interaction types and mean flux values.</br></br>
After preprocessing the data, we proceed to implement machine-learning models to predict patient-specific disease states in the folder 'ML model'. The first step involves training multiple classifiers and interpreting their results using SHAP values to understand feature importance. In the <strong>6ML_classifier</strong> script, we train several models, including Support Vector Classifier (SVC), K-Nearest Neighbors (KNN), Gradient Boosting (GB), Decision Tree, Random Forest, and Gradient Boosting classifiers. The SHAP values are calculated to interpret the model predictions. This code calculates the SHAP values for each model, visualizes the distribution of SHAP values for the top features, and categorizes them based on their importance.</br></br>
In the next step, run the <strong> Ensemble_model </strong> script to implement an ensemble model to improve prediction accuracy. The ensemble model combines the predictions of multiple classifiers using a voting mechanism. The following code snippet defines the ensemble model using SVC, KNN, and gradient-boosting classifiers.
</p>
<h2>Plots</h2>
Main plots and biological interpretation plots are plotted in R.

<h2>Project Structure</h2>
<img src="https://github.com/user-attachments/assets/092603ad-1a72-4316-b78a-2694c30bb448" class="inline" width="40%" style="display: block; margin: auto;" />

Note: Python 3.6.x is required, a check is specifically put into the code before it continues.</br>
Jupyter notebook server is required</br>
Ensure all pip dependencies are installed as listed in requirements.txt
    

