import xgboost as xgb
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import LabelEncoder

metadata = pd.read_csv("../Capper_metadata_colors.csv")
metadata.head

metadata['Class'].value_counts()
metadata['Entity'].value_counts()

## count entities and filter to minimum number of cases
entity_counts = metadata['Entity'].value_counts()
entity_counts = entity_counts[entity_counts >= 21]
entity_keeps = entity_counts.index.tolist()

## filter metadata to minimum entities
metadata = metadata[metadata['Entity'].isin(entity_keeps)]
samp_keeps = metadata['Sentrix.ID...idat.'].value_counts().index.tolist()

beta_mat = pd.read_csv("../Capper_ssNoob_BetaValues.csv", index_col=0)
## Remove GSM* prefix from names and only leave IDAT
idats = [name.split('_', 1)[1] if '_' in name else name for name in 
beta_mat.columns]

beta_mat.columns = idats
beta_mat = beta_mat[samp_keeps]
beta_mat = beta_mat.transpose()
beta_mat.shape

## TODO summarize beta matrix to gene promoters using mean methylation

y_string = metadata['Entity'].tolist()
le = LabelEncoder()
y = le.fit_transform(y_string)

num_classes = len(entity_keeps)

# filter probes to most variable ones for modelling
probe_vars = beta_mat.var().sort_values(ascending = False)
top_probes = probe_vars[0:10000].index
beta_mat = beta_mat[ top_probes]

# Split data into training and test sets
X_train, X_test, Y_train, Y_test = train_test_split(beta_mat, y, 
test_size=0.3, random_state=54321)

# Specify the parameters for XGBoost
param = {
    'max_depth': 5,  # the maximum depth of each tree
    'eta': 0.1,  # the training step for each iteration
    'objective': 'multi:softprob',  # error evaluation for multiclass 
training
    'num_class': num_classes}  # the number of classes that exist in this 
dataset
num_round = 100  # the number of training iterations

# Transform data into DMatrix format for XGBoost
dtrain = xgb.DMatrix(X_train, label=Y_train)
dtest = xgb.DMatrix(X_test, label=Y_test)

# Train the XGBoost model
bst = xgb.train(param, dtrain, num_round)

# Save the model
bst.save_model('xgboost_model.json')

## Calculate Shapley values for feature importance
import shap

# Initialize the explainer
explainer = shap.Explainer(bst)

# Compute SHAP values 
shap_values = explainer.shap_values(X_test)

print(shap_values)
