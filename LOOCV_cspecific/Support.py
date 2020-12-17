import pandas as pd
import xgboost as xgb

def runcv(splitter, dframe, goodcols, responsevar):
    dtrain=dframe.iloc[splitter[0]]
    dtest=dframe.iloc[splitter[1]]

    ddtrain=dtrain[goodcols]
    ddtest=dtest[goodcols]

    ctype=set(dtest["CancerType"])

    xgb1=xgb.XGBRegressor(max_depth=5, learning_rate=0.05, n_estimators=2000, silent=True, random_state=123)
    xgb1.fit(ddtrain, list(dtrain[responsevar]))
    return(ctype, dtest[responsevar], xgb1.predict(ddtest))