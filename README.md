# HARModel

An implementation of the Heterogeneous AutoRegressive model from Corsi(2009)

If you find any bugs or think improvements can be made, please contact me at:

emilsjoerup@live.dk

TODO:
Implement HARJ , HARQ, FULL-HARQ, and Semivariance HAR.

Nicer show and print methods

Add compatibility for using a secondary vector for rolling forecasting

Known bug(s):

1: In HARForecast, iNRoll cannot be higher than half the length of the input vector, which is not correctly captured in the documentation nor code - will be made more clear in the next update. (is fixed - atleast in code)

2: Pressing the expand button on the Observations and vForecastComp vectors in a HARforecast object provides an error I can't seem to figure out. Luckily the vectors are accessible through `HARforecast@Data[["Observations"]]` and `HARforecast@Data[["ForecastComparison"]]` respectively.
