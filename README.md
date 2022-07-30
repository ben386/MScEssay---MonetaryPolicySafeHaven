# MScEssay---MonetaryPolicySafeHaven

"How Large is the Effect of Monetary Policy on Safe Haven Asset Prices?" by Benjamin O'Sullivan.

+------------------------------------------+

I estimate the time-varying effects of unanticipated monetary policy shocks on gold prices and a long-term US Treasury yield spread. I use a 
time-varying parameter vector autoregressive model. A high-frequency identification scheme is used to non-exactly identify the structural monetary policy
shocks. Both the model and the non-exact identification scheme are proposed by Paul (2020). Computing the relative impulse response functions of the 
structural model, I find that gold and the long-term US Treasury spread experience time-varying effects under a contractionary monetary policy shock. 
The long term Treasury spread further widens under a monetary tightening before and during the great recession, whereas gold has continually become 
more sensitive to monetary tightening over time.

+-------------------------------------------+

***IMPORTANT:*** I have used the replication code for "The Time-Varying Effect of Monetary Policy on Asset Prices" from [Pascal Paul](http://www.pascalpaul.de/replication-codes-varx/). I have made the necessary adjustments to carry out the estimation procedure I have laid out in my paper. 

**Citation:** *Paul, Pascal. 2020. “The Time-Varying Effect of Monetary Policy on Asset Prices”,
Review of Economics and Statistics, Vol. 102(4), pp. 690-704.*

**Instructions:**
1. Download the replication code
2. Under the "Time-Varying Parameter VAR" folder replace the *Run_TVP_VAR.m* and *figure_settings.m* scripts with the ones in this repo
3. Under the path "Data/Macro Time Series" add the *macrodata.xls* to the "Macro Time Series" folder
