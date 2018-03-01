# Portfolio Models
A small project tracking the performance of three portfolios. A more detailed report is upcoming which will summarize the projects results. Here is a summary:

Highlights and Details:
- Fama-French 3 Factor model used to estimate asset returns
- Portfolio tracks up to 20 stocks (predetermined and fixed) from the S&P500, representing all sectors
- Real adjusted stock prices from 2012 to 2015 used. 
- Investment periods last 6 months, starting in 2013 and ending in 2015 for a total of 6 periods. Portfolio rebalancing and adjustments occur at the end of each period.
- Transaction costs are based on 0.5% of trading volume. These costs are not deducted from the portfolio value rather are tracked in a seperate account.
- Three portfolio methods are: MVO solved from standard constrained convex optimization; MVO with cardinality constraints, solved with mixed integer programming techniques offered by Gurobi; Black Litterman model with my own adaptive views.
- Initial investment of $100.
