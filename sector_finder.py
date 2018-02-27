import urllib3
import pytz
import csv
import pandas as pd
import pandas_datareader.data as web

from alpha_vantage.timeseries import TimeSeries
from bs4 import BeautifulSoup
from datetime import datetime

SITE = "http://en.wikipedia.org/wiki/List_of_S%26P_500_companies"
START = datetime(1900, 1, 1, 0, 0, 0, 0, pytz.utc)
END = datetime.today().utcnow()


def scrape_list(site):
    http = urllib3.PoolManager()
    response = http.request('GET', SITE)
    soup = BeautifulSoup(response.data, "lxml")

    table = soup.find('table', {'class': 'wikitable sortable'})
    sector_tickers = dict()
    for row in table.findAll('tr'):
        col = row.findAll('td')
        if len(col) > 0:
            sector = str(col[3].string.strip()).lower().replace(' ', '_')
            ticker = str(col[0].string.strip())
            if sector not in sector_tickers:
                sector_tickers[sector] = list()
            sector_tickers[sector].append(ticker)
    return sector_tickers

if __name__ == '__main__':
    sector_tickers = scrape_list(SITE)
    target_tickers = ['F', 'CAT', 'DIS', 'MCD', 'KO', 'PEP', 'WMT', 'C', 'WFC', 'JPM', 'AAPL', 'IBM', 'PFE', 'JNJ', 'XOM', 'MRO', 'ED', 'T', 'VZ', 'NEM']
    sectors = {}
    
    for i in range(0, len(target_tickers)):
        for sector in sector_tickers.keys():
            if target_tickers[i] in sector_tickers[sector]:
                if sector in sectors:
                    sectors[sector].append(target_tickers[i])
                else:
                    sectors[sector]= [target_tickers[i]]

    print(sectors)



    

    
