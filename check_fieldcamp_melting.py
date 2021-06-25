# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     formats: ipynb,py:percent
#     notebook_metadata_filter: all,-language_info
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.7.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Summer camp melting

# %% trusted=true
import marutils

# %% trusted=true
import pyproj


# %% trusted=true
mar_me = marutils.open_dataset('/flash/tedstona/MAR-v3.11.5-ERA5-10km/*.nc', transform_func=lambda ds: ds.ME)
crs = mar_me.rio.crs
p = pyproj.Proj(crs)
ll = (-47.235548, 66.98181)
camp_ss20_x, camp_ss20_y = p(*ll)
mar_me = mar_me.sel(x=camp_ss20_x, y=camp_ss20_y, method='nearest')

# %% trusted=true
# Forecast
mar_forecast = marutils.open_dataset('/flash/tedstona/MAR_forecast/MARv3.11-20km-daily-NCEP-NCARv1-2021.nc.1').ME
mar_forecast = mar_forecast.sel(x=camp_ss20_x, y=camp_ss20_y, method='nearest')

# %% trusted=true
avg = mar_me.groupby(mar_me['time.dayofyear']).mean(dim='time').cumsum()
store = {}
store['1980-2020'] = avg.to_pandas().squeeze()
for y in [2012, 2013, 2016, 2019, 2020]:
    d = mar_me.sel(time='%s' %y).cumsum()
    doy = d['time.dayofyear']
    as_pd = d.to_pandas().squeeze()
    as_pd.index = doy
    store[y] = as_pd


# %% trusted=true
# %matplotlib widget
import matplotlib.pyplot as plt
import pandas as pd
all_years = pd.DataFrame.from_dict(store)
all_years.plot()

f = mar_forecast.cumsum()
doy = f['time.dayofyear']
as_pd = f.to_pandas().squeeze()
as_pd.index = doy
as_pd.name = '2021'
as_pd.plot(linewidth=3)
plt.legend()

import datetime as dt
doy_now = int(dt.datetime.now().strftime('%j'))
plt.axvline(doy_now, linestyle=':', color='gray')
plt.xlim(140, 250)
plt.axvspan(199, 215, color='grey', alpha=0.2)

# NCEP-GFS divide - climato.be says that model is run with GFS from today -5d until today +7d.
plt.axvline(as_pd.index[-1] - 7 - 5, linestyle='--', color='tab:pink')
# Model run day:
plt.axvline(as_pd.index[-1] - 7 , linestyle='-', color='tab:pink')

plt.ylabel('Cumulative ME (mm w.e.)')
plt.title('{ll} / ({:d},{:d})'.format(int(camp_ss20_x), int(camp_ss20_y), ll=ll))
plt.xlabel('Day of Year')

# %% trusted=true
as_pd.index[-1]

# %%
