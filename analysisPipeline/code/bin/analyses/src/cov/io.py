#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def parseDupInfo(filepath):
    colNames = ['id', 'idCount', 'dupId']
    df = pd.read_csv(filepath, sep='\t', header=None, names=colNames)
    df.drop(columns=['idCount'], inplace=True)

    ## Unnest the list of IDs into individual rows
    df['dupId'] = df['dupId'].str.split(',')
    df  = df.explode('dupId')

    ## Fix up some formatting
    df['dupId'] = df['dupId'].str.replace('^\\s+', '')
    df['dupId'] = df['dupId'].str.replace('\\s+$', '')
    df['dupId'] = df['dupId'].str.replace('\\s', '_')
    df['dupId'] = df['dupId'].str.replace('hCoV-19/', '')
    df['dupId'] = df['dupId'].str.replace('/[0-9]+\\|', '|')
    df.reset_index(drop=True, inplace=True)

    ## Update IDs without an EPI ID
    ## These are (typically) ones from EVAg
    cond = (~df['dupId'].str.contains('\\|'))
    df.loc[cond, 'dupId'] = df.loc[cond, 'dupId'] + '|' + df.loc[cond, 'dupId']

    f  = lambda x: ','.join(x)
    df = df.groupby(['id'])['dupId'].apply(f)
    df = df.reset_index()
    return df

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
