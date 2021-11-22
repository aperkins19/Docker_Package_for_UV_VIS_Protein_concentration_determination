import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from seaborn.matrix import heatmap

import Bio
from Bio import SeqIO, SearchIO, Entrez
from Bio import Seq
from Bio.SeqUtils import GC
from Bio.Blast import NCBIWWW
from Bio.Data import CodonTable

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import regex as re

# this script converts the a280 uv absorbance valves to concentrations for the 10 elutions for 8x purifications

well_list = ['A','B','C','D','E','F','G','H']

################################################# enter the protein names, e.g. MetRS, here in a list from well A to Well H
protein_list = ['RF2', 'RF1', 'EFTs', 'EFTu', 'EFG', 'IF3', 'IF2', 'IF1']

def beer_lambert_conversion(absorbance, path_length, extinction_coefficient, dilution_factor):

    """ Beer Lambert equation: A=l.c.ec
        Solve for c = A/(l.ec)

        Takes in raw absorbance, does the conversion, incorporates dilution correction"""


    # calculate concentration
    concentration = absorbance/(extinction_coefficient/path_length)

    #dilution correction
    concentration = concentration * dilution_factor

    return concentration

def process_dataset(collection, filename, crunched, wavelength_list, elution_list):

    rawdata = pd.read_csv(collection + filename, header= None)

    #how many wavelengths are there? 96 well plate dimentions assumed.
    number_of_wavelengths_tested = rawdata.shape[0]/9

    #check that the number given matches reality
    if int(number_of_wavelengths_tested) != len(wavelength_list):
        raise Exception('Number of wavelengths given does not equal number of wavelengths measured in ' + filename + ". Given: "+ str(len(wavelength_list))+ " vs Measured: "+ str(int(number_of_wavelengths_tested)))

    #divide up data file into discrete wavelength packets stored in a dictionary

    wavelength_dict = {}

    index_counter = 0
    base = index_counter

    for i in wavelength_list:
        index_counter += 9

        # make data subset its own object called data

        data = rawdata.iloc[base:index_counter,:]

        #add entry to dictionary and names it using the wavelength_list
        wavelength_dict[i] = data

        #resets base to index counter to move the indexing on for the next one.
        base = index_counter

    #print("wavelength dict")
    #print(wavelength_dict)

    # data trimming
    for i in wavelength_dict:

        # assign first column as rownames and delete
        wavelength_dict[i].index = wavelength_dict[i].iloc[:,0]
        del wavelength_dict[i][0]

        # assign elution list as colnames
        wavelength_dict[i].columns = elution_list

        #delete first row
        wavelength_dict[i] = wavelength_dict[i].drop('<>')


        # drop all rows not indexed as A, B or C
        wavelength_dict[i] = wavelength_dict[i].loc[well_list,:]

    #print('wavelength dict')
    #print(wavelength_dict)

####################################################################################################################################

    # convert a280s using beer Lambert
    #import the PURE system proteins sequences
    protein_sequences_df = pd.read_csv('proteinsequences.csv')

    ######## assign the relevant protein to it's well.


    #subset protein_sequences_df to only include those who include the protein_list
    # and tidy up
    selected = protein_sequences_df.set_index('Name').loc[protein_list].reset_index(inplace=False)
    selected = selected.drop('Recommended MWCO', axis=1)
    selected = selected.drop('Unnamed: 0', axis=1)

    # input wells into selected
    selected['Well'] = well_list

    # check if protein contains trypophan or phenylanaline

###############################################################################

    # loop through selected 'Sequence'
    # assign protein sequence to placeholder variable X
    # calculate extinction_coefficients
    # assign coefficients to lists

    reduced = []
    oxidised = []
    for p in selected['Sequence']:
        X = ProteinAnalysis(p)
        epsilon_prot = X.molar_extinction_coefficient()
        reduced.append(epsilon_prot[0])
        oxidised.append(epsilon_prot[1])

    selected['Reduced extinction coefficient M-1cm-1'] = reduced
    selected['Oxidised extinction coefficient M-1cm-1'] = oxidised

#####################################################################################

    ################################### convert absorbances into Molarity with Beer Lambert

    #path length assumptions.
    ### assuming volume of 100ul. Path length for 350ul is 1cm so 100/350 = 0.2857 cm
    # does not take into account miniscous at this time.
    path_length = 0.2857

    # copy the raw data
    raw_a280 = wavelength_dict['a280'].copy()

    moles_row_list = []

    # loop through the row
    for r in raw_a280.iterrows():
        #unpack this weird tuple thing
        well = r[0]
        absorbance_values = r[1]

        #get the blanks from the first row
        blank_avg = sum(raw_a280.iloc[0,[10,11]].values[[0,1]])/2
        #drop blanks
        absorbance_values = absorbance_values.drop(['blank1','blank2'])
        #blank subtract the rows
        absorbance_values = absorbance_values - blank_avg

        # get the relative extinct_coefficient
        respective_record = selected.loc[selected['Well']== well]
        reduced_extinct = respective_record['Reduced extinction coefficient M-1cm-1'].values[0]

        # convert with beer_lambert_conversion to get Moles per liter.
        moles_row = beer_lambert_conversion(absorbance_values, path_length, reduced_extinct,5)
        moles_row_list.append(moles_row)

    #make the list back into a dataframe
    protein_moles_df = pd.DataFrame(moles_row_list)

############################################################################################################

    ######### iterate through protein_moles_df rows and convert molarities to mgs/ml using individual molecular_weights

    mg_per_ml_list = []

    # loop through the row
    for m in protein_moles_df.iterrows():
        #unpack this weird tuple thing
        well = m[0]
        molarity_values = m[1]

        # get the relative molecular weight in KD
        respective_record = selected.loc[selected['Well']== well]
        MW_kd = respective_record['Protein MW kD'].values[0]

        # formula
        # M = g_per_ml / MW_kd
        # Solved for g_per_ml
        # g_per_ml = M * MW_kd

        g_per_ml = molarity_values * MW_kd
        mg_per_ml = g_per_ml * 1000

        mg_per_ml_list.append(mg_per_ml)

    #make the list back into a dataframe
    protein_mg_per_ml_df = pd.DataFrame(mg_per_ml_list)

    ###########################################################################################################

    # Heatmap

    plot_df = protein_mg_per_ml_df.copy()
    plot_df['Proteins'] = protein_list

    plot_df = plot_df.set_index(plot_df['Proteins'])
    plot_df = plot_df.drop('Proteins', axis=1)


    sns.heatmap(plot_df, annot=False, fmt="g", cmap='viridis', cbar_kws={'label': 'mg/ml'})

    plt.title('Elution Concentrations')
    plt.ylabel('Protein')
    plt.yticks(rotation=0)
    plt.xlabel('Elution Fraction')

    plt.savefig("elutions_mg_per_ml.png")

    ###################################################################################################################

    #clears the plot before starting the bar chart
    plt.clf()

    # in order to do the grouped barplot, I need to reformat the data to make the protein names and elution fraction numbers metadata in columns
    plot_df_transposed = plot_df.transpose()
    # grab the colnames and index names
    protein_columns = plot_df_transposed.columns
    elution_indexes = plot_df_transposed.index

    # initalise empty dataframe
    bar_plot_df = pd.DataFrame(columns=["Protein","Mg/Ml","Elution Fraction"])

    #iterate through the columns of the data frame
    for name, values in plot_df_transposed.iteritems():
        v = values
        # move the elution fractions from the index to a column
        v = v.reset_index()
        # rename the index column name
        v.columns = v.columns.str.replace('index', 'Elution Fraction')
        # assign protein
        v['Protein'] = str(name)
        v.columns = v.columns.str.replace(str(name), 'Mg/Ml')
        #append it on
        bar_plot_df = bar_plot_df.append(v, ignore_index = True)

    #print(bar_plot_df)

    sns.set_style('darkgrid')
    sns.set_palette('Dark2_r')

    sns.barplot(data=bar_plot_df, x="Elution Fraction", y="Mg/Ml", hue="Protein")
    plt.title('Elution Concentrations')
    plt.savefig("grouped_bar.png")

##################################################################################################################################

    ################################# subplots
    scaling_factor = 0.9
    # get the max elution value for the ylimit
    max_conc = bar_plot_df['Mg/Ml'].max()
    # add 10%
    max_conc = max_conc + (max_conc/10)

    # 2 rows, 4 columns, first plot
    fig, axes = plt.subplots(2, 4,figsize=(14*scaling_factor,7*scaling_factor))

    # grab the unique values from the protein columns
    protein_uniques = bar_plot_df.Protein.unique()

    # uses a row counter
    row_counter = 0
    for column_counter, p in enumerate(protein_uniques):
        # if index is greater than 4 then reset it to 0 and add one to the row counter
        if column_counter >= 4:
            column_counter = column_counter - 4
            row_counter = 1

        # from the big data frame, grab the rows where the Protein column equals the entry in protein_uniques
        individual_protein_slice = bar_plot_df[bar_plot_df['Protein']==p]

        # make a bar plot using that slice. use the axis counters to assign it to the correct subplot. set the title to the name
        sns.barplot(data=individual_protein_slice, x="Elution Fraction", y="Mg/Ml", ax=axes[row_counter,column_counter]).set(title=p)
        # set the ylim of that subplot using the max_conc caluclated earlier
        axes[row_counter,column_counter].set_ylim(0, max_conc)


    #sns.barplot(data=bar_plot_df, x="Elution Fraction", y="Mg/Ml", ax=axes[0,0])
    plt.title('Elution Concentrations')
    # make the space between the subplot rows better
    plt.subplots_adjust(hspace = 0.3)
    plt.tight_layout()
    plt.savefig("subplots_bar.png")

######################################################################################################################
    plot_df.to_csv('mg_per_ml.csv')



wavelength_list = ["a260", "a280", "a320", "a900", "a1000"]
elution_list = ["e1","e2","e3","e4","e5","e6","e7","e8","e9","e10","blank1", "blank2"]

crunched = []

os.getcwd()
collection = "data/"
for i, filename in enumerate(os.listdir(collection)):

    process_dataset(collection, filename, crunched, wavelength_list, elution_list)
