#!/usr/bin/env python3

"""Converts otu tables from diatom analysis into a suitable format for uploading to BIOSYS.

Additionally it will create a "otus_all.xlsx" for community analysis.
If extraction sheets are provided in an appropriate format then it will also calculate 
the most similar sample and provide the Baroni–Urbani–Buser coefficient similarity 
between the 2 samples."""

import sys
import os
import glob
import argparse
import pandas as pd
import send2trash
import re
import time
import multiprocessing as mp
from math import sqrt
import statistics as stat

version = ("\033[1;34m" + "\nBiosys Version: " + str(os.path.basename(__file__))[6:9] + "\033[0m")

known_issues = '''\033[93m This program is still in WIP stage \033[0m

current issues:
    -Multicore processing not supported.
    -Remove all instances of global values.
    -Baroni-Urbani-Buser similarity not at 100% functional due to:
        -Takes too long to perfrom for entire dataset. -solved?
        -Cut off for distance not determined.
        -Similarity doesnt take into account repeat samples.
'''

messages = [version,known_issues]

def version_info(messages):
    for message in messages:
        print(message)

#global values: (yes this is dumb but it works)

no_value = "no_value"

community_analysis_to_keep_list = ["all"]

control_regions = ["Blanks", "Positives", "Gblocks", "NTCs", "Unknowns", "TR"]

batch_num_dict = {  "Run_1":"25-01-19",
                    "Run_2":"15-01-19",
                    "Run_3":"18-01-19",
                    "Run_4":"21-01-19",
                    "Run_5":"18-01-19",
                    "Run_6":"24-01-19",
                    "Run_7":"28-01-19",
                    "Run_8":"31-01-19",
                    "Run_9":"04-02-19",
                    "Run_10":"07-02-19",
                    "Run_11":"12-02-19",
                    "Run_12":"15-02-19",
                    "Run_13":"23-02-19",
                    "Run_14":"28-03-19",
                    "Run_15":"19-04-19",
                    "Run_16":"23-04-19",
                    "Run_17":"17-05-19"
                    }

#Class Objects:

class FormatError(Exception):
    '''Formating of file is incompatible with this program.'''
    pass

class MathsError(Exception):
    '''A Maths or Logic error has occured in the code.'''
    pass

class Diatom_Sample:
    """A slice of an OTU table, and associated metadata for a diatom sample."""
    def __init__(self, sampleid, siteid, area, region, prn, sitename, sampledate, barcode, folder):
        if folder:
            try:
                self.folder = str(int(folder))
            except ValueError:
                self.folder = str(folder)
        else:
            self.folder = no_value

        if sampleid:
            try:
                self.sampleid = str(int(str(sampleid).replace(".","")))
            except ValueError:
                self.sampleid = str(sampleid)
        else:
            self.sampleid = "F" + str(self.folder)

        if siteid:
            self.siteid = str(int(str(siteid).replace(".","")))
        else:
            self.siteid = no_value

        if area:
            self.area = get_capital(str(area))
        else:
            self.area = no_value

        if region:
            self.region = get_capital(str(region))
            if region == "nan":
                print("Sample " + self.folder + " has no region, check input metadata file.")
                self.region = no_value
        else:
            self.region = no_value

        if self.region or self.area == no_value:
            pass
        else:
            self.reg_area = get_initials(region) + "-" + get_capital(area)

        if prn:
            try:
                self.prn = int(prn)
            except ValueError:
                self.prn = prn
        else:
            self.prn = no_value

        if sitename:
            try:
                regex = re.compile('[^a-zA-Z ]')
                self.sitename = regex.sub('', sitename)
            except AttributeError:
                self.sitename = no_value
            except TypeError:
                self.sitename = no_value
        else:
            self.sitename = no_value

        if sampledate:
            self.sampledate = sampledate
        else:
            self.sampledate = no_value

        if barcode:
            self.barcode = barcode
        else:
            self.barcode = no_value

        self.batch_num = no_value
        self.count = 0
        self.pass_fail = "Unsuccessful"
        self.analysis_date = no_value
        self.otu_tab = None
        self.sim = 0
        self.sim_sample = no_value
        self.note = ""
        self.plate_loc = no_value

    def set_analysis_date(self):
        '''Sets the date of analysis to the date of the MiSeq run.'''
        if self.batch_num == no_value:
            self.analysis_date = no_value
        else:
            try:
                date = batch_num_dict[self.batch_num]
            except KeyError:
                date = "Run metadata has not been set"
                print(date + " for sample: " + str(self.folder) + " " + str(self.batch_num))
            self.analysis_date = date

    def assign_results(self, otus, batch_num):
        '''Assigns otu table results to the sample.'''
        self.otu_tab = otus
        try:
            count = self.otu_tab[str(self.folder)].sum()
            self.count = count
        except KeyError:
            self.count = 0
            print("Seq count for " + str(self.folder) + " has been set to 0.")
        if self.count >= 3000:
            self.pass_fail = "Successful"
        if batch_num:
            self.batch_num = str(batch_num).split(".")[0]
            try:
                date = batch_num_dict[self.batch_num]
            except KeyError:
                date = "Run metadata has not been set"
                print(date + " for sample: " + str(self.folder) + " " + str(self.batch_num))
            self.analysis_date = date
        else:
            self.batch_num = no_value
            self.analysis_date = no_value

    def sort_control(self):
        if self.region == "Control":
            s_loc = str(self.folder).rfind("S")
            if s_loc == -1:
                self.folder = self.folder
            else:
                self.folder = str(self.folder)[0:s_loc]
            self.sampleid = self.folder + "_" + str(self.batch_num)
            if self.folder.lower()[0] == "b":
                self.region = "Blanks"
            elif self.folder.lower()[0] == "n":
                self.region = "NTCs"
            elif self.folder.lower()[0] == "p":
                self.region = "Positives"
            elif self.folder.lower()[0] == "g":
                self.region = "Gblocks"
            elif self.folder.lower()[0] == "t":
                self.region = "TR"
            else:
                self.region = "Unknowns"
                self.sampleid = "F" + str(self.folder)
            self.folder = self.folder.upper()

    def assign_surrounding_samples(self, sur_coords, row, col, sheet_name):
        '''Assigns the plate, coordinate, and surrounding coordinates for similaity analysis.'''
        if sur_coords:
            self.sur_samples = sur_coords
        if row and col:
            self.plate_loc = [row,col]
        if sheet_name:
            self.plate = str(sheet_name)

    def assign_most_sim_sample(self, sim, sample):
        '''Assigns the most similar surrounding sample and the similarity.'''
        self.sim = sim
        self.sim_sample = sample

    def amend_sample_note(self, note):
        if self.note == "":
            self.note = note
        else:
            self.note = self.note + "," + note

#Global functions:

def community_analysis_export(samples_otus, keep_list, control_regions):
    otus_all_writer = pd.ExcelWriter("otus_all.xlsx")
    print("Exporting all otus for community analysis")
    if keep_list[0] == "all":
        keep_list = ["Anglian", "Midlands", "South West", 
                     "Southern", "North West", "North East",
                     "Thames", "Unknowns", "Blanks",
                     "Positives", "Gblocks", "NTCs",
                     "TR", "Aberdeen", "Perth",
                     "Eurocentrl", "Dingwall", "Dumfries",
                     "Galashiels", "Bowlblank" ]

    for sample in samples_otus:
        if sample.count >= 1:
            if sample.region in keep_list:
                if sample.region in control_regions:
                    try:
                        sample.otu_tab.columns = ["PrefTaxon", sample.sampleid]   
                        interim_df = sample.otu_tab
                        main_df = pd.merge(main_df, interim_df, on="PrefTaxon")
                    except NameError:
                        main_df = sample.otu_tab
                else:
                    try:
                        interim_df = sample.otu_tab
                        main_df = pd.merge(main_df, interim_df, on="PrefTaxon")
                    except NameError:
                        main_df = sample.otu_tab
    df = format_df(main_df)
    df = df.transpose()
    df.to_excel(otus_all_writer, sheet_name="comunity_analysis", header=True, index=True)
    otus_all_writer.save()


def add_biosys_headers(df, biosys_dict):
    biosys_df = pd.DataFrame.from_dict(biosys_dict, orient='columns')
    biosys_df = biosys_df.rename(index={0:'siteid',1:'sampleid'})
    df = pd.concat([biosys_df, df], sort=True)
    return df
    
def format_df(df):
    df = df[~df.PrefTaxon.str.contains("batch_num")]
    df = df.set_index(list(df)[0])
    df = df.loc[(df!=0).any(axis=1)]
    return df

def filter_otus_by_region(region, samples_otus, writer, control_regions):
    print(region)
    if region == no_value:
        no_reg = open("samples_with_no_region_values.text", "w")
        no_reg.write("Region for the below samples is " + region + "\n")
        no_reg.write("Folder_id\tCounts\tSample_id\tSite_id\tPRN\n")
        for sample in metadata_list:
            if sample.region == region:
                no_reg.write(i.folder + "\t" + str(i.count) + "\t" + i.sampleid + "\t" + i.siteid + "\t" + i.prn + "\n")
        no_reg.close()
    elif region == "TR":
        biosys_siteid_dict = {}
        for sample_tr in samples_otus:
            if sample_tr.region == "TR":
                if sample_tr.count >= 3000:
                    if sample_tr.folder[2] == "3":
                        sample_original_fn = sample_tr.folder[2:8]
                    elif sample_tr.folder[2] == "4":
                        sample_original_fn = sample_tr.folder[2:9]
                    else:
                        sample_orignal_fn = sample_tr.folder
                    try:
                        interim_df = sample_tr.otu_tab
                        main_df = pd.merge(main_df, interim_df, on="PrefTaxon")
                    except NameError:
                        main_df = sample_tr.otu_tab
                    for sample_og in samples_otus:
                        if sample_og.folder == sample_original_fn:
                            if sample_og.region != region:
                                if sample_og.count > 1:
                                    interim_df = sample_og.otu_tab
                                    main_df = pd.merge(main_df, interim_df, on="PrefTaxon")
        try:
            df = format_df(main_df)
            df = add_biosys_headers(df, biosys_siteid_dict)
            df.to_excel(writer, sheet_name=region, header=True, index=True)
        except UnboundLocalError:
            print("    Region " + region + " had no passing samples.")
    elif region in control_regions:
        for sample in samples_otus:
            if sample.region == region:
                if sample.count >= 3000:
                    try:
                        interim_df = sample.otu_tab
                        main_df = pd.merge(main_df, interim_df, on="PrefTaxon")
                    except NameError:
                        main_df = sample.otu_tab
        try:
            df = format_df(main_df)            
            df.to_excel(writer, sheet_name=region, header=True, index=True)
        except UnboundLocalError:
            print("    Region " + region + " had no passing samples.")
    else:
        biosys_siteid_dict = {}
        for sample in samples_otus:
            if sample.region == region:
                if sample.count >= 3000:
                    biosys_siteid_dict[sample.folder] = [str(sample.siteid), str(sample.sampleid)]
                    try:
                        interim_df = sample.otu_tab
                        main_df = pd.merge(main_df, interim_df, on="PrefTaxon")
                    except NameError:
                        main_df = sample.otu_tab
        try:
            df = format_df(main_df)            
            df = add_biosys_headers(df, biosys_siteid_dict)
            df.to_excel(writer, sheet_name=region, header=True, index=True)
        except UnboundLocalError:
            print("    Region " + region + " had no passing samples.")

def get_region_list(sample_list):
    regions_all = []
    for sample in sample_list:
        regions_all.append(sample.region)
    region_unique = []
    for region in regions_all:
        if region not in region_unique:
            region_unique.append(region)
    region_unique = [x for x in region_unique if str(x) != 'nan']
    return region_unique

def delete_file(file_in):
    file_exists = os.path.isfile(file_in)
    if file_exists == True:
        send2trash.send2trash(file_in)

def save_sample_info(sample_list, writer):
    interim = "inter.text"
    delete_file(interim)
    file_interim = open(interim, "w")
    header = ("SampleID,SiteID,Region,Area,PRN,Site Name,Sample Date,Barcode,Date of analysis,"
              "BatchNum,Sequence Counts,Pass/Fail,FolderNumber,Notes\n")
    file_interim.write(header)
    for sample in sample_list:
        if sample.region == "Control":
            pass
        else:
            line = ( str(sample.sampleid) + "," + str(sample.siteid) + "," + str(sample.region) + "," +
                    str(sample.area) + "," + str(sample.prn) + "," + str(sample.sitename) + "," +
                    str(sample.sampledate) + "," + str(sample.barcode) + "," + 
                    str(sample.analysis_date) + "," + str(sample.batch_num) + "," +
                    str(sample.count) + "," + str(sample.pass_fail) + "," +
                    str(sample.folder) + "," + str(sample.note) + "\n" )
            file_interim.write(line)
    file_interim.close()
    sample_info_df = pd.read_csv(interim, delimiter = ",")
    sample_info_df.to_excel(writer, sheet_name="Sample batch information", header=True, index=False)

def baroni_urbani_buser_coefficient(df):
    a = 0 #Num species present in both samples
    b = 0 #Num species present in only sample 0
    c = 0 #Num species present in only sample 1
    d = 0 #Num species absent in both samples
    df = df.drop(columns="PrefTaxon")
    for row_num in range(0, df.shape[0]):
        #col 0 is sample 1, col 1 is sample 2
        #Create binary data
        count_0 = df.iloc[row_num][0]
        if count_0 > 1:
            count_0 = 1
        count_1 = df.iloc[row_num][1]
        if count_1 > 1:
            count_1 = 1
        #compare binary data
        if count_0 == count_1:
            if count_0 == 1:
                a=a+1
            else:
                d=d+1
        elif count_0 > count_1:
            b=b+1
        elif count_0 < count_1:
            c=c+1
        else:
            MathsError
    #calculate coefficient
    bub = float((sqrt(a*d) + a) / (sqrt(a*d) + a + b + c))
    return bub

#dif
dif = {"PrefTaxon":["EC","IE"],"A":[1,0],"B":[0,1]}
same = {"PrefTaxon":["EC","IE"],"A":[1,0],"B":[1,0]}
assert baroni_urbani_buser_coefficient(pd.DataFrame(data=dif)) == 0
assert baroni_urbani_buser_coefficient(pd.DataFrame(data=same)) == 1

def perform_similarity_checks(sample_list, writer):
    print("Performing Baroni–Urbani–Buser coefficient similarity checks...")
    sim_all = []
    for sample_cent in sample_list:
        status = sample_cent.folder
        count = sample_list.index(sample_cent)
        filled_len = int(round(60 * count / float(len(sample_list))))
        percents = round(100.0 * count / float(len(sample_list)), 5)
        bar = '=' * filled_len + '-' * (60 - filled_len)
        sys.stdout.write('[%s] %s%s Analysing sample: %s\r' % (bar, percents, '%', status))
        sys.stdout.flush()
        cent_df = sample_cent.otu_tab
        if sample_cent.pass_fail != "Successful":
            continue
        highest_sim = 0.0
        most_sim_sample = ""
        for sample_sur in sample_list:
            if sample_cent == sample_sur:
                pass
            else:
                if sample_sur.pass_fail != "Successful":
                    continue
                try:
                    if sample_cent.plate != sample_sur.plate:
                        pass
                    else:
                        if sample_sur.plate_loc in sample_cent.sur_samples:
                            sur_df = sample_sur.otu_tab
                            if not isinstance(sur_df, pd.DataFrame):
                                continue
                            df = pd.merge(cent_df, sur_df, on="PrefTaxon")
                            similarity = baroni_urbani_buser_coefficient(df)
                            sim_all.append(similarity)
                            if similarity > highest_sim:
                                highest_sim = similarity
                                most_sim_sample = sample_sur.folder
                        else:
                            sur_df = sample_sur.otu_tab
                            if not isinstance(sur_df, pd.DataFrame):
                                continue
                            df = pd.merge(cent_df, sur_df, on="PrefTaxon")
                            similarity = baroni_urbani_buser_coefficient(df)
                            sim_all.append(similarity)
                except AttributeError:
                    pass
        sample_cent.assign_most_sim_sample(highest_sim, most_sim_sample)
    print("\n")
    folder_nums = []
    sim_samps = []
    bub_cos = []
    plates = []
    plate_locs = []
    for sample in sample_list:
        folder_nums.append(sample.folder)
        sim_samps.append(sample.sim_sample)
        bub_cos.append(sample.sim)
        plate_locs.append(sample.plate_loc)
        try:
            plates.append(sample.plate)
        except AttributeError:
            plates.append(" ")
    print("Mean similarity of all samples: " + str(stat.mean(sim_all)))
    print("Standard deviationm of similarity of all samples: " + str(stat.stdev(sim_all)))
    sim_dict = {"FolderNumber":folder_nums,"Most Similar Sample":sim_samps, "BUB_Co":bub_cos, "Plate":plates}
    sim_df = pd.DataFrame.from_dict(sim_dict)
    sim_df.to_excel(writer, sheet_name="BUB_coefficient", header=True, index=False)

def get_surrounding_coords(row, col):
    if row == 0:
        surround_rows = [row, row+1]
    elif row == 7:
        surround_rows = [row-1 ,row]
    else:
        surround_rows = [row-1,row,row+1]

    if col == 1:
        surround_cols = [col,col+1]
    elif col == 12:
        surround_cols = [col-1,col]
    else:
        surround_cols = [col-1,col,col+1]
    sur_coords = []

    for sur_row in surround_rows:
        for sur_col in surround_cols:
            sur_coords.append([sur_row,sur_col])
    return sur_coords

def import_extraction_sheets(data_dir, xl_file, samples):
    print("Importing extraction sheets...")
    dir_abspath = os.path.abspath(data_dir)
    xl_abspath = os.path.join(dir_abspath, xl_file)
    xl = pd.ExcelFile(xl_abspath)
    sheet_names = xl.sheet_names
    sample_groups = {}
    for sheet_name in sheet_names:
        sheet = xl.parse(sheet_name=sheet_name)
        if sheet.empty:
            pass
        else:
            sheet = sheet.set_index(list(sheet.columns.values)[0])
            print(sheet.head())
            for row in range(0,sheet.shape[0]):
                for col in range(1,sheet.shape[1]+1):
                    try:
                        samplesheet_id = str(int(sheet.iloc[row][col])).upper().strip()
                    except ValueError:
                        samplesheet_id = str(sheet.iloc[row][col]).upper().strip()
                    sur_coords = get_surrounding_coords(row,col)
                    for sample in samples:
                        if str(sample.folder) == samplesheet_id:
                            try:
                                if sample.plate:
                                    sample.amend_sample_note("Sample also found on " + sheet_name)
                            except AttributeError:
                                sample.assign_surrounding_samples(sur_coords, row, col, sheet_name)
                  

def import_otus(file_name):
    otus = pd.read_csv(file_name, delimiter = "\t")
    otus = otus.rename(index=str)
    return otus    

def import_otu_tables_main(directory, sample_list):
    dir_abspath = os.path.abspath(directory)
    file_paths = glob.glob(str(dir_abspath) + "/*.tsv")
    files = []
    for file_path in file_paths:
        path,file_name = os.path.split(file_path)
        files.append(file_name)
    for file_name in files:
        otus = import_otus(path + "/" + file_name)
        headers = list(otus.columns.values)
        headers.pop(0)
        new_headers = {} #changes made to accomadate formatting of re-demux
        for header in headers:
            try:
                header_split = header.split(".")
                if header_split[0][0].lower() in ["b","n","p","g","t","u"]:
                    new_header = header_split[0] + header_split[1]
                    new_header = str(new_header)
                else:
                    new_header = header_split[0]
                new_headers[header] = new_header
            except IndexError:
                print("Header " + str(header) + " has been assumed to have been changed manually.")
                new_headers[header] = header
        otus = otus.rename(columns=new_headers)
        headers = list(otus.columns.values)
        headers.pop(0)
        headers_added = []
        for header in headers:
            for sample in sample_list:
                if header == str(sample.folder):
                    if sample.count > 3000:
                        headers_added.append(header)
                    else:   
                        df = otus[["PrefTaxon"]].copy()
                        df[header] = otus[header].copy()
                        sample.assign_results(df, file_name)
                        headers_added.append(header)
        for header in headers:
            if header in headers_added:
                pass
            else:
                sample = Diatom_Sample(None, None, None, "Control", None, None, None, None, header)
                df = otus[["PrefTaxon"]].copy()
                df[header] = otus[header].copy()
                sample.assign_results(df, file_name)
                sample.sort_control()
                sample.otu_tab.columns = ["PrefTaxon", sample.sampleid]
                sample_list.append(sample)
    return sample_list

def get_initials(string):
    xs = (string)
    words_list = xs.split()
    initials = ""
    for word in words_list:
        initials = initials + word[0].upper()
    return initials

def get_capital(string):
    xs = (string)
    word_list = xs.split()
    words_num = len(word_list)
    words_return = ""
    for word in word_list:
        word_return = word[0].upper() + word[1:].lower()
        words_return = words_return + word_return + " "
    return words_return.strip() 

def import_metadata_sepa(inputxl, directory):
    dir_abspath = os.path.abspath(directory)
    xl_abspath = os.path.join(directory, inputxl)
    print("Metadata input file: " + xl_abspath)
    xl = pd.ExcelFile(xl_abspath)
    samples_reg = xl.parse(sheet_name=0)
    samples_list = []
    samples_info = samples_reg[["Region", "S_SAMPLING_PT_DESC", "SAMPLE_NUMBER", "SAMPLED_DATE"]]
    samples_info.columns = ["region", "site_name", "sepa_num", "sample_date"] 
    for row in samples_info.itertuples():
        try:
            sepa_num = str(int(row[3]))
        except ValueError:
            sepa_num = str(row[3])
        region = get_capital(str(row[1]))
        site_name = str(row[2])
        try:
            sample_date = str(row[4])
        except ValueError:
            print("Value for SampleDate not accepted for sample: " + sepa_num + " orignal value: " + str(row[10]))
            SampleDate = "01/01/18"
        diatom_sample_info = Diatom_Sample(sepa_num, "1", "area", region, sepa_num, site_name, sample_date, sepa_num, sepa_num)
        samples_list.append(diatom_sample_info)
    return samples_list

def import_metadata_ea(inputxl, directory):
    dir_abspath = os.path.abspath(directory)
    xl_abspath = os.path.join(directory, inputxl)
    print("Metadata input file: " + xl_abspath)
    xl = pd.ExcelFile(xl_abspath)
    samples_reg = xl.parse(sheet_name=0)
    samples_list = []
    samples_info = samples_reg[["Region", "Area", "BIOSYS site ID", "Water body", "Site/Station Name", "Sample Id", "Barcode received", "PRN", "Folder", "Sample Date"]]
    samples_info.columns = ["Region", "Area", "SiteID", "WaterBody", "SiteName", "SampleId", "Barcode_R", "PRN", "Folder", "SampleDate"] 
    for row in samples_info.itertuples():
        region = row[1]
        area = row[2]
        siteid = row[3]
        site_name = row[5]
        sampleid = row[6]
        barcode = row[7]
        prn = row[8]
        folder = row[9]
        sample_date = row[10]
        diatom_sample_info = Diatom_Sample(sampleid, siteid, area, region, prn, site_name, sample_date, barcode, folder)
        samples_list.append(diatom_sample_info)
    return samples_list

def get_args():
    parser = argparse.ArgumentParser(description="Processes diatom data into regions.")
    parser.add_argument("--input_xl", help="Semi-Optional: input excel file in .xlsx format containing information about the samples.", default="infoEA.xlsx", required=False)
    parser.add_argument("--input_dir", help="Semi-Optional: input directory for all input files files.", default="Data", required=False)
    parser.add_argument("--area", help="Semi-Optional: EA or SEPA, defaults to EA.", default="EA", required=False)
    parser.add_argument("--plate_info", help="Semi-Optional: gives name of extraction plate file, defaults to Plates.xlsx.", default="Plates.xlsx", required=False)
    parser.add_argument("--similarity", help="Semi-Optional: If True perfroms similarity checks.", default=False, required=False)
    args = parser.parse_args()
    return args

def main():

    options = get_args()
    
    writer = pd.ExcelWriter("output.xlsx")

    print("Area: " + str(options.area))
    if str(options.area).upper() == "SEPA":
        samples_meta_list = import_metadata_sepa(options.input_xl, options.input_dir)
    else:
        samples_meta_list = import_metadata_ea(options.input_xl, options.input_dir)
    
    samples_otus = import_otu_tables_main(options.input_dir, samples_meta_list)

    if options.similarity != False:
        try:
            import_extraction_sheets(options.input_dir, options.plate_info, samples_otus)
            perform_similarity_checks(samples_otus, writer)
        except FileNotFoundError:
            print("Extraction sheets could not be found, skipping similarity checks.")
    else:
        print("Not perfroming similarity checks, set --similarity True to perform this test.")
    
    save_sample_info(samples_otus, writer)

    region_list = get_region_list(samples_otus)
    
    for region in region_list: 
        filter_otus_by_region(region, samples_otus, writer, control_regions)

    community_analysis_export(samples_otus, community_analysis_to_keep_list, control_regions)

    delete_file("inter.text")

    writer.save()

if __name__ == "__main__":
    version_info(messages)
    main()


